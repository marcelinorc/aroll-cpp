/* 
 * File:   DegradedUnrollingKit.cpp
 * Author: elmarce
 * 
 * Created on March 7, 2016, 11:36 AM
 */

#include "opto/degradedUnrollingKit.hpp"
#include "opto/graphKit.hpp"
#include "libadt/dict.hpp"
#include "utilities/bitMap.hpp"

//Performs the DFS staring in the node given as parameter
//start:  Node from which the DFS will start
//pre: Current preorder counter in the DFS

//TODO: FIND HOW TO USE THE INDEX OF THE NODES
#define CAN_VISIT(nn) low[nn->_idx] == 0 //If low == 0 is not visited yet

//******************************************************************************
//------------------------------------------------------------------------------
//------------- APPROXIMATED UNROLLING -----------------------------------------
//Finds two paired
//Indicates if the loop is a signal loop.
//This happens when inside the loop values are assigned to an array
//being indexed by a function of the loop controlling variable

Node** DegradedUnrollingKit::find_array_storages(IdealLoopTree *loop) {
    CountedLoopNode *head = loop->_head->as_CountedLoop();

    uint storagesExpected = head->unrolled_count();
    //We want exactly 2 storages
    if (storagesExpected != 2) return NULL;

    //Initialize the Find Storage DFS with the number of expected storages
    FindStorageDFS dfs(storagesExpected);

    Arena *a = loop->_phase->arena();

    //Index of the increment node of the loop
    Node *incNode = head->incr();
    //Number of stores actually found
    uint storesFound = 0;

    //STEP 1: Find the right number of array storages
    //We need to find exactly a number of storages equal to the loop->unrolling_count() 
    //which can be accessed from the loop->incr() node.
    //This means that there is an array being indexed by a function of i.
    //It maybe more. In that case, reject to be conservative.
    //It may be also zero, in which case this is not a signal loop
    uint i = 0;

#ifdef TRACE_DEGRADE  
    tty->print_cr("STORAGES EXPECTED: %d", storagesExpected);
#endif 

    while (i < loop->_body.size() &&
            dfs._storageFound <= storagesExpected) {
        StoreNode *snn = loop->_body.at(i)->isa_Store();

        if (snn && snn->raw_adr_type() != NULL &&
                (snn->raw_adr_type()->base() == Type::Array ||
                snn->raw_adr_type()->base() == Type::AryPtr)) {
            //We found an an array storage:
            //STEP 2: See if the Array storages we found is connected to the INCREMENT
            //Search the INCREMENT node. If we can reach it 
            //from here, then increment the number of stores found
            //TODO: Can be done much more efficiently if we travel up the
            //      nodes representing the <Base + i * SizeOf(type)> to address
            //      the storage
#ifdef TRACE_DEGRADE  
            tty->print_cr("ARRAY STORAGES FOUND");
#endif    
            //TODO: If no need for DFS later, go back to original code
            dfs.search(snn, incNode, a);
        }
        i++;
    }

#ifdef TRACE_DEGRADE  
    tty->print_cr("STORAGES FOUND: %d - EXPECTED: %d",
            dfs._storageFound, storagesExpected);
#endif

    if (dfs._storageFound != storagesExpected) return NULL;
    return dfs._storages;
}

bool DepthFirstSearch::search(Node* begin, Node* end, Arena *a) {
    //Memory stuff
    VectorSet visited(a);
    Node_List stack(a); //Improve this?
    uint endIdx = end->_idx;
    stack.push(begin);
    bool keepSearch = true;
    while (keepSearch && stack.size() > 0) {
        Node *nn = stack.pop();
#ifdef TRACE_DEGRADE
        nn->dump();
#endif
        visited.set(nn->_idx);
        onVisited(nn);
        if (nn->_idx == endIdx) {
            onEndNode(nn);
            keepSearch = false;
        } else
            for (uint k = 0; k < nn->len(); k++) {
                Node *in = nn->in(k);
                if (in != NULL && !visited.test(in->_idx)) stack.push(in);
            }
    }
    return keepSearch;
}

//Says whether a projection can be approximable or not.
//A projection is approximable if it fulfills the following:
//1. It projects control output from a call
//2. The call data output has been deleted by DCE elimination, i.e., 
//   has only ctrl output and it used to have data output but not anymore
//3. Is the LAST 

bool DegradedUnrollingKit::can_approximate_projection(Node *n,
        Dict *callsOutCnt) {

    if (!n) return false;
    ProjNode *proj = n->isa_Proj();
    if (!proj) return false;
    CallNode *call = NULL;
    if (proj->in(0)) call = n->in(0)->isa_Call();

    //if true, upwards is approximable
    bool result = call &&
            //The projection is a control projection
            proj->len() == 1 && proj->is_CFG() &&
            //The call has lost one output and remains only control
            (intptr_t) (*callsOutCnt)[call] != 1 && call->outcnt() == 1;

    //We can ONLY remove a projection that does not gives control to a another
    //approximable projection, this ensures that we don't run into the case of
    //passing control upwards to a node that is still alive
    if (result && proj->outcnt() == 1 && proj->raw_out(0))
        call = proj->raw_out(0)->isa_Call();
    else return result;

    if (call && ((intptr_t) (*callsOutCnt)[call] != 1) &&
            call->outcnt() == 1 && call->raw_out(0))
        proj = call->raw_out(0)->isa_Proj();
    else return result;

    if (proj && proj->is_CFG()) return false;
    else return result;
}

//******************************************************************************
//Performs the degraded approximation of the loop 

bool DegradedUnrollingKit::approximate(IdealLoopTree *loop,
        PhaseIterGVN &igvn) {

#ifdef TRACE_DEGRADE  
    tty->print_cr("ABOUT TO LOOK FOR STORAGES");
#endif

    //Find the Storage nodes
    Node **stores = find_array_storages(loop);
    if (!stores) {
#ifdef TRACE_DEGRADE  
        tty->print_cr("NO STORAGES");
#endif
        return false;
    }

#ifdef TRACE_DEGRADE  
    tty->print_cr("STORAGES FOUND");
#endif

    //The very simplistic thing for the one loop of the example,
    //The Store with the biggest index is always the even and we want the in
    //to the even in the odd Store.
    //TODO: Don't do the approximation if the Input of the store have more than 
    //      one value output. There will be no perf. improvement.

    PhaseIdealLoop *phase = loop->_phase;
    Compile *C = phase->C;
    C->print_method(PHASE_TRACE_POINT_1, 1);


    Node *accurated;
    Node *approximated;
    if (stores[0]->_idx < stores[1]->_idx) {
        approximated = stores[0];
        accurated = stores[1];
    } else {
        approximated = stores[1];
        accurated = stores[0];
    }
    //Input node that will possibly die after the approximation
    Node *dead = approximated->in(StoreNode::ValueIn);
    //Can't do DCE if there will be no DCE
    //If no chances of DCE after the switch, there is no point in approximated
    if (dead->outcnt() > 1) return false;

    igvn.replace_input_of(approximated, StoreNode::ValueIn,
            accurated->in(StoreNode::ValueIn));

    Arena *a = phase->arena();

    //A stack to avoid recursion
    Node_List projStack(a);

    //Count the outputs of every call before the DCE removal
    //This is done to know what calls have lose outputs and therefore
    //potentially become approximable
    Dict *callOutCnt = new Dict(cmpkey, hashkey, a);
    for (uint k = 0; k < loop->_body.size(); k++) {
        Node *nn = loop->_body.at(k);
        if (nn->is_Call())
            callOutCnt->Insert(nn, (void*) (intptr_t) nn->outcnt(), true);
    }

    //After the replacement of outputs, a part of the graph is dead.
    //Iteratively clean this dead part
    bool progress = true;
    while (progress) {
        progress = false;

#ifdef TRACE_DEGRADE
        tty->print_cr("DEAD NODES:");
#endif
        //STEP1: Remove dead nodes
        uint k = 0;
        while (k < loop->_body.size()) {
            Node *nn = loop->_body.at(k);
            //Dead nodes that aren't already removed (req > 0)
            if (nn->outcnt() == 0) {
                //If the node is not already removed from the graph
                if (!(nn->len() <= 1 && nn->in(0) == NULL)) {
#ifdef TRACE_DEGRADE  
                    nn->dump();
#endif
                    igvn.remove_dead_node(nn);
                }
                loop->_body.remove(k);
                progress = true;
            } else k++;
        }

        //STEP2: Destroy non-void calls with just one exit (i.e. its output 
        //       is approximated)
        //TODO: This is NOT always safe. Only in pure side-effect free methods
        //      FIND A WAY TO CONFIGURE AGGRESSIVE vs conservative
        //Have into consideration the set_ctrl nodes
        //Have into consideration to give control to a live node that will stay
        //that way
        //IF AGGRESSIVE:
#ifdef TRACE_DEGRADE  
        tty->print_cr("DELETED PROJS: ");
#endif
        k = 0;
        bool proj_progress = true;
        while (k < loop->_body.size()) {
            proj_progress = false;
            //If the node is an approximated projection:
            ProjNode *proj = loop->_body.at(k)->isa_Proj();
            bool recomputeDepth = false;
            if (can_approximate_projection(proj, callOutCnt)) {
                proj->dump();
                //Pass control upwards to the next projection or loop head
                //TODO: must be the next control projection
                Node *inInCtrl = proj->in(0);
                inInCtrl = inInCtrl->in(0);
                if (inInCtrl == NULL) inInCtrl = loop->_head;

                Node *projIdom = phase->idom(proj);
                uint dd = phase->dom_depth(proj);

                //Keep the PhaseIdealLoop control structures
                for (uint j = 0; j < loop->_body.size(); j++) {
                    Node *n = loop->_body.at(j);
                    if (!n->is_CFG()) {
                        if (phase->has_node(n) && phase->get_ctrl(n) == proj)
                            loop->_phase->set_ctrl(n, inInCtrl);
                    }/* else if (proj == phase->idom_no_update(n)) {
                        phase->set_idom(n, projIdom, dd);
                        recomputeDepth = true;
                    }*/
                }

                //Remove node by making it dead
                tty->print_cr("DISCONECTED OUTPUTS:");
                while (proj->outcnt() > 0) {
                    Node *ot = proj->raw_out(0);
                    if (ot->in(0) == proj) ot->set_req(0, inInCtrl);
                    ot->dump();
                    for (uint i = 0; i < ot->req(); i++) {
                        if (ot->in(i) == proj) {
                            assert(i > 0, "This should not happen!");
                            ot->del_req(i);
                        }
                    }
                }
                igvn._print_dead = true;
                igvn.remove_dead_node(proj);
                igvn._print_dead = false;

                loop->_body.remove(k);
                proj_progress = true;
                progress = true;
                //if (recomputeDepth) phase->recompute_dom_depth();
            } else k++;
            if (k >= loop->_body.size() && proj_progress) {
                proj_progress = false;
                k = 0;
            }
        }
    }
    phase->Dominators();
    tty->print_cr("DELETE PROJS END ");
    return true;
}

//******************************************************************************
//Find which nodes can be erased using the Tarjan's algorithm
//to finds when the nodes are in a cycle (Strongly connected )
//Return the erasable nodes
//
// Tarjan's algorithm adapted from here:
// http://algs4.cs.princeton.edu/42digraph/TarjanSCC.java.html

/*void DegradedUnrollingKit::find_erasable_nodes() {


  //INIT STUFF
  
  Arena *a = _loop->_phase->arena();
  uint pre = 0; // preorder number counter
  //Resulting erasable nodes
  _erasableNodes = new Node_List(a);
  uint nodeCount = _loop->_phase->C->unique();
  uint *low = NEW_ARENA_ARRAY(a, uint, nodeCount);
  memset(low, 0, nodeCount * sizeof (uint));
  Node_List stack(a);
  VectorSet inCycle(a);
  //INIT STUFF END

  //Perform the Deep First Search only for the nodes inside the loop
  for (uint i = 0; i < _loop->_body.size(); i++) {
      Node * n = _loop->_body.at(i);
      if (CAN_VISIT(n)) do_DFS_search(stack, n, inCycle, pre, low);
      //A node can be erased if is not in a cycle inside the graph
      if (!inCycle.test(n->_idx)) _erasableNodes->push(n);
  }
#ifdef TRACE_DEGRADE
  tty->print("find_erasable_nodes. END");tty->cr();
#endif
}


//******************************************************************************
//Recursively destroy all nodes that remains useles in the event of
void DegradedUnrollingKit::do_DFS_erase(VectorSet &visited, Node *start) {


#ifdef TRACE_DEGRADE      
  level++;
  tty->cr();
  for (uint k = 0; k < level; k++) tty->print("-");
  tty->print("Class: %s ID:%2d", start->Name(), start->_idx);
#endif
    
  visited.set(start->_idx);
    
  //All nodes taking this node as input are are to be deleted as well
  for (DUIterator i = start->outs(); start->has_out(i); i++) {
      Node * out = start->out(i);
      //Perhaps not if they are Phi Nodes?
      //Perhaps not if they are connected to the CountedLoopEnd?
      if (!visited.test(out->_idx) && out->len() > out->req() - 1) { 
          do_DFS_erase(visited, out);
      } else {
#ifdef TRACE_DEGRADE                  
          tty->cr();
          for (uint k = 0; k < level; k++) tty->print("-");            
          tty->print("DEL OUT: %s ID:%2d - INS: %2d OUTS: %2d", out->Name(), out->_idx, out->len(), out->outcnt() );
#endif            
          int index = out->find_prec_edge(start);
          assert(index != -1, "Precedence Edge should exist");
          out->del_req(index);
      }
  }

  //All inputs to this node 
  while (start->len() > 0) {
      Node * in = start->in(start->len() - 1);
      start->del_req(in->_idx);
      //If the node is left with no output, then it must be deleted
      if (!visited.test(in->_idx) && in->outcnt() == 0) {
#ifdef TRACE_DEGRADE                  
          tty->cr();
          for (uint k = 0; k < level; k++) tty->print("-");            
          tty->print("DEL OUT: %s ID:%2d - INS: %2d OUTS: %2d", in->Name(), in->_idx, in->len(), in->outcnt() );
#endif            
          do_DFS_erase(visited, in);
      }
  }

#ifdef TRACE_DEGRADE                  
  tty->cr();
  for (uint k = 0; k < level; k++) tty->print("-");            
  tty->print("DESTROYING: %s ID:%2d - INS: %2d OUTS: %2d", 
          start->Name(), start->_idx, start->len(), start->outcnt() );
  level--;
#endif
  //Finally destroy the current node
  start->destruct();
     
}*/


/*
 void DegradedUnrollingKit::do_DFS_search(Node_List &stack, Node *start,
        VectorSet &inCycle, uint &pre, uint *low) {

    low[start->_idx] = ++pre;
    uint min = pre;
    stack.push(start);

#ifdef TRACE_DEGRADE   
    level++;
    tty->cr();
    for (uint k = 0; k < level; k++) tty->print("-");
    tty->print("%2d Low:%2d", start->_idx, low[start->_idx]);
#endif

    for (DUIterator_Fast imax, i = start->fast_outs(imax); i < imax; i++) {
        Node * nn = start->fast_out(i);
        if (!nn->is_Root()) { //Root node closes an SCC artifically. It does not play.
            if (CAN_VISIT(nn)) do_DFS_search(stack, nn, inCycle, pre, low);
            if (low[nn->_idx] < min) min = low[nn->_idx];
        }
    }
    if (min < low[start->_idx]) {
        low[start->_idx] = min;
#ifdef TRACE_DEGRADE  
        tty->cr();
        for (uint k = 0; k < level; k++) tty->print("*");
        tty->print("%2d Low:%2d", start->_idx, low[start->_idx]);
        level--;
#endif        
        return;
    }

#ifdef TRACE_DEGRADE  
    tty->cr();
    tty->print("%2d %2d/%2d \n", start->_idx, pre, min);
    tty->print("UNSTACKING:\n");
    uint counter = 1;
    level--;
#endif

    uint componentSize = stack.size() - 1; // Size of the connected component
    Node *nPop = NULL;
    do {
        nPop = stack.pop();
        //The node is in a cycle if the stack was bigger than one 
        //OR the top node is not the start node
        if (nPop->_idx != start->_idx || componentSize != stack.size())
            inCycle.set(nPop->_idx);
        low[nPop->_idx] = max_juint; //Don take this node as small again.
#ifdef TRACE_DEGRADE          
        tty->print("%2d. NOD: %2d/%2d - In cycle: %2d \n", counter, nPop->_idx, low[nPop->_idx], inCycle.test(nPop->_idx));
        counter++;
#endif        
    } while (nPop != NULL && (nPop->_idx != start->_idx));
}
 
 */
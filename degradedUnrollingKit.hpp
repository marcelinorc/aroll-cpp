/* 
 * File:   DegradedUnrollingKit.hpp
 * Author: elmarce
 *
 * Created on March 7, 2016, 11:36 AM
 */

#ifndef DEGRADEDUNROLLINGKIT_HPP
#define	DEGRADEDUNROLLINGKIT_HPP

#define TRACE_DEGRADE

#include "utilities/bitMap.hpp"
#include "utilities/globalDefinitions.hpp"
#include "opto/loopnode.hpp"
#include "opto/node.hpp" 


/*
 * The algorithm does the following:
 * 1 - Find the Storage to know whether this is a signal loop
 * 2 - Find nodes not in a cycle in the Graph
 * 3 - Create the LOAD node copying values from a[i] to a[i + 1]
 * 4 - Erase nodes not in a cycle.
 */

//Set of functions to perform the degraded unrolling
class DegradedUnrollingKit {    
private:
    
    enum SignalStatus { unknown_status, signal, not_signal };
    
    //Performs the DFS search for erasables starting in the node given as parameter
    //visited: Bitmap of visited nodes    
    //ids: Ids of the nodes
    
   
    //stack: Stack of visited nodes
    //start:  Node from which the DFS will start
    //pre:
    //void do_DFS_search(Node_List &stack, Node *start, 
    //        VectorSet &inCycle, uint &pre, uint *low);
    
    //Do the DFS search to erase nodes.
    //void do_DFS_erase(VectorSet &visited, Node *start);

    
    //SignalStatus _is_signal_loop;

    //Storage node where the signal calculations sinks.
    //Node *_signalStore;
    
    //Nodes selected as erasables
    //Node_List *_erasableNodes;
    
        //Indicates if the loop is a signal loop.
    //This happens when inside the loop values are assigned to an array.
    //In this case the function will return 1.
    //It can also be the case is a HARD signal loop, which happens
    //when the array storage is related to the increment node of the loop
    Node** find_array_storages(IdealLoopTree *loop);
    
    //Find which nodes can be erased using an slightly modified Tarjan's algorithm
    //that finds only whether the nodes is in a cycle (Strongly connected )
    //Return false if there is no erasable nodes or if the loop is not a signal loop
    //
    // Tarjan's algorithm adapted from here:
    // http://algs4.cs.princeton.edu/42digraph/TarjanSCC.java.html
    void find_erasable_nodes();
    
    //Says whether a projection can be approximable or not.
    //A projection is approximable if it fulfills the following:
    //1. It projects control output from a call
    //2. The call data output has been deleted by DCE elimination, i.e., 
    //   has only ctrl output and it used to have data output but not anymore
    bool can_approximate_projection(Node *n, Dict *callsOutCnt);
    
    
#ifdef TRACE_DEGRADE                
    //Tracing purposes
    uint level;        
#endif
    
public:
    DegradedUnrollingKit() { };
   
    DegradedUnrollingKit(const DegradedUnrollingKit& orig) {}
    
    virtual ~DegradedUnrollingKit() {};
    
    //Approximate the loop if it can be approximated
    //Returns true if the algorithm was able to approximate the loop
    bool approximate(IdealLoopTree *loop, PhaseIterGVN &igvn);
};

class DepthFirstSearch {

public:
    //Search for the node end, starting in begin
    virtual bool search(Node* begin, Node* end, Arena *a);
    //What to do when visiting a node while searching
    virtual void onVisited(Node *n) = 0;
    //What to do when visiting the end node
    virtual void onEndNode(Node *n) = 0;    
};

class FindStorageDFS : public DepthFirstSearch {
public:
    uint _storageFound;
    uint _storagesExpected;
    Node **_storages;
    FindStorageDFS(uint storagesExpected) : _storageFound(0),
        _storagesExpected(storagesExpected), _storages(NULL) {} 
    
    //Search for the node end, starting in begin
    bool search(Node* begin, Node* end, Arena *a) {
        if ( _storages == NULL )
            _storages = NEW_ARENA_ARRAY(a, Node*, _storagesExpected);
        _storages[_storageFound] = begin;
        return DepthFirstSearch::search(begin, end, a);
    }
    void onVisited(Node *n) {}
    void onEndNode(Node *n) { 
        _storageFound++;
        #ifdef TRACE_DEGRADE  
            tty->print_cr("CONNECTED STORAGES FOUND") ;
        #endif
    }
};

#endif	/* DEGRADEDUNROLLINGKIT_HPP */


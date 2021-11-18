/*  See accompanying .h for licence and documentation   */

/*----------------------------------------------------------------------------------------------+
|                                 PARETO FRONT DATA STRUCT                                      |
|   This structure records a list of cost/profits (cst1,val1), (cst2,val2), (cst3,val3).. s.t.  |
|      cst1<cst2<cst3.... and val1<val2<val3...                            (*)                  |
|   Any new insertion returns 0 if the new pair can not verify this property, if it is          |
|   dominated by the pair that that would stand before it. Otherwise, addIfHigherVal inserts the|
|   new pair and returns 1 (any subsequent pairs that become dominated is deleted).             |
|                                                                                               |
|   Complexity: log(nr of entries)                                                              |
|                                                                                               |
|   -  This is a pareto front of non-dominated solutions for: min cost and max profit (value)   |
|   -  You can think of this structure as a list of non-dominated states in a knapsack problem: |
|            never increase the cost (weight) without increasing the profit (value)             |
|   -  This is not a (priority) queue, because it needs some lookup operations                  |
|   -  This is not a vague set, because it needs some previous/next operations                  |
+----------------------------------------------------------------------------------------------*/

//The implementation from frontpareto1 is based on the c++ set binary search data
//structure enriched to record pairs of values and to remove dominated nodes. 
//The implementation from frontpareto2 is based on my personal Binary Search Tree. 
//This one is faster if you use small pareto fronts without much balancing needed


#ifndef FRONTPARETO2_H
#define FRONTPARETO2_H
#include<cstdlib>

//The use of this option may slow down things a bit, but it can be useful
//if you need attach some additional information to each pareto front element,
//e.g., a precedence relation to reconstruct an optimal solution in the end
#define ATTACH_INFO_TO_PAIRS

typedef int cst_t;       //What kind of cost type you have ?
                          //If you can use ints (even by multiplying double values), 
                          //the proposed Pareto frontier is faster
typedef double val_t   ;    //What kind of values you want to record ?

typedef struct bst_{      //BST = Binary Search Tree
    cst_t cst;
    val_t val; 
    int  type;              //LEFT or RIGHT or ROOT NODE
    struct bst_* left;
    struct bst_* right;
    struct bst_* up;
    #ifdef ATTACH_INFO_TO_PAIRS
    void* info = NULL;
    #endif
} bst;

class frontpareto2{
   public:
   //return 1 if the insertion is successful, 0 otherwise
   int addIfHigherVal(cst_t cst, val_t v);
   //The 8 functions below return 1 if the required pair was successful fetched
   int first(cst_t &cst, val_t& v);
   int next(cst_t  &cst, val_t& v);
   int first(val_t& v);
   int next(val_t& v);
   int last(cst_t &cst, val_t& v);
   int prev(cst_t  &cst, val_t& v);
   int last(val_t& v);
   int prev(val_t& v);
   void print_tree();
   void freeMem();


#ifdef ATTACH_INFO_TO_PAIRS
   bst* lastAdded = NULL;
   void putInfoOnLastAdded(void *);
   void* getInfoCurrElem(); 
   //scan all pairs and apply the provided attached info eraser on each one
   void freeMem(void(*newCallbackFreeAttached)(void*));
   void freeMem(void(*newCallbackFreeAttached)(void*),bst*);
   //If you call freeMem() with no arguments, you had better also use func below
   //so to free attached info memory yourself even by calling freeMem(). Example:
   //void eraseTransitionInfo(void* ptr)
   //{
   //     delete static_cast<transition*>(ptr);
   //}
   //<objname>.setCallbackFreeAttachInfo(&eraseTransitionInfo);
   //<objname>.freeMem() //calls eraseTransitionInfo on each pair
   void setCallbackFreeAttachInfo(void(*newCallbackFreeAttached)(void*));
#endif


   private:
   bst* root = NULL;
   bst*iterator;         //for first/next, as well as last/prev
   bst* new_leaf(cst_t w, val_t p, bst* up, int type);
   bst* dive_left(bst* r);
   bst* dive_right(bst* r);
   bst* next(bst* r);
   bst* prev(bst* r);
   void remove_node(bst* r);
   void filter(bst* r);
   int addIfHigherVal(bst* r, cst_t cst, val_t v);
   void liste_print(bst* r, char* prefix, int left_child);
   void freeMem(bst*);

};
#endif

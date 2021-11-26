/*See file LICENSE at the root of the git project for licence information*/


/*----------------------------------------------------------------------------------------------+
|                                 PARETO FRONT DATA STRUCT                                      |
|   This structure holds a list of cost/profits (cst1,val1), (cst2,val2), (cst3,val3).. s.t.    |
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

#ifndef FRONTPARETO1_H
#define FRONTPARETO1_H

#define ATTACH_INFO_TO_PAIRS

typedef int cst_t;        //What kind of cost type you have ?
                          //If you can use ints (even by multiplying double values), 
                          //the frontpareto1 is faster
typedef double val_t;     //What kind of values you want to record ?
class dpointer;           //this d-pointer hides internals (see cheshire cat), tests show that 
                          //eliminating it does not necessarily speed-up the code significantly

class frontpareto1{
   public:
   //return 1 if the insertion is successful, 0 otherwise
   int addIfHigherVal(cst_t cst, val_t v);
   int first(cst_t &cst, val_t& v);
   int next(cst_t  &cst, val_t& v);
   int first(val_t& v);
   int next(val_t& v);
   void freeMem();
   void setDeleteCallback(void(*newCallbackOnDelete)(cst_t, val_t));
       //use function above if you want to act when a pair is deleted. Here is an example.
       //void test(int a, int b){
       //    cout<<"Confirm delete at cst="<<a<<","<<"val="<<b<<endl;
       //}
       //...
       //<objname>.setDeleteCallback(&test);
   void printPairs();

   frontpareto1();
   ~frontpareto1();
#ifdef ATTACH_INFO_TO_PAIRS
   void putInfoOnLastAdded(void *);
   void* getInfoCurrElem(); 
   //scan all pairs and apply the provided attached info eraser on each one
   void freeMem(void(*newCallbackFreeAttached)(void*));
   //If you call freeMem() with no arguments, you had better also use func below
   //so to free attached info memory yourself even by calling freeMem(). Example:
   //void eraseTransitionInfo(void* ptr)
   //{
   //     delete static_cast<transition*>(ptr);
   //}
   //<objname>.setCallbackFreeAttachInfo(&eraseTransitionInfo);
   //<objname>.freeMem() //calls eraseTransitionInfo on each pair
   void setCallbackFreeAttachInfo(void(*newCallbackFreeAttached)(void*));
   void keepDeleted(void*currInfo);        //an internal: moves info to deleted list
#endif
   private:
   dpointer*dptr;
};
#endif

/*//Try code below:
#include<iostream>
using namespace std;
void test(cst_t a, val_t b){
   cout<<"Confirm delete at cst="<<a<<","<<"val="<<b<<endl;
}
int main(){
    //You could also keep the values in a classical table/array. Observe that
    //addIfHigherVal(a,b) returns 1 if the insertion (a,b) can be made. In
    //this case you can state table[a] = b and find the value table[a]=b later.
    //You may need a sparse set to keep all values like above a, values for
    //which table[a] means something. For this, check Briggs and Torczon's  
    //1993 paper "An Efficient Representation for Sparse Sets", or easier 
    //here: research.swtch.com/sparse  . For general tree considerations, 
    //see "Sec. 2.2 A Mixed Tree-Array Data Structure for Fast Recording of
    //Non-Dominated States" in the ARP paper: "Besides using an array of size
    //Q to record all states of v, we use a logarithmic-time red-black tree
    //that maintains a sorted list of states that verify above conditions."

    frontpareto1 rl;
    frontpareto1 rl2;
    rl.setDeleteCallback(&test);
    rl.addIfHigherVal(30,9);
    rl.addIfHigherVal(28,49);
    rl2.addIfHigherVal(1,4);
    rl.addIfHigherVal(10,7);
    rl.addIfHigherVal(30,60);
    rl.addIfHigherVal(28,49);
    rl2.addIfHigherVal(1,49);
    rl2.addIfHigherVal(1,50);
    #ifdef ATTACH_INFO_TO_PAIRS
    int* testInfo = new int ;
    *testInfo=7777;
    rl2.putInfoOnLastAdded((void*)testInfo);
    #endif
    rl2.addIfHigherVal(2,51);
    rl.addIfHigherVal(28,49);
    rl.addIfHigherVal(30,9);
    rl2.addIfHigherVal(11,8);
    rl.addIfHigherVal(1,8);
    rl2.addIfHigherVal(1,9);
    rl2.addIfHigherVal(1,10);

    rl.printPairs();
    cout<<"****"<<endl;
    rl2.printPairs();
    cout << '\n';

    cst_t iii; double vvv;
    //cout<<"Success of getting first of rl:"<<rl.first(iii,vvv)<<endl;
    //cout<<iii<<":"<<vvv<<endl;
    //cout<<"Success of getting second of rl:"<<rl.next(iii,vvv)<<endl;
    //cout<<iii<<":"<<vvv<<endl;

    for(int cont=rl.first(iii,vvv);cont;cont=rl.next(iii,vvv)){
        cout<<iii<<":"<<vvv<<endl;
    }
    rl.freeMem();    //move it above the for loop to see it works
    cout<<"****"<<endl;
    for(int cont=rl2.first(iii,vvv);cont;cont=rl2.next(iii,vvv)){
        cout<<iii<<":"<<vvv;
        #ifdef ATTACH_INFO_TO_PAIRS
        if(rl2.getInfoCurrElem()!=NULL){
            cout<<"Info:"<<*(int*)rl2.getInfoCurrElem();
        }
        #endif
        cout<<endl;
    }
    rl2.addIfHigherVal(1,100);
    rl2.freeMem();    //move it above the for loop to see it works
    
    cout<<"**** list 3 ****"<<endl;
    frontpareto1 rl3;
    rl3.setDeleteCallback(&test);
    rl3.addIfHigherVal(0.1,4);
    rl3.addIfHigherVal(10.2,7);
    rl3.addIfHigherVal(28,49);
    rl3.addIfHigherVal(38,5);
    rl3.addIfHigherVal(10.1,60);
    rl3.printPairs(); 
    return 1;
}
//*/

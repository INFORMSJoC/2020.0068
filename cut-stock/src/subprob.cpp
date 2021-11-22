/*--------------------------------------------------------------------------------------+
| Functions to solve random base intersection (projection) and separation sub-problems  |
--------------+---------------------------------------------------------------+---------+
              | Author: @ Daniel Porumbel 2021                                |
              |License: Any person obtaining a copy of this code is free to   |
              |         use it in any manner, subject to two conditions:      |             
              |           1) no profit may ever be made from using this code  |
              |           2) these 5 lines of text shall be included          |
              +--------------------------------------------------------------*/
#include "subprob.h"

//Below one can choose between two frontpareto implementations.
//By default, we use frontpareto version 2 (FPARETO2) when no compilation option is given
#if !defined (FPARETO1) && ! defined (FPARETO2) 
    #define FPARETO2
#endif 

#ifdef FPARETO1
    #include "frontpareto1.h"
    #define PARETOCLASS frontpareto1
#endif

#ifdef FPARETO2
    #define PARETOCLASS frontpareto2
    #include "frontpareto2.h"
#endif


#include<climits>
#include<cassert>
#include<algorithm>            //for sort
#include<cmath>
#include<cstdlib>
#include<iostream>
#include<iomanip>
using namespace std;


//#define DP_SCAN_ALL_W_RANGE  //Use this to scan all range [0..extC/C] in DP, even for 
                               //weight values not associated to patterns of that weight.
#define EPSILON 1.0e-6

//put them inside namespace
double* bst_xbase = NULL;      //the best feasible solution found so far, truncated
double bst_xbase_val = 0;
namespace{

PARETOCLASS* states;       

int*         prev=NULL;        //when not all range [0..extC] is used, you scan the range
int          last;             //by moving from prev in prev starting with last
struct transition{             //transitions between states:
    int article;               //article used to perform a transition to current state
    transition * prec;         //prec is the precedent state
};


void eraseTransitionInfo(void* ptr)
{
     delete static_cast<transition*>(ptr);
}
void correlatePrevs(int lastNew, int*prevNew)
{
    int i = last;
    int j = lastNew;
    if(lastNew>last){
        i       = lastNew;
        prev[i] = last;
        last    = lastNew;
    }
    while(true){
        assert(j<=i);
        if(i==j)
            j = prevNew[j];
        while(prev[i]<j){
            prev[j] = prev[i];
            prev[i] = j;          //i  ---> prev[i]=j  ---->prev[i]
            i       = j;          //i=prev[i]
            j       = prevNew[j];
        }
        if(j==-1)
            break;
        while(j<=prev[i])
            i = prev [i];
    }
}
/*---------------------------------------------------------------------------------------+
|                                                                                        |
|                       Random base intersection sub-problem                             |
|                                                                                        |
+---------------------------------------------------------------------------------------*/
#define USE_TOUCHED                             /*can bring a minimal speed-up like 1-2%*/

//              c_a - x^Ta
//Return t=min ----------- , over all patterns a where c_a is the ptn cost
//               y^T a
//We'll have (x^T+ty^T )a <= c_a for all patterns a
double generalInter(double *x, double* y, double * newCut, double&rHnd)
{
    clog<<"*********         Start gen inter alg            *********\n";
    #ifndef NDEBUG
    for(int i=0;i<n;i++){
        if(!(abs(x[i]-TRUNC(x[i]))<EPSILON)){
            clog<<x[i]<<"vs truncated version "<<TRUNC(x[i])<<endl;
            clog<<"not right. I need truncated input. Now need exit"<<endl;
            exit(1);
        }
    }
    #endif


    /*------------                       Init Data                 --------------*/
    transition* tNew;
    states = new PARETOCLASS [(int)extC+1]();  //() not necessary, default constructor called 
                                              //by default on not built-in types
    #ifndef DP_SCAN_ALL_W_RANGE
    static int*prevNew = new int[(int)extC+1];
    int lastNew;
    if(prev==NULL)
        prev = new int[(int)extC+1];
    last   = 0;
    prev[last] = -1;
    #endif
    assert(extC==floor(extC));
    assert(extC==ceil(extC));
    #ifdef USE_TOUCHED
    static int*touched=new int[(int)extC+1];
    for(int i=0;i<(int)extC+1;i++)
        touched[i]=0;
    #endif

    //initial state
    states[0].addIfHigherVal(TRUNC_MLT(EXT_F(0)), 0);
    tNew           = new transition;
    tNew->article  = -1;
    tNew->prec     = NULL;
    states[0].putInfoOnLastAdded((void*)tNew);


    /*------------                MAIN DP scheme calc              --------------*/
    for(int i=0;i<n;i++)
      if( (x[i]!=0) || (y[i]!=0) )
        for(int mult=0;mult<b[i];mult++){
            #ifdef DP_SCAN_ALL_W_RANGE
            double tmp;
            for(int basew = extC; basew>=0; basew--)
                if (states[basew].first(tmp))                     {
            #else
            int j = -1;
            lastNew = 0;
            for(int basew = last; basew >= 0; basew = prev[basew]){
            #endif
                int neww = basew+w[i];
                if(neww<=extC){
                      //Generate all new states at neww weight
                      double baseProf=-1; 
                      int    baseCost=-1;
                      int    deltaCost = -TRUNC_MLT(EXT_F(basew/C)) 
                                         +TRUNC_MLT(EXT_F(neww/C))-TRUNC_MLT(x[i]);
                      int    cont;              //signals if below for can continue
                      for(cont=states[basew].first(baseCost,baseProf); 
                                  cont; cont=states[basew].next(baseCost,baseProf)){
                            int newCost = baseCost + deltaCost;
                            #ifndef NDEBUG
                            if(newCost<0)
                                clog<<"Problem newCost="<<newCost<<endl;
                            assert(newCost>=0);
                            #endif
                            double newProf = baseProf+y[i];
                            if(states[neww].addIfHigherVal(newCost,newProf)){
                                 tNew = new transition;
                                 tNew->article = i;
                                 tNew->prec = (transition*)states[basew].getInfoCurrElem();
                                 states[neww].putInfoOnLastAdded((void*)tNew);
                            }
                      }
                      #ifndef DP_SCAN_ALL_W_RANGE
                          if(lastNew == 0){
                              lastNew = neww;
                              j       = lastNew;
                          }else
                              prevNew[j] = neww;
                          //can use touched to skip j=neww, reducing prevNew
                          #ifdef USE_TOUCHED
                          if((!touched[neww]))
                          #endif
                              j = neww;
                      #endif /* ndef DP_SCAN_ALL_W_RANGE */
                } /*end neww<=extC */
            } /* end scan basew */

            #ifndef DP_SCAN_ALL_W_RANGE
                assert(j>=0);
                prevNew[j] = -1;
                correlatePrevs(lastNew,prevNew);
                #ifdef USE_TOUCHED
                    for(int basew = lastNew; basew >= 0; basew = prevNew[basew])
                        touched[basew] = 1;
                #endif
            #endif
        } /* end scan all multiplicities of article i */


    /*---------------              find optimum state             ---------------*/
    transition* tranBest=NULL;
    double bestProf  = INT_MIN;
    int    bestCost  = 1;
    double profLocl  = INT_MIN;
    int    costLocl  = INT_MIN;
    int    bestw     = -1;
    #ifdef DP_SCAN_ALL_W_RANGE
    for(int basew = extC; basew>=0; basew--)              {
    #else
    for(int basew = last; basew >= 0; basew = prev[basew]){
    #endif
        for(int cont=states[basew].first(costLocl,profLocl); 
                                cont; cont=states[basew].next(costLocl,profLocl)) {
                assert(costLocl>=0);
                if(profLocl<=EPSILON)//or <=EPSILON?   //tStar = infty for current state
                    if(bestProf==INT_MIN){
                        bestProf = profLocl;
                        bestCost = costLocl;
                        tranBest = (transition*)states[basew].getInfoCurrElem();
                        bestw    = basew;
                    }
                if(profLocl>EPSILON){//or >EPSILON
                    if(costLocl==0){
                        bestProf = profLocl;
                        bestCost = costLocl;
                        tranBest = (transition*)states[basew].getInfoCurrElem();
                        bestw    = basew;
                        goto found_best_state;         //tStar = 0, can not be lower
                    }else{
                        if( (bestProf<=EPSILON) //or <=EPSILON?
                                         ||(bestProf/bestCost < profLocl/costLocl)){
                            bestProf = profLocl;
                            bestCost = costLocl;
                            tranBest = (transition*)states[basew].getInfoCurrElem();
                            bestw    = basew;
                        }
                    }
                }
            //clog<<basew<<":"<<profLocl<<"   ";
        } /* end for loop scanning states of basew */
    }
    //clog<<endl;
    found_best_state:
    clog<<"bestCost (multiplied)="<<bestCost<<", bestProf="<<bestProf<<" at bestw="<<bestw<<endl;


    /*===         Fill newCut using precedence relations between states       ===*/
    for(int i=0;i<n;i++)
        newCut[i] = 0;
    rHnd = EXT_F(bestw/C);

#ifdef NDEBUG
    while(tranBest->article!=-1){               //article=-1 only in states[0]
        newCut[tranBest->article]++;
        tranBest         = tranBest->prec;
    }
#else
    double verifyProfit = 0;
    double verifyCost   = rHnd;
    clog<<"Using articles: ";
    while(tranBest->article!=-1){               //article=-1 only in states[0]
        newCut[tranBest->article]++;
        clog<<","<<tranBest->article;
        verifyProfit += y[tranBest->article];
        verifyCost   -= x[tranBest->article];
        tranBest         = tranBest->prec;
    }
    clog<<endl;
    clog<<"   verify profit="<<verifyProfit<<" calculated profit="<<bestProf<<endl;
    assert(abs(verifyProfit-bestProf)<EPSILON);
    clog<<"   verify cost non trunc mult="<<(verifyCost)<<":";
    clog<<"   verify cost multiplied="<<TRUNC_MLT(verifyCost)<<endl;
    assert(abs(TRUNC_MLT(verifyCost)-bestCost)<EPSILON);
#endif /*NDEBUG*/


    /*-----------------                 free all mem         --------------------*/
    #ifdef DP_SCAN_ALL_W_RANGE
    for(int basew = extC; basew>=0; basew--)              
        states[basew].freeMem(&eraseTransitionInfo);
    #else
    for(int basew = last; basew >= 0; basew = prev[basew])
        states[basew].freeMem(&eraseTransitionInfo);
    #endif
    delete[] states;
    if(bestProf<=0)                                  //open direction, quite strange
        return INT_MAX;                              //for this problem
    return ((double)bestCost)/(bestProf*TRUNC_FACT);
}

/*---------------------------------------------------------------------------------------+
|                                                                                        |
|                  Separation sub-problem, extended (mult-len) knapsack                  |
|                                                                                        |
+---------------------------------------------------------------------------------------*/

//Almost everything is global below, don't use ::
//cost is the cost of the pattern material, arising in rHnd of constraint
//returns the best profit including cost, i.e., profit/value of articles - pattern cost
double extendedKnapskDP(double *p, double * newCut, double&rHnd)
{
    clog<<"---------         Start sep alg            ---------\n";
    /*============                       Init Data                 ==============*/
    transition* tNew;
    states = new PARETOCLASS [(int)extC+1]();//() not necessary, default constructor 
                                            //called by default on not built-in types
    #ifndef DP_SCAN_ALL_W_RANGE
    int lastNew;
    int*prevNew;
    prev   = new int[(int)extC+1];
    prevNew= new int[(int)extC+1];
    last   = 0;
    prev[last] = -1;
    #endif
    assert(extC==floor(extC));
    assert(extC==ceil(extC));

    //initial state
    states[0].addIfHigherVal(TRUNC_MLT(EXT_F(0)),0-EXT_F(0));
    tNew           = new transition;
    tNew->article  = -1;
    tNew->prec     = NULL;
    states[0].putInfoOnLastAdded((void*)tNew);

    /*------------                MAIN DP scheme calc              --------------*/
    for(int i=0;i<n;i++)
      if(p[i]!=0)
        for(int mult=0;mult<b[i];mult++){
            #ifdef DP_SCAN_ALL_W_RANGE
            double tmp;
            for(int basew = extC; basew>=0; basew--)
                if (states[basew].first(tmp))                     {
            #else
            int j = -1;
            lastNew = 0;
            for(int basew = last; basew >= 0; basew = prev[basew]){
            #endif
                int neww = basew+w[i];
                if(neww<=extC){
                    double baseProf=-1; states[basew].first(baseProf);
                    int    newCost = TRUNC_MLT(EXT_F(neww/C));
                    double newProf = baseProf+EXT_F(basew/C)+p[i]-EXT_F(neww/C);
                    if(states[neww].addIfHigherVal(newCost,newProf)){
                         tNew = new transition;
                         tNew->article = i;
                         tNew->prec    = (transition*)states[basew].getInfoCurrElem();
                         states[neww].putInfoOnLastAdded((void*)tNew);
                    }
                    #ifndef DP_SCAN_ALL_W_RANGE
                    if(lastNew == 0){
                        lastNew = neww;
                        j       = lastNew;
                    }else{
                        prevNew[j] = neww;
                    }
                    //can use touched to skip j=neww, reducing prevNew
                    j = neww;
                    #endif
                }
            }
            #ifndef DP_SCAN_ALL_W_RANGE
            assert(j>=0);
            prevNew[j] = -1;
            correlatePrevs(lastNew,prevNew);
            #endif
        }

    /*---------------              find optimum state             ---------------*/
    transition* tranBest;
    int bstWeight = 0;
    double bestProf = INT_MIN;
    double profLocl = INT_MIN;
    #ifdef DP_SCAN_ALL_W_RANGE
    for(int basew = extC; basew>=0; basew--)              {
    #else
    for(int basew = last; basew >= 0; basew = prev[basew]){
    #endif
        if(states[basew].first(profLocl)){
            if(profLocl>bestProf){
                bstWeight  = basew;
                bestProf   = profLocl;
            }
            //clog<<basew<<":"<<profLocl<<"   ";
        }
    }
    //clog<<endl;
    clog<<"bstStateIdx="<<bstWeight<<"of profit-cost (rHand)"<<bestProf<<endl;
    assert(profLocl>=INT_MIN);

    /*------      Fill newCut using precedence relations between states    ------*/
    for(int i=0;i<n;i++)
        newCut[i] = 0;
    rHnd = EXT_F(bstWeight/C);
    states[bstWeight].first(profLocl);
    tranBest = (transition*)states[bstWeight].getInfoCurrElem();
    double verifyProfit = 0;
    clog<<"Using articles: ";
    while(tranBest->article!=-1){               //article=-1 only in states[0]
        newCut[tranBest->article]++;
        clog<<","<<tranBest->article;
        verifyProfit += p[tranBest->article];
        tranBest         = tranBest->prec;
    }
    clog<<"\nverify profit-cost="<<verifyProfit-rHnd<<endl;
    assert(abs(verifyProfit-rHnd-bestProf)<EPSILON);

    /*-----------------                 free all mem         --------------------*/
    #ifdef DP_SCAN_ALL_W_RANGE
    for(int basew = extC; basew>=0; basew--)              
        states[basew].freeMem(&eraseTransitionInfo);
    delete[] states;
    return bestProf;
    #else
    for(int basew = last; basew >= 0; basew = prev[basew])
        states[basew].freeMem(&eraseTransitionInfo);
    delete[] states;
    delete[] prev;
    delete[] prevNew;
    return bestProf;
    #endif
}

/*-----------------------+------------------------------------+--------------------------
                         |            SORT FUNCTIONS          |
                         +-----------------------------------*/
//x[i] = x[f[i]] 
template<typename T>
inline void reshuffle(T*x, int*f)
{
    T *xcpy = new T[n];
    for(int i=0;i<n;i++)
        xcpy[i] = x[i];
    for(int i=0;i<n;i++)
        x[i] = xcpy[f[i]];
    delete[] xcpy;
}
//x[f[i]] = x[i] 
template<typename T>
inline void reshuffleinv(T*x, int*f)
{
    T *xcpy = new T[n];
    for(int i=0;i<n;i++)
        xcpy[i] = x[i];
    for(int i=0;i<n;i++)
        x[f[i]] = xcpy[i];
    delete[] xcpy;
}

int*    order     = NULL;
//when used with std::sort, it will sort in descendent order using comp below
bool myComp (int i,int j) { return (w[i]/(1+bst_xbase[i])>w[j]/(1+bst_xbase[j]));}

}//namespace

//The linux kernel style allows function-like macros in a do-while that defines a block
//It also states: "macros resembling functions may be named in lower case"
#define setOrdre()                    \
  do{                                 \
      if(order==NULL)                 \
          order = new int[n];         \
      for(int i=0;i<n;i++)            \
          order[i] = i;               \
      sort(order, order+n, myComp);   \
  }while(0);                           
#define generalInterShuffle(tStar, query_bs,ydirect,newCut,rHnd)            \
  do{                                                                       \
    setOrdre();                                                             \
    reshuffle(ydirect,order);                                               \
    reshuffle(query_bs,order);                                              \
    reshuffle(w,order);                                                     \
    reshuffle(b,order);                                                     \
    tStar = generalInter(query_bs,ydirect,newCut,rHnd);                     \
    reshuffleinv(ydirect,order);                                            \
    reshuffleinv(query_bs,order);                                           \
    reshuffleinv(w,order);                                                  \
    reshuffleinv(b,order);                                                  \
    reshuffleinv(newCut,order);                                             \
  }while(0);

/*---------------------------------------------------------------------------------------+
|                                                                                        |
|                        Functions that call intersection/separation                     |
|                                                                                        |
+---------------------------------------------------------------------------------------*/

//The linux kernel style allows function-like macros in a do-while that defines a block
//It also states: "macros resembling functions may be named in lower case"
#define clogPrintQueryPnt(p)                                 \
  do{                                                        \
    clog<<"query point (curr opt sol, an ub)"<<":";          \
    int printedVals = 0;                                     \
    for(int i=0;i<n;i++)                                     \
        if(p[i]!=0){                                         \
            clog<<p[i]<<"/"<<i<<" ";                         \
            printedVals++;                                   \
            if(printedVals==5){                              \
                clog<<"  ...  ";                             \
                if(i<n/2){                                   \
                    clog<<p[n/2]<<"/"<<n/2<<" ";             \
                    clog<<"  ...  ";                         \
                    i = n-i;                                 \
                }                                            \
            }                                                \
        }                                                    \
    clog<<endl;                                              \
  }while(0);

bool sepClassicalCalcBounds(double *p, double * newCut, double&rHnd)
{
#ifndef NDEBUG
    clogPrintQueryPnt(p);
#endif
    double prof_min_rHnd = extendedKnapskDP(p,newCut,rHnd);
    clog<<"prof-rHnd of classical DP knapsack="<<prof_min_rHnd<<endl;

    //Lagrangian bound for (multiple len) csp
    double minRedCost = -prof_min_rHnd; 
    ::lowerBound = upperBound/(1-minRedCost*1.0/MIN_PATT_COST);
    //MIN_PATT_COST is the minimum non-zero pattern cost (1 for std cut stock)

    return(minRedCost<-EPSILON);             //separation successful
}

//solve intersection, record results in newCut and rHnd
//returns true of it is possible to separate point p
bool sepByIntersectCalcBounds(double *p, double * newCut, double&rHnd,int iter)
{
#ifndef NDEBUG
    clogPrintQueryPnt(p);
#endif
    static double* query_bs  = new double[n];             //query xbase
    static double* ydirect   = new double[n];             //ydirection
    static double* xbase     = new double[n];
    if(bst_xbase==NULL){
        bst_xbase   = new double[n];
        for(int i=0;i<n;i++)
            bst_xbase[i] = 0;
    }
    
    //seems faster that putting the conditions (iter%12...) before the loop
    for(int i=0;i<n;i++){
        query_bs[i] = bst_xbase[i]; 
        if(iter%12>=3)   
            query_bs[i]=TRUNC_DOWN(bst_xbase[i]*0.5);
        if(iter%12>=6)
            query_bs[i] = 0;
        ydirect[i] = p[i] - query_bs[i];
        //clog<<ydirect[i]<<"/"<<p[i]<<" ";
    }
    //clog<<endl;

    double tStar;
    //Main call to intersection sub-problem after first well shuffling 
    generalInterShuffle(tStar, query_bs,ydirect,newCut,rHnd);


    //Calculate current lower bound, and base of best objective value
    ::lowerBound = 0;
    double xbase_val = 0;
    for(int i=0;i<n;i++){
        xbase[i]      = query_bs[i]+tStar*ydirect[i];
        ::lowerBound += xbase[i] * b[i];
        xbase[i]      = TRUNC_DOWN(xbase[i]);
        xbase_val    += xbase[i] * b[i];
    }

    //update best feasible solution found so far.
    if(xbase_val>=bst_xbase_val){
        clog<<"update new base of val "<<xbase_val<<"\n";
        bst_xbase_val = xbase_val;
        //take xbase as bst_xbase
        delete[] bst_xbase;
        bst_xbase     = xbase;
        xbase         = new double[n];
    }

#ifndef NDEBUG
    clog<<"tStar="<<tStar<<"   ";
    clog<<"lowerBound="<<::lowerBound<<", ";
    static int nonZeros=0;
    for(int i=0;i<n;i++)
        nonZeros+=(query_bs[i]>0);
    clog<<"total nonZeros of all query points till here="<<nonZeros<<endl;
#endif
    return (tStar<=1-EPSILON);
}

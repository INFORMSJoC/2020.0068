/*----------------------------------------------------------------------------------------------+
|               Solving random base projection and separation sub-problems                      |
------------------------------------------------------------------------------------------------+
see accompanying .cpp info for licence information.*/

#ifndef SUBPROB_H_INCLUDED
#define SUBPROB_H_INCLUDED

extern double  C;                    //capacity
extern double  extC;                 //extended capacity (eg, for elastic vers, or multi-size bp)
extern int     n;                    //the number of dual variables
extern int*    b;                    //demands
extern int*    w;                    //weights
extern double  lowerBound;           //attention: defined in CuttingPlanesEngine.h
extern double  upperBound;

//try to separate point and fill newCut and rHnd, return true if success
bool sepClassicalCalcBounds(double *point, double * newCut, double&rHnd);

//calculates the query points (base and direction) and
//solves intersection, recording results in newCut and rHnd
//returns true of it is possible to separate p
//At some iterations it shoots from zeros, at others from the best base, etc.
bool sepByIntersectCalcBounds(double *point, double * newCut, double&rHnd, int iter);

/*--------------------+------------------------------------------+-----------------------
                      |FUNCTIONS BELOW EXIST IN UNNAMED NAMESPACE|
                      +------------------------------------------+

//              c_a - x^Ta
//Return t=min ----------- , over all patterns a where c_a is the ptn cost
//               y^T a
//We'll have (x^T+ty^T )a <= c_a for all ptns a
double generalInter(double *x, double* y, double * newCut, double&rHnd);
//solve extended knapsack for separation
//max p^Ta - c_a
double extendedKnapskDP(p,newCut,rHnd);
*/

/*-----------------------+------------------------------------+--------------------------
                         |  MULTI LEN CUT STOCK DEFINITION    |
                         +-----------------------------------*/
//1. Classical cutting stock : function f and maximum extension 1
//#define MIN_PATT_COST 1.0
//#define EXT_F(X) (((double)X)<=1 ? 1 : INT_MAX)
//#define EXT_MAX (1)
//2. Multiple Length Cut Stock, IPOPT is a multiple of 0.2
#define EQUAL_BOUNDS(X,Y) (ceil(5*X)==ceil(5*Y))
#define MIN_PATT_COST 0.6
#define EXT_MAX (1)
#define EXT_F(X) (((double)X)<=1 ? (((double)X)<=0.7 ? 0.6 : 1 ) : INT_MAX)
//3. Multiple Length Cut Stock Second Version
//#define EQUAL_BOUNDS(X,Y) (ceil(5*X)==ceil(5*Y))
//#define MIN_PATT_COST 0.6
//#define EXT_MAX (1)
//#define EXT_F(X) (((double)X)<=1 ?                                         \
//                  (((double)X)<=0.7? (((double)X)<=0.5 ? 0.4 : 0.6 ) : 1 ) \
//                  :                                                        \
//                  INT_MAX   /*infeasible*/                                 \
//                 )

/* other problems can be tested
    //4. Elastic/Extendable Bin Packing Quadratic
    //#define EXT_MAX (2)
    //#define EXT_F(x) ( ((double)x)<=1 ? 1 : ((double)x)*((double)x) )
    //#define EXT_F(x) ( ((double)x)<=1 ? 1 : ((double)x)*((double)x)*((double)x) )
    //5. Variabled Sized Bin Packing
    //#define EXT_MAX (2)
    //#define EXT_F(x) ( ceil(10.0*(double)x)/10 )
    //#define EXT_MAX (2)
    //#define EXT_F(x) ( ceil(2.0*(double)x)/2 )
*/

/*-----------------------+------------------------------------+--------------------------
                         |        TRUNCATION MACROS           |
                         +-----------------------------------*/
//Below I use the fact that the closest integer to x is floor(x+0.5)
#define TRUNC_FACT    5.0
//The frontpareto is faster when using integer costs corresponding to double 
//costs that are multiples of 0.2
#define TRUNC_MLT(x)  ((int)floor((x)*TRUNC_FACT+0.5))
#define TRUNC(x)      (floor((x)*TRUNC_FACT+0.5)/TRUNC_FACT)
#define TRUNC_DOWN(x) (floor((x)*TRUNC_FACT)/TRUNC_FACT)
//#define TRUNC_FACT    4.0     //Fails on hard, better on triplets (rtime, not iters)
                                //but it is more elegant to take 1/TRUNC_FACT multiple of 0.2

#endif

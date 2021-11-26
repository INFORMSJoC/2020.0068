/*----------------------------------------------------------------------------------------------+
|               Cutting Planes Engine to solve an LP adding constraints one by one              |
|                   - using the cheshire cat technique, no need to include cplex                | 
|                     .h headers in all files that include CuttingPlanesEngine.h,               |
|                     significantly speeding-up compilation                                     |    
---------+-------------------------------------------------------------------------+------------+
         | See file LICENSE at the root of the git project for licence information |
         +------------------------------------------------------------------------*/

#ifndef CUTPLANESENGINE_H_INCLUDED
#define CUTPLANESENGINE_H_INCLUDED
#include<fstream>

/*-----------------------+------------------------------------+--------------------------
                         |  General Macros for Configuration  |
                         +-----------------------------------*/


//#define CPLOG(x) ;            //definitively suppress any printing without consuming any CPU
#define CPLOG(x) cplog<<x       //allow some printing to cplog that is anyway disabled by default. 
                                //However, you could enable cplog via activateLog() or send it 
                                //to other stream via printProgressMsg(...)
//#define CPLOG(x) cplog<<__LINE__<<x

#define NO_CPLEX_OUTPUT         //disable terminal output for cplex routines

#define THREADS 1               //The number of threads to use in cplex solver. Let 0 the
                                //default and it will choose the nr of CPUs or hyper-threads of the machine

#define EPS 1.0e-6              //A specific epsilon

#define EPS_CUT_COEFS 1.0e-7    //When adding new cuts, the separator might provide very small coefficients.
                                //Any coefficient less than below value is considered zero

#define TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS (n/120 +1)  //The number of seconds before stopping a call to solve() and take the current
                                //non-optimal solution. The cutting-planes can continue and maybe even the 
                                //current non-optimal solution can be cut. Otherwise, we multiply this timeout
                                //by 100 and try a full solve. If it fails, stop.

/*--------------+----------------------------------------------------------+-------------
                |typedefs for SEPARATOR FUNCTIONS (see USAGE EXAMPLE below)|
                +---------------------------------------------------------*/

//given current solution primals (dimension nrVars), find a new cut newCut and the
//associated rightHandVal. Writing a=newCut et b=rightHandVal, the added cut is:
//  maximize=1   =>              a^T x <= b   => Return value  b - a^T x
//  maximize=0   =>              a^T x >= b   => Return value: a^T x - b
//The returned violation is a negative value if x violates newCut, positive otherwise
//If the cut separator realizes the gap is closed, it returns INT_MAX and cutPlanes stops
//Attention: one can modify primals, but it is not recommended
typedef double (*cutSeprtSimple_t) (const int nrVars, double*primals, 
                                  double * newCut, double&newRightHand);
//Same as above, but it sends to the separator the current iteration it and time tm 
typedef double (*cutSeprtSolver_t) (const int nrVars, double*primals, 
                                 double * newCut, double&newRightHand
                                 ,int it, double tm);
//A separator that can return multiple cuts
typedef double (*cutSeprtExtended_t) (const int nrVars, double*primals, 
                                  double * newCut, double&newRightHand,
                                  int it, double tm,
                                  double **newCutMore, double*newRightHandMore,
                                  int & newMore, int maxMoreLen);

extern double lowerBound;         //global variable visible in the cut separator
extern double upperBound;         //global variable visible in the cut separator
extern int    switchToIntVarsNow ;//put this 1 when the cut separator decides to switch

/*-----------------------+------------------------------------+--------------------------
                         |          MAIN CLASS                |
                         +-----------------------------------*/

class CuttingPlanesEngine{
    public:
    //constructs object using nrVars unbounded variables
    CuttingPlanesEngine(int nrVars, cutSeprtSimple_t pSimple);
    CuttingPlanesEngine(int nrVars, cutSeprtSolver_t pComplex);
    CuttingPlanesEngine(int nrVars, cutSeprtExtended_t pExtended, int maxMoreConstrToReturn);
    //classical destructor
    ~CuttingPlanesEngine();

    //by default, this makes it show progress (debug) message to clog
    void activateLog();
    //change the log stream (by default it is a disabled (with fail bit) clog)
    void printProgressMsg(std::ostream& cstream);
    //the lower and upper bound on each variable is defined by *arrays* varLB and resp varUb
    void setVarBounds(double*varLB, double*varUB);
    //the lower and upper bound on each variable is defined by doubles varLB and resp varUb
    void setVarBounds(double varLB, double varUB);
    //the lower bound on each variable is defined varLB , impose no upper bounding
    void setVarLowerBoundsOnly(double varLB);
    //Set the objective function for minimization
    //All cuts have form ...>=rightHand
    void setObjCoefsMinimize(double* coefs);
    void setObjCoefsMinimize(int*    coefs);
    //Set the objective function for maximization
    //All cuts have form ...<=rightHand
    void setObjCoefsMaximize(int*    coefs);
    void setObjCoefsMinRandomizedSolving(double* coefs,int seed);//report different opt sols
    void setObjCoefsMaxRandomizedSolving(double* coefs,int seed);//report different opt sols
    void setObjCoefsMinRandomizedSolving(int* coefs,int seed);//report different opt sols
    void setObjCoefsMaxRandomizedSolving(int* coefs,int seed);//report different opt sols
    void setObjCoefsMaximize(double* coefs);
    void setTimeoutSolve(double timeout);
    //Ask to always turn to integer at the end
    void alwaysTurnIntegerInTheEnd();
    //Gets the objective value (it does not verify that the model is fully optimized)
    double getObjVal();
    //Add to the model a new cut 
    //Return value: the number (id) of the new constraint
    //Insert a new cut of the form coefs^T x <= rightHand
    int modelAddCut(double * coefs, double rightHand);
    int modelAddCut(int*     coefs, double rightHand);
    //sense: 1 for <=, -1 for >=, or 0 for ==
    int modelAddWithSense(double * coefs, double rightHand, int sense);
    //Add to the model an equality, return id of the new constraint
    int modelAddEquality(int * coefs, double rightHand);

    //deletes constraint of nr/id i
    void modelDelCut(int i);
    //the coefficient (dual value) of the constraint of nr/id i
    double getCutDualVal(int i);
    //the righthand of cut nr/id i
    double getCutRightHand(int i);
    //returns the CPU time needed only to re-optimize after each cut (excluding the separation time)
    double getTmOnlySolve();
    //returns the number of recorded cuts
    int getNbCuts();
    //returns the number of integer variables
    int nbIntVars();
    //makes variable i integer
    void turnVarInteger(int i);
    //make all decision variables discrete, do nothing if it already the case
    void turnAllVarsInteger();
    //the contribution in the reduced cost of coef at position idx
    double coefCutContrib(int coef, int idx);
    //given newCutCoefs (coefficients of a new column), calc. \sum newCutCoefs[i]*primals[i]
    double leftHandSum(int * newCutCoefs);
    //given newCutCoefs (coefficients of a new cut) and right hand, check feasibility
    int    violatedCut(int * newCutCoefs, double rightHand);
    //solve over current columns and return the value of the objective function (restricted model)
    //set lowerBound = INT_MAX or upperBound = INT_MIN if error and return INT_MAX/INT_MIN if error
    //(not enough time as set by //setTimeoutSolve())
    double solve();
    //puts the dual variables in *yStar (if non null)
    void getPrimals(double*ystar);
    //puts xstart in the primal variables that are a starting point
    void setPrimals(double*xstart, double objValInit);
    //saves the model in a file
    void exportModel(char* filename);
    //releases/free all array memory 
    void freeData(double*&a,double*&b, double**&c,int nr);

    //Solve column generator using a cutSeprt routine called 'cutSeprt' 
    //Params itMax and tmMax indicate the maximum allowed time and iterations
    //Last 2 args (filled by the routine) indicate the nr. of iterations and time required to finish
    //Returns EXIT_SUCCESS (0) in case of success or -1 otherwise
    int runCutPlanes(const int itMax, const double tmMax, int& it, double&tm);
    int runCutPlanes(int& it, double&tm);

    void removeVar(int outIdx);
    void setVarBounds(int i, double lb, double ub);
    void addVar(double outIdx, double lb, double ub);
    void addVar(double outIdx, double lb, double ub, int andSwapLastTwoAfterAdd);
    long getTotalNrCoefs();

    private:
    std::ostream cplog;             //by default the log is disabled, use activateLog
    int n;                          //the total number of decision variables
    int intVars ;                   //the number of integer decision variables (<=n)
    int noRows;                     //number of rows (actually columns/configs)
    double * primals;
    double currObj;                 //the optimum of the restricted program 
    int * intStatus;                //indicate which variables are integer
    class cplexCheshireData & d;    //d-reference (see d-pointer, compiler firewall or cheshire data)
                                    //need no cplex include outside cplex.cpp, faster compilation
                                    //At the file bottom,you find the stack-overflow post used to it
    long totalNrCoefs;              //The total number of non-zero coefficients added by the cuts
    int maximize ;
    int primalsSetByUser ;          //If you call setPrimals, this will be set
    void init();
    cutSeprtSimple_t internalCutSeprtSimple;
    cutSeprtSolver_t internalCutSeprtSolver;
    cutSeprtExtended_t internalCutSeprtExtended;
    int maxMoreConstr;
    int turnIntegerEnd;
    double timeoutSet;
    double tmOnlySolve;
#ifdef TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS
    int suboptimal;
#endif
    //Internals:
    double runSelectedCutSeprt(const int nrVars, double*yyy, double * newRow,double&newRightHand, int it, double tm, double ** newCutMore,double*newRightHandMore, int&newMore, int maxMoreConstr);
    //an internal to set a param
    void setToleranceParamToEpsilon();
};
#endif


/*-----------------------+------------------------------------+--------------------------
                         |          USAGE EXAMPLE             |
                         +-----------------------------------*/
/*
ATTENTION : undefine USE_CPU_TIME_FUNCTION on top of this file unless you can provide the getCPUTime function
#include <iostream>
#include <climits>
#include <algorithm>
#include "CuttingPlanesEngine.h"                               
using namespace std;

//given current solution primals (dimension nrVars), find a new cut newCut and
//the associated rightHandVal. Writing a=newCut et b=rightHandVal, the added cut is:
//                               a^T x >= b
//Return value: a^T x - b, i.e., a negative value if x violates newCut, positive otherwise
//This cut example implements the cut family a^T x >=1, for all a with 2 non-zero values
double separator (const int nrVars, double*x, double * newCut, double&rightHandVal)
{
    int indexMin = 0;
    //cout<<"Current decision variables x";
    //for(int i=0;i<nrVars;i++)
    //    cout<<x[i]<<" ";
    //cout<<endl;

    for(int i=1;i<nrVars;i++)
        if(x[i]<x[indexMin])
            indexMin = i;
    int indexSecondMin=0;
    if(indexMin==0)
        indexSecondMin=1;
    for(int i=indexSecondMin+1;i<nrVars;i++)
        if(i!=indexMin)
            if(x[i]<x[indexSecondMin])
                indexSecondMin = i;

    //Finished computing the lowest two values of decision variables x
    for (int i=0;i<nrVars;i++)
        newCut[i]=0;
    newCut[indexMin]       = 1;    
    newCut[indexSecondMin] = 1;    
    rightHandVal           = 1;

    double violationDegree=x[indexMin]+x[indexSecondMin]-rightHandVal;
    cout<<"The separator returns the following cut:";
    cout<<"x["<<indexMin<<"]+x["<<indexSecondMin<<"]>=1\n";
    return violationDegree;
}

int main()
{
    int n=30;                       //number of variables
    double *b;                      //objective value
    b = new double[n];
    for (int i=0;i<n;i++)
        b[i] = 10; 

    CuttingPlanesEngine cutPlanes(n,separator);      
    cutPlanes.setVarBounds(0,100);
    cutPlanes.setObjCoefsMinimize(b);
    //cutPlanes.turnVarInteger(0);  //x[0] is integer
    //cutPlanes.turnVarInteger(1);  //x[0] is integer
    double * cutExample = new double[n];
    for(int i=0;i<n;i++) cutExample[i] = 1;
    cutPlanes.modelAddCut(cutExample,n/2);

    //cutPlanes.activateLog();       //-> see messages of cutting planes
    //cutPlanes.turnAllVarsInteger();//-> the variables are discrete, as in a Benders decomposition
    
    int itersUsed;
    double timeUsed;                 //Times might be reported by problematic clock() (wrap around issue in man 3 clock), 
                                     //when USE_CPU_TIME_FUNCTION is undefined in CuttingPlanesEngine.cpp and no getCPUTime() function is known

    if(cutPlanes.runCutPlanes(itersUsed, timeUsed)==EXIT_FAILURE)
        cout<<"\n\n ATTENTION: NOT enough time or iterations to optimize to the full!";

    double finalObj = cutPlanes.getObjVal();
    cout<<"\nFinal obj val="<<finalObj<<" obtained after "<<itersUsed<<" iterations.\n";
    double*x = new double[n]; 
    cutPlanes.getPrimals(x); 
    for(int i=0;i<n;i++)
        if(x[i]>EPS)
            cout<<"x["<<i<<"]="<<x[i]<<"; ";
    cout<<endl; 
}


//A MAXIMIZATION EXAMPLE THAT SOLVES THE CUTTING-STOCK TRIPLETS EXAMPLE, as a dual
//Any three items can make a cutting pattern. In the duals you have
//x[i]+x[j]+x[k]<=1, for any i,j,k
ATTENTION : undefine USE_CPU_TIME_FUNCTION on top of this file unless you can provide the getCPUTime function
#include <iostream>
#include <climits>
#include <algorithm>
#include "CuttingPlanesEngine.h"                               
using namespace std;

//given current solution primals (dimension nrVars), find a new cut newCut and
//the associated rightHandVal. Writing a=newCut et b=rightHandVal, the added cut is:
//                               a^T x <= b
//Return value: b- a^T x , i.e., a negative value if x violates newCut, positive otherwise
//This cut example implements the cut family a^T x >=1, for all a with 2 non-zero values
double separatorTriplets (const int nrVars, double*x, double * newCut, double&rightHandVal){ 
    int indexMax = 0;
    //cout<<"Current decision variables x";
    //for(int i=0;i<nrVars;i++)
    //    cout<<x[i]<<" ";
    //cout<<endl;
    for(int i=1;i<nrVars;i++)
        if(x[i]>x[indexMax])
            indexMax = i;
    int indexSecondMax=0;
    if(indexMax==0)
        indexSecondMax=1;
    for(int i=indexSecondMax+1;i<nrVars;i++)
        if(i!=indexMax)
            if(x[i]>x[indexSecondMax])
                indexSecondMax = i;
    int indexThirdMax=0;
    if((indexMax==0)||(indexSecondMax==0)){
        indexThirdMax=1;
        if((indexMax==1)||(indexSecondMax==1))
            indexThirdMax = 2;
    }
    for(int i=indexThirdMax+1;i<nrVars;i++)
        if((i!=indexMax)&&(i!=indexSecondMax))
            if(x[i]>x[indexThirdMax])
                indexThirdMax = i;
    //Finished computing the lowest two values of decision variables x
    for (int i=0;i<nrVars;i++)
        newCut[i]=0;
    newCut[indexMax]       = 1;    
    newCut[indexSecondMax] = 1;    
    newCut[indexThirdMax]  = 1;    
    rightHandVal           = 1;

    double violationDegree=rightHandVal-x[indexMax]-x[indexSecondMax] - x[indexThirdMax];
    cout<<"The separator returns the following cut:";
    cout<<"x["<<indexMax<<"]+x["<<indexSecondMax<<"]+x["<<indexThirdMax<<"]<=1\n";
    return violationDegree;
}
int main(){
    int n=30;             //number of variables
    double *b;            //objective value
    b = new double[n];
    for (int i=0;i<n;i++)
        b[i] = 10; 

    CuttingPlanesEngine cutPlanes(n,separatorTriplets);      
    cutPlanes.setVarBounds(0,100);
    cutPlanes.setObjCoefsMaximize(b);
    //cutPlanes.turnVarInteger(0);  //x[0] is integer
    //cutPlanes.turnVarInteger(1);  //x[0] is integer
    double * cutExample = new double[n];
    for(int i=0;i<n;i++) cutExample[i] = 1;
    cutPlanes.modelAddCut(cutExample,n/2);
    //cutPlanes.activateLog();       //-> see messages of cutting planes
    //cutPlanes.turnAllVarsInteger();//-> the variables are discrete, as in a Benders decomposition
    
    int itersUsed;
    double timeUsed;                 //Times might be reported by problematic clock() (wrap around issue in man 3 clock), 
                                     //when USE_CPU_TIME_FUNCTION is undefined in CuttingPlanesEngine.cpp and no getCPUTime() function is known
    if(cutPlanes.runCutPlanes(itersUsed, timeUsed)==EXIT_FAILURE)
        cout<<"\n\n ATTENTION: NOT enough time or iterations to optimize to the full!";
    double finalObj = cutPlanes.getObjVal();
    cout<<"\nFinal obj val="<<finalObj<<" obtained after "<<itersUsed<<" iterations.\n";
    double*x = new double[n]; 
    cutPlanes.getPrimals(x); 
    for(int i=0;i<n;i++)
        if(x[i]>EPS)
            cout<<"x["<<i<<"]="<<x[i]<<"; ";
    cout<<endl; 
}





*/
	
/// d-pointers are one implementation, among many, of the pimpl pattern. It is also one of the early
/// implementations: "The name 'd-pointer' stems from Trolltech's Arnt Gulbrandsen, who first introduced
/// the technique into Qt, making it one of the first C++ GUI libraries to maintain binary compatibility
/// even between bigger release." Source
/// 
/// One advantage of using macros is the option of changing some implementation details of the pattern
/// implementation in a central place at compile time. You could for example design your macros to leave
/// you the option of switching to the fast pimpl implementation at a later time without changing tons
/// of code (hopefully you won't need this if you are using pimpl :-)). Provided that you made no
/// mistakes in your macro design/implementation...
/// 
/// However, I would personally recommend avoiding macros for your pimpl implementation as they are
/// cryptic for any newcomer to your source tree. Macros create magical dialects that are often
/// error-prone and not as meaningful as the original source code. They also come with all the problems
/// associated with the C Pre Processor; it's unaware of the underlying language.
/// 
/// Personally I like to use what I call a d-reference. Instead of a pointer, you use a reference and
/// you don't have to d-reference. 8-) It looks something like this:
/// 
/// // MyClass.h
/// 
/// class MyClass
/// {
/// public:
///     MyClass();
///     ~MyClass();
/// 
///     // implementation methods
/// 
/// private:
///     class MyClassPrivate& d;
/// };
/// 
/// // MyClass.cpp
/// 
/// struct MyClassPrivate
/// {
///     int x;
/// };
/// 
/// MyClass::MyClass()
/// : d(*new MyClassPrivate)
/// {
/// 
/// }
/// 
/// MyClass::~MyClass()
/// {
///     delete &d;
/// }
/// 
/// // In methods use d.x
/// 
/// 


/*----------------------------------------------------------------------------------------------+
|               Cutting Planes Engine to solve an LP adding constraints one by one              |
|                   - using the cheshire cat technique, no need to include cplex                | 
|                     .h headers in all files that include CuttingPlanesEngine.h,               |
|                     significantly speeding-up compilation                                     |    
---------+-------------------------------------------------------------------------+------------+
         | See file LICENSE at the root of the git project for licence information |
         +------------------------------------------------------------------------*/

#include "CuttingPlanesEngine.h"
#include "general.h"                //only useful for getCPUTime()
#include <cstdlib>
#include <vector>
#include <cmath>
#include <ilcplex/ilocplex.h>                    
using namespace std;

#define RANDOMIZE_MASTER_SOLVE

//The cples parameters change from version to version. The macros below 
//may prevent such problems, based on a CPLEXVER that can be defined by the
//compilation line using option -D. It should contain 3 digits
#if !defined CPLEXVER ||  (CPLEXVER>=129)
    #define TIME_LIMIT_PARAM            IloCplex::Param::TimeLimit
    #define MAX_THREADS_PARAM           IloCplex::Param::Threads
    #define NUMERICAL_EMPHASIS_PARAM    IloCplex::Param::Emphasis::Numerical
    #define FEASIBLE_TOLERANCE_PARAM    IloCplex::Param::Simplex::Tolerances::Feasibility
#else
    #define TIME_LIMIT_PARAM            IloCplex::TiLim
    #define MAX_THREADS_PARAM           IloCplex::Threads
    #define NUMERICAL_EMPHASIS_PARAM    IloCplex::NumericalEmphasis
    #define FEASIBLE_TOLERANCE_PARAM    IloCplex::EpRHS
#endif



////TESTING A CERTAIN FEASIBLE PRIMAL IS REALLY FEASIBLE
//#define TESTING_A_FEASIBLE_PRIMAL_zz
#ifdef TESTING_A_FEASIBLE_PRIMAL_zz
int zz[1100];
#endif

double lowerBound;                             //a global visible anywhere
double upperBound;                             //a global visible anywhere
int    switchToIntVarsNow;                     //a global usable from cutSeprt
//#define TMP_CUT_ERASER
#ifdef TMP_CUT_ERASER
int recordedConstr = 0;
#endif
class cplexCheshireData{                       //From other files, one can refer to this struct
    public: 
    IloEnv         env;                        //without needing to include cplex headers
    IloModel       model;
    IloNumVarArray vars;
    IloCplex       cplex ;
    IloNumArray    lb ;
    IloNumArray    ub ;
    IloRangeArray cuts;                     
    #ifdef RANDOMIZE_MASTER_SOLVE
    int itersRandomized;
    double* objCoefsCopy = NULL;               //the randomizing is done by modifying original
                                               //objective function. We need a copy of it.
    #endif
    #ifdef TMP_CUT_ERASER
    IloRange rrecConstr;
    #endif
    cplexCheshireData(int n):
            env(), model(env), vars(env,n,INT_MIN,INT_MAX), cplex (env),
            lb(env,n), ub (env,n), cuts(env){
    };

};

void CuttingPlanesEngine::removeVar(int outIdx)
{
    d.vars[outIdx].end();
    d.vars.remove(outIdx);
    d.lb.remove(outIdx);
    d.ub.remove(outIdx);
    n--;
    #ifdef TMP_CUT_ERASER
    if(outIdx==recordedConstr){
        CPLOG("Erasing"<<d.rrecConstr<<endl);
        d.rrecConstr.end();
        //d.model.remove(d.rrecConstr);
    }
    #endif

}
void CuttingPlanesEngine::setVarBounds(int i, double lbnd, double ubnd)
{
    d.lb[i] = lbnd;
    d.ub[i] = ubnd;
    d.model.add(d.vars[i]>=d.lb[i]);
    d.model.add(d.vars[i]<=d.ub[i]);
}
void CuttingPlanesEngine::addVar(double objCoef, double min, double max, int swapOK)
{
    IloObjective obj = d.cplex.getObjective();
    d.vars.add(obj(objCoef));
    n++;
    d.lb.add(min);
    d.ub.add(max);

    d.model.add(d.vars[n-1]>=d.lb[n-1]);
    d.model.add(d.vars[n-1]<=d.ub[n-1]);
    if(swapOK){
        IloNumVar tmp = d.vars[n-1];
        d.vars[n-1] = d.vars[n-2];
        d.vars[n-2] = tmp; 
        double  tmpd= d.lb[n-1];
        d.lb[n-1]   = d.lb[n-2];
        d.lb[n-2]   = tmpd;
        tmpd= d.ub[n-1];
        d.ub[n-1]   = d.ub[n-2];
        d.ub[n-2]   = tmpd;
    }
    //delete[] primals;
    if(primals!=NULL)
        delete[] primals;
    primals = new double[n];
}
void CuttingPlanesEngine::addVar(double objCoef, double min, double max)
{
    return addVar(objCoef,min,max,0);
}
void CuttingPlanesEngine::init()
{
    //by default there is no log output
    cplog.setstate(std::ios_base::failbit);
    ::upperBound = INT_MAX;
    ::lowerBound = INT_MIN;
    ::switchToIntVarsNow = 0;
    timeoutSet = -1;
    #ifdef TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS
    suboptimal = 0;
    #endif

    ////TESTING A CERTAIN FEASIBLE PRIMAL IS REALLY FEASIBLE
    #ifdef TESTING_A_FEASIBLE_PRIMAL_zz
        for(int ij=0;ij<1100;ij++)zz[ij]=0;
        zz[0]=374;zz[9]=285;zz[12]=198;zz[14]=521;zz[19]=108;zz[44]=628;zz[45]=389;zz[58]=1122;zz[62]=399;zz[79]=314;zz[83]=309;zz[86]=367;zz[90]=879;zz[96]=174;zz[97]=226;zz[108]=54;zz[109]=212;zz[112]=187;zz[114]=90;zz[115]=178;zz[118]=348;zz[122]=297;zz[141]=342;zz[145]=503;
    #endif
    d.cplex.extract(d.model);

    //The number of threads cplex can use (at maximum)
    d.cplex.setParam(MAX_THREADS_PARAM,THREADS);

    //feasibility tolerance: Specifies the feasibility tolerance, the degree  
    //to which the basic variables of a model may violate their bounds.
    d.cplex.setParam(FEASIBLE_TOLERANCE_PARAM, 1e-9);
    
    d.cplex.setParam(NUMERICAL_EMPHASIS_PARAM, true);

    #ifdef NO_CPLEX_OUTPUT
    d.cplex.setOut(d.env.getNullStream());
    d.cplex.setError(d.env.getNullStream());
    //This says when it is MIP
    d.cplex.setWarning(d.env.getNullStream());
    #endif 
    primals   = new double[n];
    intStatus = new int[n];
    noRows           = 0;
    tmOnlySolve      = 0;
    intVars          = 0;
    totalNrCoefs     = 0;
    maximize         = 0;
    primalsSetByUser = 0;
}
int  CuttingPlanesEngine::nbIntVars()
{
    return intVars;
}
void CuttingPlanesEngine::setToleranceParamToEpsilon()
{
     //the first below si an old version (<Cplex 2.6) of second below
     #ifndef CPLEXVER
        d.cplex.setParam(IloCplex::EpInt, EPS);
     #else
        //CPLEXVER passed in command line via -D. It should contain 3 digits
        //see the beginning of the file
        #if (CPLEXVER>=126)
            d.cplex.setParam(IloCplex::Param::MIP::Tolerances::Integrality,EPS);//you said cplex >=12.6
        #endif
     #endif
}
void CuttingPlanesEngine::turnVarInteger(int i)
{
     if(intStatus[i]==1){
        CPLOG("Variable " << i <<" is already integer, can not re-turn it integer.");
        return;
     }
     CPLOG("Turning variable " << i <<" integer...");
     d.model.add(IloConversion(d.env, d.vars[i], ILOINT));
     setToleranceParamToEpsilon();
     intVars ++;
     intStatus[i] = 1;
     CPLOG("Done\n");
}
void CuttingPlanesEngine::turnAllVarsInteger()
{
    if(intVars==0){
        d.model.add(IloConversion(d.env, d.vars, ILOINT));
        setToleranceParamToEpsilon();
        intVars =n ;
        for(int i=0;i<n;i++)
             intStatus[i] = 1;
        return ;
    }
    if(intVars<n){
        CPLOG("\n\nMaking all variables DISCRETE******************************************************** \n\n\n");
        //d.model.add(IloConversion(d.env, d.vars, ILOINT));
        for(int i=0;i<n;i++)
            turnVarInteger(i);
    }
}

void CuttingPlanesEngine::activateLog()
{
    cplog.clear();                   //it clear fail bit, if previously set
}
void CuttingPlanesEngine::printProgressMsg(ostream& cstream)
{
    cplog.rdbuf(cstream.rdbuf());
}
CuttingPlanesEngine::CuttingPlanesEngine(int nrVars, cutSeprtSimple_t cutSeprtSimple): 
                           cplog(clog.rdbuf()),//disabled by default via failbit
                           n(nrVars), currObj(0), d(*new cplexCheshireData(n)),
                           internalCutSeprtSimple(cutSeprtSimple),
                           internalCutSeprtSolver(NULL),
                           internalCutSeprtExtended(NULL),
                           maxMoreConstr(0),
                           turnIntegerEnd(0)
{
    init();
}
CuttingPlanesEngine::CuttingPlanesEngine(int nrVars, cutSeprtSolver_t cutSeprtSolver): 
                           cplog(clog.rdbuf()),//disabled by default via failbit
                           n(nrVars), currObj(0), d(*new cplexCheshireData(n)),
                           internalCutSeprtSimple(NULL),
                           internalCutSeprtSolver(cutSeprtSolver),
                           internalCutSeprtExtended(NULL),
                           maxMoreConstr(0),
                           turnIntegerEnd(0)
{
    init();
}
CuttingPlanesEngine::CuttingPlanesEngine(int nrVars, cutSeprtExtended_t cutSeprtExtended, int maxMoreConstrToReturn): 
                           cplog(clog.rdbuf()),//disabled by default via failbit
                           n(nrVars), currObj(0), d(*new cplexCheshireData(n)),
                           internalCutSeprtSimple(NULL),
                           internalCutSeprtSolver(NULL),
                           internalCutSeprtExtended(cutSeprtExtended),
                           maxMoreConstr(maxMoreConstrToReturn),
                           turnIntegerEnd(0)
{
    init();
}

CuttingPlanesEngine::~CuttingPlanesEngine()
{ 
    d.lb.end();
    d.ub.end();
    for(int i=noRows-1;i>=0;i--)
        modelDelCut(i);
    d.cuts.end();
    d.vars.end();
    d.model.end();
    d.env.end();
    delete[] primals;
    delete[] intStatus;
}
void CuttingPlanesEngine::setVarBounds(double * varLb, double * varUb)
{
    for (int i=0;i<n;i++){
        d.model.add(d.vars[i]>=varLb[i]);
        d.model.add(d.vars[i]<=varUb[i]);
        d.lb[i] = varLb[i];
        d.ub[i] = varUb[i];
    }
    d.vars.setBounds(d.lb,d.ub);
}
void CuttingPlanesEngine::setVarBounds(double varLb, double  varUb)
{
    for (int i=0;i<n;i++){
        d.lb[i] = varLb;
        d.ub[i] = varUb;
        d.model.add(d.vars[i]>=d.lb[i]);
        d.model.add(d.vars[i]<=d.ub[i]);
    }
    d.vars.setBounds(d.lb,d.ub);
}
void CuttingPlanesEngine::setVarLowerBoundsOnly(double varLb)
{
    for (int i=0;i<n;i++){
        d.lb[i] = varLb;
        d.ub[i] = IloInfinity;
        d.model.add(d.vars[i]>=d.lb[i]);
    }
    d.vars.setBounds(d.lb,d.ub);
}
#define ADD_OBJ(MinOrMax)                                    \
    IloExpr expr(d.env);                                     \
    for(int i=0;i<n;i++)                                     \
        if(abs(coefs[i])>EPS)                                \
            expr+=coefs[i]* d.vars[i];                       \
    IloObjective obj(d.env, expr, IloObjective::MinOrMax);   \
    d.model.add(obj);                                        \
    expr.end();
void CuttingPlanesEngine::setObjCoefsMaximize(int * coefs)
{
    ADD_OBJ(Maximize);
    maximize = 1;
    currObj  = INT_MAX;
}
void CuttingPlanesEngine::setObjCoefsMaximize(double * coefs)
{
    ADD_OBJ(Maximize);
    maximize = 1;
    currObj  = INT_MAX;
}
void CuttingPlanesEngine::setObjCoefsMinimize(int * coefs)
{
    ADD_OBJ(Minimize);
    maximize = 0;
    currObj = INT_MIN;
}
void CuttingPlanesEngine::setObjCoefsMinimize(double * coefs)
{
    ADD_OBJ(Minimize);
    maximize = 0;
    currObj = INT_MIN;
}
#define setObjCopyCoefs(coefs)                               \
  do{                                                        \
    d.itersRandomized = 0;                                   \
    srand(seed);                                             \
    if(d.objCoefsCopy==NULL)                                 \
        d.objCoefsCopy = new double[n];                      \
    for(int i=0;i<n;i++)                                     \
        d.objCoefsCopy[i] = coefs[i];                        \
  }while(false);

void CuttingPlanesEngine::setObjCoefsMinRandomizedSolving(double * coefs, int seed)
{
    setObjCoefsMinimize(coefs);
    #ifndef RANDOMIZE_MASTER_SOLVE
    cout<<"Please define macro RANDOMIZE_MASTER_SOLVE in CuttingPlanesEngine.cpp\n";
    exit(EXIT_FAILURE);
    #endif
    setObjCopyCoefs(coefs);
}
void CuttingPlanesEngine::setObjCoefsMaxRandomizedSolving(double * coefs, int seed)
{
    setObjCoefsMaximize(coefs);
    #ifndef RANDOMIZE_MASTER_SOLVE
    cout<<"Please define macro RANDOMIZE_MASTER_SOLVE in CuttingPlanesEngine.cpp\n";
    exit(EXIT_FAILURE);
    #endif
    setObjCopyCoefs(coefs);
}
void CuttingPlanesEngine::setObjCoefsMinRandomizedSolving(int * coefs, int seed)
{
    setObjCoefsMinimize(coefs);
    #ifndef RANDOMIZE_MASTER_SOLVE
    cout<<"Please define macro RANDOMIZE_MASTER_SOLVE in CuttingPlanesEngine.cpp\n";
    exit(EXIT_FAILURE);
    #endif
    setObjCopyCoefs(coefs);
}
void CuttingPlanesEngine::setObjCoefsMaxRandomizedSolving(int * coefs, int seed)
{
    setObjCoefsMaximize(coefs);
    #ifndef RANDOMIZE_MASTER_SOLVE
    cout<<"Please define macro RANDOMIZE_MASTER_SOLVE in CuttingPlanesEngine.cpp\n";
    exit(EXIT_FAILURE);
    #endif
    setObjCopyCoefs(coefs);
}
void CuttingPlanesEngine::alwaysTurnIntegerInTheEnd()
{
    turnIntegerEnd = 1;
}
int CuttingPlanesEngine::getNbCuts()
{
    return noRows;
}
double CuttingPlanesEngine::getCutDualVal(int i)
{
    return d.cplex.getDual(d.cuts[i]);
}
double CuttingPlanesEngine::getCutRightHand(int i)
{
    return d.cuts[i].getLB();
}
void CuttingPlanesEngine::modelDelCut(int i)
{
    d.cuts[i].end();
}
int CuttingPlanesEngine::modelAddCut(double * coefs, double rightHand)
{
    int i;
    IloExpr expr(d.env);
    for(i=0;i<n;i++)
        if(abs(coefs[i])>EPS_CUT_COEFS){
            totalNrCoefs+=1;
            expr+=coefs[i]* d.vars[i];
        }
    #ifdef TMP_CUT_ERASER
    if((!recordedConstr)&&(coefs[n-3]>0)&&(coefs[n-3]<1)){
        recordedConstr = n-3;
        d.rrecConstr = IloRange(expr<=rightHand);
        d.model.add(d.rrecConstr);
        expr.end();
        assert(1==0);
        return -1;
    }
    #endif
    ////TESTING A CERTAIN FEASIBLE PRIMAL IS REALLY FEASIBLE
    #ifdef TESTING_A_FEASIBLE_PRIMAL_zz
       double zSum=0;
       for(i=0;i<n;i++)
           if(abs(coefs[i])>EPS_CUT_COEFS)
               zSum+=coefs[i]*zz[i];
       if(zSum<rightHand){
           CPLOG("zSum="<<zSum<<endl);
           CPLOG("rightHand="<<rightHand<<endl);
           for(i=0;i<n;i++)
               if(abs(coefs[i])>EPS){
                   CPLOG("coef[i]="<<coefs[i]<<",");
                   CPLOG("zz[i]="<<zz[i]<<endl);
               }
       }
    #endif
    CPLOG("Cut Generated:");
    #ifndef NDEBUG
    for(i=0;i<n;i++)
        if(abs(coefs[i])>EPS_CUT_COEFS){
                    CPLOG(coefs[i]<<"*y["<<i<<"]+");
        }
    CPLOG("<=");
    CPLOG(rightHand<<endl);
    #endif
    d.cuts.add(expr<=rightHand);
    noRows++;
    d.model.add(d.cuts[noRows-1]);
    expr.end();    
    return noRows-1;
}

//Linux Kernel style: "macros resembling functions may be named in lower case"
#define addCutWithSense(coefs,rightHand,sense)                                \
 do{                                                                            \
   int i;                                                                       \
   IloExpr expr(d.env);                                                         \
   for(i=0;i<n;i++){                                                            \
           expr+=coefs[i]* d.vars[i];                                           \
           totalNrCoefs++;                                                      \
   }                                                                            \
   if(sense==1)                                                                 \
       d.cuts.add(expr<=rightHand);                                             \
   if(sense==0)                                                                 \
       d.cuts.add(expr==rightHand);                                             \
   if(sense==-1)                                                                \
       d.cuts.add(expr>=rightHand);                                             \
                                                                                \
   noRows++;                                                                    \
   d.model.add(d.cuts[noRows-1]);                                               \
   expr.end();                                                                  \
   return noRows-1;                                                             \
  }while(false);

int CuttingPlanesEngine::modelAddWithSense(double * coefs, double rightHand, int sense)
{
    addCutWithSense(coefs, rightHand, sense);
}
int CuttingPlanesEngine::modelAddCut(int * coefs, double rightHand)
{
    addCutWithSense(coefs, rightHand, 1);
}
int CuttingPlanesEngine::modelAddEquality(int * coefs, double rightHand)
{
    addCutWithSense(coefs, rightHand,0);
}

double CuttingPlanesEngine::getObjVal()
{
    return currObj;
}
//Linux Kernel style: "macros resembling functions may be named in lower case"
#define randomizedSolve(cplex,model,vars,objCoefsCopy,env)                               \
  do{                                                                                    \
     IloObjective obj = cplex.getObjective();                                            \
     IloExpr exprOrigObj(env);                                                           \
     for(int i=0;i<n;i++)                                                                \
          exprOrigObj=exprOrigObj+objCoefsCopy[i]*vars[i];                               \
     IloConstraint objKeepFixed = (exprOrigObj==cplex.getObjValue());                    \
     /*cout<<setprecision(12)<<cplex.getObjValue()<<endl;*/                              \
     for(int i=0;i<n;i++)                                                                \
          obj.setLinearCoef(vars[i],rand()%100);                                          \
     model.add(objKeepFixed);                                                            \
     cplex.solve();                                                                      \
     model.remove(objKeepFixed);                                                         \
     for(int i=0;i<n;i++)                                                                \
          obj.setLinearCoef(vars[i],objCoefsCopy[i]);                                    \
     exprOrigObj.end();                                                                  \
     objKeepFixed.end();                                                                 \
     cplex.solve();                                                                      \
   }while(false);                                                                        \

//set ::lowerBound = INT_MAX (when minimizing) resp INT_MIN (when maximizing) if error 
//because (not enough time as set by setTimeoutSolve())
double CuttingPlanesEngine::solve()
{
    double objVal=INT_MIN;
    #ifdef TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS
    suboptimal = 0;
    #endif
    try{
        double startTmSolve = getCPUTime();
        //Does not work on LPs, but mainly on ILPS
        //d.cplex.setParam(IloCplex::Param::RandomSeed, 62);
        d.cplex.solve();
        #ifdef RANDOMIZE_MASTER_SOLVE
            d.itersRandomized++;
            if(d.itersRandomized<10)       //enough to randomized first 10
            if(d.objCoefsCopy!=NULL)
                randomizedSolve(d.cplex,d.model,d.vars,d.objCoefsCopy,d.env);
        #endif

        tmOnlySolve+=(getCPUTime()-startTmSolve);
        objVal = d.cplex.getObjValue();   //Extract solution
        //CPLOG("Cut Planes objVal="<<objVal<<endl);
        //CPLOG("Status:"<<d.cplex.getCplexStatus()<<endl);
        //if(d.cplex.getCplexStatus() ==IloCplex::AbortTimeLim) CPLOG("Time limit\n");
        if(d.cplex.getCplexStatus() ==IloCplex::AbortTimeLim){
            #ifdef TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS
            if(timeoutSet!=-1){//-1 means disabled
            #endif
                 cerr<<"\n\n\nATTENTION: time limit "<<timeoutSet<<" secs exceeded in ConstrGener. I will report this and\n\
                       and probably stop afterwords. Did you use setTimeoutSolve()?\n\n";
                 if(maximize){
                    ::upperBound = INT_MIN;
                    currObj      = INT_MIN;
                    return  INT_MIN;
                 }else{
                    ::lowerBound = INT_MAX;
                    currObj      = INT_MAX;
                    return  INT_MAX;
                 }
            #ifdef TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS
            }
            CPLOG("\n\n\nATTENTION: macro TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS activated in Cutting Planes \n\
                       Its value is "<<TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS<<"s. The objVal is "<<objVal<<" Explanation:\n\
                       It is the number of seconds before stopping a call to solve() and take\n\
                       the current non-optimal solution. The cutting-planes can continue\n\
                       and maybe even the current non-optimal solution can be cut.\n\
                       Otherwise, I will multiply this timeout by 100 and try a full solve.\n\
                       If it fails, I'll stop. For the moment I report the primals to cutSeprt,\n\
                       I set suboptimal=1 and don't update the lower bound"<<::lowerBound<<"\n");
            if(maximize)
                objVal = INT_MIN+1;
            else
                objVal = INT_MAX-1;
            IloNumArray solution(d.env);
            d.cplex.getValues (d.vars,solution);
            for(int i=0;i<n;i++)
                primals[i] = solution[i];
            suboptimal = 1;
            return currObj;//old lower bound, I don't update it when suboptimal
            #endif     /*TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS*/
        }
    }catch (IloCplex::Exception e){
        if(d.cplex.getStatus()==IloAlgorithm::InfeasibleOrUnbounded){
            cerr<<"Attention: InfeasibleOrUnbounded. Last time I had this error, it was because\n \
                   the program you served me was infeasible, it was not possible to find the \n \
                   configurations that do satisfy all set-covering constraints\n\
                   You can think of temporarily adding the 1 1 1 ... 1 columns with cost INT_MAX\n\
                   Later, you can erase it with modelDelCut(i) where i is the return value of \n\
                   the call modelAddCut for the above temporarily column. You do this when \n\
                   this column reaches a dual value of 0, via getCutCoefSol(...)\n\
                   Model saved to trycatch.lp";
            cerr<<"I do not call exit but set currObj=INT_MAX, because this situation might be normal. You should check this.\n"; 
            d.cplex.exportModel("trycatch.lp");
            if(maximize){
                objVal = INT_MIN;
                ::upperBound = INT_MIN;
                currObj = INT_MIN;
            }else{
                objVal = INT_MAX;
                ::lowerBound = INT_MAX;
                currObj = INT_MAX;
            }
        }
    }
    if(d.cplex.getStatus()==IloAlgorithm::Infeasible){
        d.cplex.exportModel("trycatch.lp");
        cerr<<"ATTENTION on Cutting-Planes: finding infeasible solution.  Model saved to trycatch.lp. \n";
        cerr<<"I do not call exit but set currObj=INT_MAX, because this situation might be normal.  You should check this.\n"; 
        if(maximize){
            objVal = INT_MIN;
            ::upperBound = INT_MIN;
            currObj = INT_MIN;
        }else{
            objVal = INT_MAX;
            ::lowerBound = INT_MAX;
            currObj = INT_MAX;
        }
    }

    IloNumArray solution(d.env);
    d.cplex.getValues (d.vars,solution);
    for(int i=0;i<n;i++){
            primals[i] = solution[i];
            //CPLOG(primals[i]<<",");
            //assert(primals[i]>=d.lb[i]);
            //assert(d.lb[i]==0);
            //use below if you don't set NumericalEmphasis true
            if(lowerBound!=INT_MAX)           //if still feasible
            if(!(primals[i]>=d.lb[i]-EPS)){
                cout<<"wrong i="<<i<<endl;
                cout<<primals[i]<<" "<<d.lb[i]<<endl;
                exit(1);
            }
            //rounding the primals to the correct values
            //if(primals[i]<d.lb[i]+EPS)primals[i]=d.lb[i];
    }
    //CPLOG("\n");
    
    if(maximize){
        if(objVal<currObj){                    //Update bound
            currObj      = objVal;
            ::upperBound = objVal;
        }
    }else{
        if(objVal>currObj){                    //Update bound
            currObj = objVal;
            ::lowerBound = objVal;
        }
    }
    return currObj;
}//end solve
void CuttingPlanesEngine::getPrimals(double *ystar)
{
    int i;
    if(ystar!=NULL)
        for(i=0;i<n;i++)
            ystar[i] = primals[i];
}

void CuttingPlanesEngine::setPrimals(double *startPrimals, double objValInit)
{
     IloNumArray startValues(d.env);
     IloNumVarArray varCpy(d.env);
     if(maximize)
        upperBound = objValInit;
     else
        lowerBound = objValInit;
     for (int i = 0; i < n; ++i){
             startValues.add(startPrimals[i]);
             varCpy.add(d.vars[i]);
             primals[i] = startPrimals[i];
     }
     try{
        //d.cplex.addMIPStart(varCpy, startValues);
        d.cplex.setStart(startValues,NULL,varCpy,NULL,NULL,NULL);
     }catch (IloCplex::Exception e){
        cerr<<"\n\nConcert Exception in when setting initial values for variables:"<<e.getMessage()<<endl;
        exit(EXIT_FAILURE);
     }
     varCpy.end();
     startValues.end();
     primalsSetByUser = 1;
}

double CuttingPlanesEngine::coefCutContrib(int coef, int idx)
{
    if(abs(coef)<EPS) //Useless when coef is integer
        return 0;
   return (double)coef * primals[idx]; 
}

double CuttingPlanesEngine::leftHandSum(int * coefs)
{
    double totalSum = 0;
    for (int i=0;i<n;i++)
        totalSum += coefCutContrib(coefs[i],i);
    return totalSum;
}

int CuttingPlanesEngine::violatedCut(int * coefs, double rightHand)
{
    double redCost = leftHandSum(coefs);
    if(maximize)
        return (redCost > rightHand + EPS);
    else
        return (redCost < rightHand - EPS);
}
void CuttingPlanesEngine::exportModel(char* filename)
{
    d.cplex.exportModel(filename);
}
void CuttingPlanesEngine::freeData(double*&a,double*&b, double**&c,int nr)
{
    delete[] a;
    if(internalCutSeprtExtended!=NULL){
        delete[] b;
        int i;
        for(i=0;i<nr;i++)
            delete[] c[i];
        delete[] c;
    }
}

double CuttingPlanesEngine::runSelectedCutSeprt(const int n, double*primals, double * newCut,double&newRightHand,
                                     int it, double tm, double ** newCutMore,double*newRightHandMore,
                                     int&newMore, int maxMoreConstr)
{
      if(internalCutSeprtExtended!=NULL)
           return internalCutSeprtExtended(n,primals,newCut, newRightHand,it,tm,
                                   newCutMore,newRightHandMore,newMore, maxMoreConstr);
      if(internalCutSeprtSolver!=NULL)
           return internalCutSeprtSolver(n,primals,newCut, newRightHand,it,tm);
      return internalCutSeprtSimple(n,primals,newCut, newRightHand);
} 
void CuttingPlanesEngine::setTimeoutSolve(double timeOut)
{
        //before enforcing time limits, you may need call below to be sure cplex
        //thinks in terms of CPU time and not Wall Time which includes dead moments
        //d.cplex.setParam(IloCplex::Param::ClockType,1);
        d.cplex.setParam(TIME_LIMIT_PARAM,timeOut);
        timeoutSet = timeOut;
}

int CuttingPlanesEngine::runCutPlanes(const int itMax, const double tmMax, int& it, double&tm)
{
        #ifdef TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS
        if(timeoutSet==-1){//user defined timeout has priority
            //before enforcing time limits, you may need call below to be sure cplex
            //thinks in terms of CPU time and not Wall Time which includes dead moments
            //d.cplex.setParam(IloCplex::Param::ClockType,1);
            d.cplex.setParam(TIME_LIMIT_PARAM,TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS);
        }
        #endif
  try{
       it = 0;
       tm = 0;
       double newRightHand, newViolation,startTm = getCPUTime();
       double* newCut = new double[n];
   
       double** newCutMore = NULL;        //not necessarily used 
       double* newRightHandMore=NULL;     //only for cutSeprtExtended
       int newMore;
       if(internalCutSeprtExtended!=NULL){
           newCutMore = new double*[maxMoreConstr];
           newRightHandMore = new double[maxMoreConstr];
           int iii;
           for(iii=0;iii<maxMoreConstr;iii++)
                newCutMore[iii] = new double[n];
       }

       //if(!d.cplex.isDualFeasible())                        //if unbounded
       //    setVarBounds(-IloInfinity,INT_MAX/((double)n));  //Add some bounds on the variables

       if(!primalsSetByUser)
          solve();
       primalsSetByUser = 0;
       do{
           try{
               newViolation = runSelectedCutSeprt(n,primals,newCut, newRightHand,it,tm,
                                    newCutMore,newRightHandMore, newMore, maxMoreConstr);
               if(newViolation==INT_MAX){//gap closed
                     freeData(newCut,newRightHandMore,newCutMore,maxMoreConstr);
                     tm = getCPUTime() - startTm;
                     return EXIT_SUCCESS;
               }
               #ifdef TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS
               if(suboptimal){
                    CPLOG("\n\nThe cutSeprt reports the following redCost on suboptimal solution"<<newViolation<<endl);
                    CPLOG("If this is negative, the suboptimal solution is anyway cut and I go on normally\n\n");
               }
               #endif 
           }catch (IloCplex::Exception e){
                cerr<<"\n\nConcert Exception in your cutSeprt called by the Cutting-Planes:"<<e.getMessage()<<endl;
                cerr<<"I do not call exit, because this might be normal.  You should know better, maybe catch the exception in your cutSeprt.\n\n\n"; 
                if(maximize){
                    ::upperBound = INT_MIN;
                    currObj = INT_MIN;
                }else{
                    ::lowerBound = INT_MAX;
                    currObj = INT_MAX;
                }
                freeData(newCut,newRightHandMore,newCutMore,maxMoreConstr);
                tm = getCPUTime() - startTm;
                return EXIT_FAILURE;
           }
           //max: newViolation = rightHand - a^T x
           //min: newViolation = a^T x   - rightHand
           if(newViolation >= -EPS){        //no violation
               #ifdef TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS
               if(suboptimal){
                    CPLOG("Could not cut current suboptimal solution\n");
                    CPLOG("Adding generated cuts anyways\n");
                    modelAddCut(newCut,newRightHand);
                    if(internalCutSeprtExtended!=NULL){
                         int ii;
                         for(ii=0;ii<newMore;ii++)
                             modelAddCut(newCutMore[ii], newRightHandMore[ii]);
                    }
                    CPLOG("Allowing 100 times more time, i.e., "<<100*TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS<<" seconds.\n");
                    d.cplex.setParam(TIME_LIMIT_PARAM,100*TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS);
                    solve();
                    it++;
                    if(suboptimal){
                        cerr<<"Even after allowing 100 times more time, Cutting-Planes  can not finish solve. return INFEASIBLE ";
                        cerr<<"It could be possible to use CPX_PARAM_OBJLLIM to try to gain at least one\n";
                        freeData(newCut,newRightHandMore,newCutMore,maxMoreConstr);
                        tm = getCPUTime() - startTm;
                        return EXIT_FAILURE;
                    }
                    if(maximize)
                    if(::upperBound==INT_MAX){//solve sets ::upperBound=INT_MAX if infeasible
                        cerr<<"\n\n I DID try 100 times more time. Quite strange finding infeasibility here.";
                        cerr<<"\nATTENTION: THE CONSTRAINT/COL GENERATOR LP BECOME INFEASIBLE. STOP HERE\n \
                                without calling exit because this situation might be normal. Check it. \n\n\n";
                        break;
                    }
                    if(!maximize)
                    if(::lowerBound==INT_MAX){//solve sets ::lowerBound=INT_MAX if infeasible
                        cerr<<"\n\n I DID try 100 times more time. Quite strange finding infeasibility here.";
                        cerr<<"\nATTENTION: THE CONSTRAINT/COL GENERATOR LP BECOME INFEASIBLE. STOP HERE\n \
                                without calling exit because this situation might be normal. Check it. \n\n\n";
                        break;
                    }
                    CPLOG("It seems I could solve it. Hope everything goes back to normal\n");
                    d.cplex.setParam(TIME_LIMIT_PARAM,TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS);
                    continue;
               }
               else//not suboptimal, but classical case
               #endif /*TIMEOUT_BEFORE_GOING_SUBOPTIMAL_PRIMALS*/
               {
                 if(turnIntegerEnd){
                     if(intVars==n){           //all integer => exit
                          freeData(newCut,newRightHandMore,newCutMore,maxMoreConstr);
                          tm = getCPUTime() - startTm;
                          return EXIT_SUCCESS;
                     }
                     turnAllVarsInteger();
                 }else{
                      freeData(newCut,newRightHandMore,newCutMore,maxMoreConstr);
                      tm = getCPUTime() - startTm;
                      it ++;
                      return EXIT_SUCCESS;
                 }
               }
           }
           modelAddCut(newCut,newRightHand);
           if(internalCutSeprtExtended!=NULL){
                int ii;
                for(ii=0;ii<newMore;ii++)
                    modelAddCut(newCutMore[ii], newRightHandMore[ii]);
           }
           if(::switchToIntVarsNow)  //the cutSeprt can set this global variable,
                turnAllVarsInteger();//but it can't call turnAllVarsInteger() directly
           solve();
           if(maximize)
           if(::upperBound==INT_MAX){//solve sets ::upperBound=INT_MAX if infeasible
                cerr<<"\n\n\nATTENTION: THE CONSTRAINT/COL GENERATOR LP BECOME INFEASIBLE. STOP HERE\n \
                       without calling exit because this situation might be normal. Check it. \n\n\n";
                break;
           }
           if(!maximize)
           if(::lowerBound==INT_MAX){//solve sets ::lowerBound=INT_MAX if infeasible
                cerr<<"\n\n\nATTENTION: THE CONSTRAINT/COL GENERATOR LP BECOME INFEASIBLE. STOP HERE\n \
                       without calling exit because this situation might be normal. Check it. \n\n\n";
                break;
           }
           //CPLOG("  -> "<<upBound<<endl);
           it ++;
           tm = getCPUTime() - startTm;
           CPLOG("ITER="<<it<<" TIME="<<tm<<" LBND="<<currObj<<"        //message from Cutting Planes engine"<<endl);
       }while((tm<=tmMax)&& it<=itMax);
       it++;
       freeData(newCut,newRightHandMore,newCutMore,maxMoreConstr);
       return EXIT_FAILURE;                                 //not success
  }
  catch (IloCplex::Exception e){
       d.cplex.exportModel("trycatch.lp");
       cerr<<"\n\nConcert Exception Cutting-Planes:"<<e.getMessage()<<endl;
       cerr<<"Model saved to trycatch.lp. I. I think the problem is infeasible. \n";
       cerr<<"I do not call exit, because this infeasibility/error might be normal.  You should check this.\n"; 
       return EXIT_FAILURE;
  }
}

int CuttingPlanesEngine::runCutPlanes(int& it, double&tm)
{
    return runCutPlanes(INT_MAX, INT_MAX, it, tm);
}
double CuttingPlanesEngine::getTmOnlySolve()
{
    return tmOnlySolve;
}
long CuttingPlanesEngine::getTotalNrCoefs()
{
    return totalNrCoefs;
}


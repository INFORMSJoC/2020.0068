/*----------------------------------------------------------------------------------------------+
|              Solving multiple-length cutting stock using the projection sub-problem           |
|                     - solve actually ML-CSP with both std and inter methods                   |
--------------+---------------------------------------------------------------+-----------------+
              | Author: @ Daniel Porumbel 2021                                |
              |License: Any person obtaining a copy of this code is free to   |
              |         use it in any manner, subject to two conditions:      |             
              |           1) no profit may ever be made from using this code  |
              |           2) these 5 lines of text shall be included          |
              +--------------------------------------------------------------*/
#include "general.h"
#include "CuttingPlanesEngine.h"
#include "inout.h"
#include "subprob.h"
#include<iostream>
#include<cstdlib>
#include<cassert>
#include<cstring>
#include<climits>
#include<cmath>

using namespace std;
#define EPSILON 1.0e-6

/*----------------------   Project-wide Global Variables     ---------------------*/

double    C;                 //base capacity
double extC;                 //extended capacity (eg, for elastic vers, or multi-size bp)
int       n;                 //the number of dual variables
int*      b;                 //demands
int*      w;                 //weights

/*----------------------   Global Variables Used in Main     ---------------------*/

int stdMethd         = 0;
int rndmizedRunSeed  = 0; 
int iter             = 0;
double bstLowerBound = 0;
int iterLowGap       = -1;   //iteration when ub<=bstLowerBound*1.2
double tmLowGap      = -1;   //tm for above

void checkLowGap(double lb, double ub, int iter, double tm){
    if(iterLowGap>=0)       //low gap already achieved
        return;
    if(ub<=lb*1.2){
        iterLowGap = iter;
        tmLowGap   = tm;
        cout<<"REACHED LOW GAP iter,tm="<<iter<<","<<tm<<"\n";
    }
}
int equalBounds()
{
#ifdef EQUAL_BOUNDS
    return EQUAL_BOUNDS(bstLowerBound,upperBound);
#else
    return (ceil(bstLowerBound)==ceil(upperBound)) ;
#endif
}

//We return the violation rHand-neVars^T x that is negative in case of real violation
double separator (const int nrVars, double*x, double * newCut, double&rHand)
{
    static int    lstIter       = 0;               //iter when gap was closed
    static double lstTm         = 0;               //time when gap was closed
    static double startCpuTime  = getCPUTime();
    assert(::n==nrVars);

    bool separated;
    if(stdMethd)
        separated = sepClassicalCalcBounds(x,newCut,rHand);
    else
        separated = sepByIntersectCalcBounds(x,newCut,rHand,iter);

    bstLowerBound    = max(lowerBound,bstLowerBound);
    double tmElapsed = getCPUTime()-startCpuTime;
    clog<<"            LB="<<lowerBound<<" BST LB="<<bstLowerBound<<
          " it="<<::iter<<"\n              UB="<<::upperBound<<
          " Tm="<<tmElapsed<<endl;
    //cout<<"("<<::iter<<","<<lowerBound<<")\n";//good for curve drawing with pgfplots
    checkLowGap(bstLowerBound,::upperBound,::iter,tmElapsed);

#ifndef NDEBUG
    double violation = rHand;
    clog<<"Adding ";
    for(int i=0; i<(::n); i++) {
        violation -=newCut[i]*x[i];
        if(newCut[i]>EPSILON)
            clog<<"x["<<i<<"]*"<<newCut[i]<<"+";
    }
    clog<<"<="<<rHand<<" with ";
    clog<<"violation="<<violation<<endl;
    if(separated)
        assert(violation<  -EPSILON);
    else
        assert(violation>= -EPSILON);
#endif

    if( (lstIter==0)&&(equalBounds())) {
        lstIter = ::iter;
        lstTm = getCPUTime()-startCpuTime;
        clog<<"------------>Tail cut iter="<<lstIter<<" tail cut time="<<lstTm<<endl;
        return INT_MAX;
    }

    if(!separated)
        clog<<"------------>curr opt sol is globally opt at iter="<<::iter<<endl;

    ::iter++;
    if(separated)
        return -1;
    else
        return 0;
}
//The linux kernel style allows function-like macros in a do-while that defines a block
//It also states: "macros resembling functions may be named in lower case"
//Macro below generates a value for ::lowerBound 
#define separateOnce(x,newCut,rHand,cutPlanes)                                      \
  do{                                                                               \
    if(stdMethd){                                                                   \
        upperBound = 0;                                                             \
        for(int i=0;i<n;i++)                                                        \
            upperBound+=x[i];                                                       \
        sepClassicalCalcBounds(x,newCut,rHand);       /* this uses upperBound */    \
    }else                                                                           \
        sepByIntersectCalcBounds(x,newCut,rHand,iter);                              \
    cutPlanes.modelAddCut(newCut,rHand);                                            \
    clog<<endl<<"new init cut ";                                                    \
    for(int ii=0;ii<n;ii++)                                                         \
        if(newCut[ii]!=0)                                                           \
            clog<<newCut[ii]<<"/"<<ii<<" ";                                         \
    clog<<endl<<"***************END Initial sub-problem*******************"<<endl;  \
  }while(0);

string getTmStr(double tm){
    int precision = 1;
    if(tm>=1)
        precision = 2;
    if(tm>=10)
        precision = 3;
    if(tm>=1000)
        precision = 4;
    if(tm>=10000)
        precision = 5;

    std::ostringstream sBuilder;   
    if((tm>=0.001)||(tm==0))
        sBuilder <<left<<setw(6)<<setprecision(precision)<<tm;
    else
        sBuilder <<left<<setw(6)<<"<1e-3";
    return sBuilder.str();
}
int main(int argc, char**argv)
{
#ifdef NDEBUG
    cout<<"I'll be silent, ignoring clog in release more (NDEBUG)."<<endl;
    clog.setstate(ios_base::failbit);                                       //disable clog
#endif
    if(argc==1) {
        cerr<<"Usage: ./main instance zeroIndexedInstNrInFile [-si[lent]] "
              "[-std: use standard Col Gen] [-rnd : randomized run]\n"
              "       projective cutting planes used by default, use -std to change to the standard col gen\n"
              "       to change the multiple-length variant, modify lines 40-60 "
                      " in subprob.h, you can also test variable sized bin packing\n";
        return EXIT_FAILURE;
    }
    while(argv[argc-1][0]=='-') {
        int optionsFnd = 0;
        if( (argv[argc-1][0]=='-') && (!strncmp(argv[argc-1],"-si",3))) {
            cout<<"I'll be silent, ignoring clog."<<endl;
            clog.setstate(ios_base::failbit);                               //disable clog
            argc--;
            optionsFnd++;
        }
        if( (argv[argc-1][0]=='-') && (!strncmp(argv[argc-1],"-std",4))) {
            cout<<"I'll use old classical method."<<endl;
            stdMethd=1;
            argc--;
            optionsFnd++;
        }
        if( (argv[argc-1][0]=='-') && (!strncmp(argv[argc-1],"-rnd",4))) {
            if(strlen(argv[argc-1])==4)
                rndmizedRunSeed = time(NULL);
            else
                rndmizedRunSeed = atoi(argv[argc-1]+4);
            cout<<"I'll randomize Simplex with seed "<<rndmizedRunSeed<<endl;
            argc--;
            optionsFnd++;
        }
        if(optionsFnd==0) {
            cerr<<"There is an argument starting with '-' that I can not understand\n";
            cerr<<"Accepted arguments:[-si[lent]] [-cl[assicalOldMeth]]\n";
            exit(EXIT_FAILURE);
        }
    }
    if(argc>3) {
        cerr<<"You gave me 3 arguments with no '-'. I don't understand the third.\n";
        exit(EXIT_FAILURE);
    }
    if(argc==2) {
        cout<<"I take the first instance in "<<argv[1]<<endl;
        readInstNrFromFile(0,argv[1]);
    } else
        readInstNrFromFile(atoi(argv[2]),argv[1]);
    //clog<<"C="<<C<<endl ; for(int i=0;i<n;i++) clog<<w[i]<<" "<<b[i]<<endl;

    extC = C*EXT_MAX;

    CuttingPlanesEngine cutPlanes(n,separator);
    cutPlanes.setVarBounds(0,EXT_MAX);
    if(rndmizedRunSeed>0)
        cutPlanes.setObjCoefsMaxRandomizedSolving(b,rndmizedRunSeed);
    else
        cutPlanes.setObjCoefsMaximize(b);

    int itersUsed;
    double CPUtimeUsed;                     //including time of cplex threads
    time_t start = time(NULL);

    /*--------------          Start add initial constraints        --------------*/
    static double*x = new double[n];
    static double * newCut = new double[n];
    double rHand;

    for(int i=0; i<n; i++)                
        x[i] = w[i]/extC;                  
    separateOnce(x,newCut,rHand,cutPlanes);
    if( !((stdMethd) && (MIN_PATT_COST!=1)) )        //bound not correct for std method
        bstLowerBound=max(lowerBound,bstLowerBound); //in mult-len CSP 

    for(int i=0; i<n; i++)
        x[i] = b[i];
    separateOnce(x,newCut,rHand,cutPlanes);
    //if( !((stdMethd) && (MIN_PATT_COST!=1)) )      //bound ok for std method in
    bstLowerBound=max(lowerBound,bstLowerBound);     //mult-len CSP, as b is infeasible

    clog<<"                                      -> START BEST LB="<<bstLowerBound<<endl;

    /*--------------           End  add initial constraints        --------------*/

    //Launch main Cutting Planes (dual Col Gen) Engine
    if(cutPlanes.runCutPlanes(itersUsed, CPUtimeUsed)==EXIT_FAILURE)
        cerr<<"\n\n ATTENTION: NOT enough time or iters to fully optimize!";

    /*--------------               Printing final results         ---------------*/
    double finalObj = cutPlanes.getObjVal();
    cout<<"CPU Time:"<<CPUtimeUsed<<"   Real time:"<<time(NULL)-start<<endl;
    cout<<"\nFinal obj val="<<finalObj<<" obtained after "<<itersUsed<<" iterations."<<endl;

    //Print final solution to clog
    double*xx = new double[n];
    cutPlanes.getPrimals(xx);
    for(int i=0; i<n; i++)
        if(xx[i]>EPSILON)
            clog<<"x["<<i<<"]="<<xx[i]<<"; ";
    clog<<endl;
    
    /*--------------               Printing latex results         ---------------*/
    if(iterLowGap==itersUsed)
        CPUtimeUsed = tmLowGap;      //do not count time spent on printing above info
    char* inst = new char[100];      //instance name
    strcpy(inst,argv[1]);
    inst = rindex(inst,'/')+1;       //remove folders, keep filename
    char instName[300];
    if(strchr(inst,'.')!=NULL){
        while(inst[strlen(inst)-1]!='.')
            inst[strlen(inst)-1]='\0';
        inst[strlen(inst)-1]='\0'; //remove dot
    }
    sprintf(instName,"\\texttt{%s-%d}&%g&",inst,1+((argc==2)?0:atoi(argv[2])),finalObj);

    cout<<"Tabular latex data below:\n";
    if(stdMethd==0)
        cout<<setw(40)<<instName;
    else
        cout<<setw(40)<<" ";
    cout<<setw(5)<<iterLowGap           <<"&"<<
          setw(7)<<getTmStr(tmLowGap)   <<"&"<<
          setw(5)<<itersUsed            <<"&"<<
          setw(7)<<getTmStr(CPUtimeUsed);
    if(stdMethd==0)
        cout<<"&"<<endl;
    else
        cout<<"\\\\%"<<endl;
}

/*---------------+---------------------------------------------------------------+-------------+
                 | Author: @ Daniel Porumbel 2021                                |
                 |License: Any person obtaining a copy of this code is free to   |
                 |         use it in any manner, subject to two conditions:      |            
                 |           1) no profit may ever be made from using this code  |
                 |           2) these 5 lines of text shall be included          |
                 +--------------------------------------------------------------*/
#include "inout.h"
#include "subprob.h"
#include "CuttingPlanesEngine.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cassert>
#include <climits>
#include <cstring>
#include <fstream>
using namespace std;

int      n;                       //number of variables
int      m;                       //number of initial rows
int      gamma=10;
double*  obj;
double** rows;                    //rows[...][n] = rHand, rows[...][n+1]=EQUALITY,LESS_THAN_EQ
double*  lb;
double*  ub;
double*  xbase;
double*  d;                       //used to solve proj_subprob(xbase->d)
bool     runStd;
int      seed;
int      iters;                   //iters and tmCPlanes will be filled by cutPlanes.runCutPlanes(...)
double   tmCPlanes;
double   tmSort     = 0;
int      iterLowGap = -1;         //iteration when ub<=bstLowerBound*1.2
double   tmLowGap   = -1;         //tm for above
double   nominalObj, finalObj;
char*    startsol = NULL;
int      total_multi_cuts = 0;
bool     multi_cuts_per_round;
bool     multi_cuts_limited=false;//Put a limit on the number of multiple cuts
                                  //prj: return only cuts that decrease tStar
                                  //std: maximum 10 cuts per iter

//return rHand - newRow^T x. The cut that will be added is: newRow^T x <= rHand
double sub_problem_single_cut (const int nrVars, double*x, double * newRow,double&newRHand)
{
    for(int i=0;i<(::n);i++)
        if(x[i]<lb[i])        //again numerical problems should be eliminated
            x[i]=lb[i];

    if(runStd)
        return separation(x,newRow,newRHand);
    else
        return projection(x,newRow,newRHand);
}

//return rHand - newRow^T x. The cut that will be added is: newRow^T x <= rHand
double sub_problem (const int nrVars, double*x, double * newRow, double&newRHand, 
                            int it, double tm,
                            double**newRows, double*newRHands, int& newMore, int maxNewMore )
{
    if(!multi_cuts_per_round){
        newMore = 0;
        total_multi_cuts++;
        return sub_problem_single_cut(nrVars,x,newRow,newRHand);
    }

    for(int i=0;i<(::n);i++)
        if(x[i]<lb[i])        //again numerical problems should be eliminated
            x[i]=lb[i];

    double ret_val ;
    if(runStd)
        ret_val = separation_multi(x,newRow,newRHand, newRows, newRHands, newMore, multi_cuts_limited);
    else
        ret_val = projection_multi(x,newRow,newRHand, newRows, newRHands, newMore, multi_cuts_limited);

    total_multi_cuts+=(newMore+1);
    if(total_multi_cuts>=10000){
        clog<<"\nI switch to a mono-cut sub-problem because I already have too many cuts: "<<total_multi_cuts <<endl;
        multi_cuts_per_round = false;
    }
    return ret_val;
}

void checkAllParams(int& argc, char** argv)
{
    int verbose = 0;
    multi_cuts_per_round = false;
    seed = 0;                                               //no randomization
    while(argv[argc-1][0]=='-'){
        int argc_start = argc;
        if(argv[argc-1][1]=='v') {
            verbose = 1;
            argc--;
        }
        if(argv[argc-1][1]=='m') {
            multi_cuts_per_round = true;
            argc--;
        }
        if(argv[argc-1][1]=='l') {
            multi_cuts_per_round = true;
            multi_cuts_limited   = true;
            argc--;
        }
        if(argv[argc-1][1]=='r') {
            seed = time(NULL);                              //randomization
            clog<<"I will randomize run using seed="<<seed<<endl;
            argc--;
        }
        if(argv[argc-1][1]=='g') {
            gamma = atoi(argv[argc-1]+2);
            clog<<"I will use gamma=="<<gamma<<endl;
            argc--;
        }
        if(argv[argc-1][1]=='i') {
            startsol = argv[argc-1]+2;
            clog<<"Input start sol: "<<startsol<<endl;
            argc--;
        }
        if(argc_start==argc){
            cerr<<"Failed parsing option '"<<argv[argc-1]
                <<"'. Run ./main to see all options."<<endl;
            exit(EXIT_FAILURE);
        }
    }
    if(!verbose)
        clog.setstate(ios_base::failbit);
    if(argc!=3){
        cerr<<"\nUsage: ./main method(std, prj or wEXTENSION) instance "
              "[-m[ultiple_cuts_per_round]] [-l[limited multi cuts] "
              "[-v[erbose]] [-gGAMMA, eg -g50] [-r[nd]] [-iINPUTSTARTSOL] \n"
              "           std=standard meth, prj=projective meth, "
              "w=write feasib sol to instance.EXTENSION\n"
              "           Use -v to enable printing log messages\n"
              "           Use -r[nd] to randomize algo \n"
              "           Best methods/switches: -m for prj and -l for std\n";
        exit(EXIT_FAILURE);
    }
}

void buildFeasibSol(char* outfile)
{
    double** rrows = new double*[::m];
    for(int i=0;i<(::m);i++){
        rrows[i] = new double[::n+2];
        for(int j=0;j<(::n)+2;j++)
            rrows[i][j] = rows[i][j];
    }
    //srand(9);
    for(int i=0;i<(::m);i++)
         if(rows[i][n+1]==LESS_THAN_EQ){
             for(int j=0;j<n;j++)
                //if(rand()%10>3)
                    rrows[i][j]+=0.02*absVal(rrows[i][j]);
             rrows[i][n]-=0.00150;
         }
    //double tt=2;
    CuttingPlanesEngine cutPlanes(::n,sub_problem_single_cut); 
    for(int i=0; i < (::n); i++){
        assert(lb[i]<=ub[i]);
        cutPlanes.setVarBounds(i,lb[i],ub[i]);
    }
    for(int i=0; i < (::m); i++){
        cutPlanes.modelAddWithSense(rrows[i],rrows[i][::n],rrows[i][::n+1]);
        delete[] rrows[i];
    }
    delete[] rrows;
    //for(int i=0;i<n;i++) obj[i]=0;
    //for(int i=0;i<n;i++) obj[i]=rand()%20;
    cutPlanes.setObjCoefsMinimize(obj);

    nominalObj = cutPlanes.solve();
    if(nominalObj!=INT_MAX)
        cout<<"I generated a feasible solution."<<endl;
    else{
        cout<<"I can not generate feasible solution. Exit.\n";
        exit(EXIT_FAILURE);
    }

    xbase = new double[n]; 
    cutPlanes.getPrimals(xbase); 
    cout<<"Loaded solution to xbase.\n";
    if(outfile!=NULL){
        cout<<" Writing it to "<<outfile<<"  ..   ";
        ofstream out(outfile);
        for(int i=0;i<n;i++)
            out<<setprecision(21)<<xbase[i]<<" ";
        out.close();
        cout<<"Done"<<endl;
    }
}


int main(int argc, char** argv)
{
    checkAllParams(argc,argv);
    runStd = 1;
    if(strcmp(argv[1],"std"))                             //if first arg different from std
        runStd = 0;                                       //run projective
    readInstance(argv[2]);

    /*-----------          GENERATE FEASIBLE SOLUTION MODE        -----------*/
    if(argv[1][0]=='w'){
        char* outsol = new char[strlen(argv[2])+10];
        strcpy(outsol,argv[2]);                          //instance
        strcat(outsol,".");
        if(strlen(argv[1])==1)
            strcat(outsol,"start");
        else 
            strcat(outsol,argv[1]+1);
        buildFeasibSol(outsol);
        return EXIT_SUCCESS;
    }

    /*-----------        INIT CUT PLANES OBJ and CONSTRAINTS       -----------*/
    CuttingPlanesEngine cutPlanes(::n,sub_problem,m*2); //m*2=maximum m*2 cuts per iter
                                                        //two cuts per row maxi,
                                                        //useful for projection
    //cutPlanes.activateLog();
    for(int i=0; i < (::n); i++){
        assert(lb[i]<=ub[i]);
        cutPlanes.setVarBounds(i,lb[i],ub[i]);
    }
    if(seed>0)
        cutPlanes.setObjCoefsMinRandomizedSolving(obj,seed);
    else
        cutPlanes.setObjCoefsMinimize(obj);
    for(int i=0; i < (::m); i++)
        cutPlanes.modelAddWithSense(rows[i],rows[i][::n],rows[i][::n+1]);

    /*-----------------   INPUT SOL FOR PROJECTIVE METHOD   ------------------*/
    if(!runStd){
        if(startsol==NULL){
            startsol = new char[strlen(argv[2])+10];
            strcpy(startsol,argv[2]);
            strcat(startsol,".start");
        }
        ifstream in(startsol);
        if(in.good()){
            xbase = new double[::n];
            d     = new double[::n];
            for(int i=0;i<(::n);i++)
                in>>xbase[i];
            clog<<"I loaded start solution from "<<startsol<<endl;
        }else{
            d     = new double[::n];
            cerr<<"Can not open input feasible solution '"<<startsol
                <<"'. Use 'wEXTENSION' instead of "<<argv[1]<<" to generate a "
                  "solution and write it to "<<argv[2]<<".EXTENSION\n"
                <<"I will try to generate a solution now!"<<endl;
            buildFeasibSol(NULL);
        }
    }

    /*-----------------           RUN CUT PLANES           ------------------*/
    nominalObj = cutPlanes.solve();
    //int cuts_start = cutPlanes.getNbCuts();
    if(cutPlanes.runCutPlanes(100000, 150000, iters, tmCPlanes)==EXIT_FAILURE){
        finalObj = cutPlanes.getObjVal();
        cout<<"\nFinal obj val (with timeout or error)="<<setprecision(11)
            <<finalObj <<" obtained after "<<iters<<" iters of CutPlanes and "
            <<tmCPlanes<<" secs.\n";
        cout<<"Infeasible\n";
        return EXIT_FAILURE;
    }

    /*-----------------           REPORT OPT SOL           ------------------*/
    double*opt = new double[n]; 
    cutPlanes.getPrimals(opt); 
    //for(int i=0;i<n;i++){
    //    clog<<opt[i]<<" ";
    //}
    //clog<<endl;
    
    finalObj = cutPlanes.getObjVal();
    cout<<"\nFinal obj val ="<<setprecision(12)<<finalObj<<" obtained after "
        <<iters <<" iters of CutPlanes and " <<tmCPlanes<<" secs.\n";
    cout<<"Tabular data below for latex inclusion:\nRatio:"<<setprecision(4)<<setw(7)<<(finalObj-nominalObj)/absVal(nominalObj)*100;
    if(iterLowGap!=-1)
        cout<<" "<<setw(9)<<iterLowGap<<" "<<setw(9)<<tmLowGap;
    cout<<" ITERS"<<setw(5)<<iters<<" TIME"<<setw(9)<<tmCPlanes;
    //if(iterLowGap!=-1)
    //    cout<<"P"<<setw(6)<<setprecision(3)<<100.0*tmLowGap/tmCPlanes;
    if(total_multi_cuts>0)
        cout<<" MULTICUTS "<<setw(5)<<total_multi_cuts;   //or cutPlanes.getNbCuts()-cuts_start;
    //cout<<"| "<<setw(6)<<setprecision(3)<<100.0*tmSort/tmCPlanes;
    cout<<endl;
    //fast exit
    exit(EXIT_SUCCESS);
}

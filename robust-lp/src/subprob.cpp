/*---------------+---------------------------------------------------------------+-------------+
                 | Author: @ Daniel Porumbel 2021                                |
                 |License: Any person obtaining a copy of this code is free to   |
                 |         use it in any manner, subject to two conditions:      |            
                 |           1) no profit may ever be made from using this code  |
                 |           2) these 5 lines of text shall be included          |
                 +--------------------------------------------------------------*/
#include "general.h"
#include "subprob.h"
#include "CuttingPlanesEngine.h"
#include <algorithm>
#include <cassert>
#include <climits>
#include <iostream>
#include <iomanip>
#include <set>
using namespace std;

#define EPSILON    1.0E-6
#define DEV_NOM    0.01      //it can deviate 1% of the nominal value

double*  slack = NULL;           
double*  maxDev  = NULL;    //maximum deviation in a robust cut using gamma terms
int*     order = NULL;
double*  profits;           //sortProfits below sorts putting max profits first
double*  toCut;
int      iter = 0;
extern double tmSort;
extern int total_multi_cuts;

/*-----------------------+------------------------------------+--------------------------
                         |VARIOUS GENERAL ROUTINES (eg., sort)|
                         +-----------------------------------*/
bool comparator (int i,int j) { return (profits[i]>profits[j]); }
struct classcomp {
  bool operator() (const int& lhs, const int& rhs) const
  {return profits[lhs]>profits[rhs];}
};
void sortProfits(int* order)
{
    double tmStart = getCPUTime();
    multiset<int,classcomp> bestGamma;
    int added = 0;
    for(int i=0;i<n;i++)
        if(profits[i]>0){
            bestGamma.insert(i);
            added++;
            if(added>gamma)
                bestGamma.erase(--bestGamma.end());
        }
    int ii = 0;
    int* selected = new int[n];
    for(int i=0;i<n;i++)
        selected[i] = 0;
    for (multiset<int>::iterator it=bestGamma.begin(); it!=bestGamma.end(); ++it){
        assert(ii>=0);
        assert(ii<n);
        order[ii] = *it;
        selected[*it] = 1;
        ii++;
    }
    int jj= 0;
    for(;ii<n;ii++){
        while(selected[jj])
            jj++;
        order[ii] = jj;
        jj++;
    }
    delete[] selected;
    //for(int i=0;i<(::n);i++)
    //    order[i]=i;
    //std::sort(order, order+::n,comparator);                                
    tmSort+= (getCPUTime()-tmStart);
}

double absVal(double z){
    if(z>=0)
        return z;
    return -z;
}
double scalprod(double *x, double* y)
{
    double prod=0;
    for(int i=0;i<(::n);i++)
        prod += x[i] * y[i];
    return prod;
}


/*-----------------------+------------------------------------+--------------------------
                         | Routines for projection algorithms |
                         +-----------------------------------*/
//Given (implicit) parameters xbase and d, solve proj-subprob(xbase->td) with
//regards to row i only, fill resulting constraint in newRow and return tStar
//If xbase+td can not be separated, return tStar=t; 
double tStarForRow(int i, double t, double* newRow)
{
    if(profits==NULL){
        profits = new double[::n];
        order   = new int[::n];
        toCut = new double[n];
    }
    for(int j=0;j<(::n);j++)
        toCut[j] = xbase[j] + t*d[j];
    for(int j=0;j<(::n);j++)
        profits[j] = absVal(rows[i][j]*toCut[j]);
    sortProfits(order);
    for(int j=0,jj=order[j]; j<(::gamma); ++j,jj=order[j])
        if(rows[i][jj]*toCut[jj]>=0)
            newRow[jj] = rows[i][jj] * (1+DEV_NOM);
        else
            newRow[jj] = rows[i][jj] * (1-DEV_NOM);
    for(int j=::gamma,jj=order[j];j<(::n);++j,jj=order[j])
        newRow[jj] = rows[i][jj];

    double base    = scalprod(newRow,xbase);
    double advance = scalprod(newRow,d) * t;
    if(advance<EPSILON*t)                                    //null advance
        return t;

    //Below code can check if xbase is really feasible. Numerical probls not fatal
    //if((base>rows[i][::n])&&(base<rows[i][::n]+EPSILON*100))
        //base = rows[i][::n];
    //if(base> rows[i][::n]) {
    //    clog<<"diff at "<<i<<"="<<base-rows[i][::n]<<"base="<<base<<","<<rows[i][::n]<<endl;
    //    for(int i=0;i<n;i++)
    //        if(newRow[i]*xbase[i]!=0)
    //            clog<<" "<<newRow[i]<<"/"<<xbase[i]<<" ";
    //}
    //assert(base <= rows[i][::n]);
    if( base+advance <= rows[i][::n] + EPSILON)              //can not cut toCut
        return t;
    double newtStar = t*((rows[i][::n]-base)/advance);
    return newtStar;
}

//linux coding style: "macros resembling functions may be named in lower case." 
#define checkBoundsAndPrint(xbase,d,tStar)                                    \
  do{                                                                         \
    double ubNow   = scalprod(obj,xbase)+tStar*scalprod(obj,d);               \
    double lbNow   = ::lowerBound; /*from CuttingPlanesEngine.cpp*/           \
/*  double ubDelta = ubNow - nominalObj;                                      \
    double lbDelta = lbNow - nominalObj;                                   */ \
    static double bestUb = INT_MAX;                                           \
    if(ubNow<bestUb)                                                          \
        bestUb = ubNow;                                                       \
/*  clog<<"      (iter, ub) (iter, lb) below, with ub or lb expressed ";      \
    clog<<"in terms of bndVal/nomVal*1000 or similar if negatives:\n";        \
    clog<<"      ("<<setw(4)<<iter<<","<<setw(6)<<setprecision(5)             \
                   <<1000.0*ubDelta/absVal(nominalObj)+1000<<")   ("          \
                   <<setw(4)<<iter<<","<<setw(6)<<setprecision(5)             \
                   <<1000.0*lbDelta/absVal(nominalObj)+1000 <<") "<<endl;  */ \
    if(iterLowGap==-1)                                                        \
        if(bestUb-lbNow<0.01*min(absVal(lbNow),absVal(bestUb)) ){             \
            /*cout<<bestUb<<endl<<absVal(lbNow)<<endl<<iter<<endl;*/          \
            iterLowGap = iter;                                                \
            tmLowGap   = getCPUTime() - startCpuTime;                         \
        }                                                                     \
  }while(false);

//return -1 if violated cut found; the cut that will be added is: newRow^T x <= rHand
double projection (double*x, double * newRow,double&rHand)
{

    //1. VARIABLES AND ITERATION/TIME COUNTERS
    static double startCpuTime  = getCPUTime();
    iter++;
    clog<<"Mono-cut projection at iteration "<<iter<<":\n";
    for(int i=0;i<(::n);i++)
        d[i] = x[i]-xbase[i];
    double tStar  = 1;
    double tStarLast = 1;
    int    bstRow = 0;

    //if(iter%100==0){
    //    char* outfile = new char[200];
    //    sprintf(outfile, "stocfor3.txt.%d.start",iter);
    //    ofstream out(outfile);
    //    for(int i=0;i<n;i++)
    //        out<<setprecision(21)<<xbase[i]<<" ";
    //    out.close();
    //}
    //for(int i=0;i<(::m);i++)
    //    if(rows[i][n+1]==LESS_THAN_EQ)
    //    if(tStarForRow(i,tStar,newRow, newRowProd_d)<0.00005)
    //        clog<<"rrows["<<i<<"][n]-=tt;"<<endl;
    //exit(1);

    //2. FIND tStar ANALYZING ROW BY ROW
    for(int i=0;i<(::m);i++)
        if(rows[i][n+1]==LESS_THAN_EQ){
            double tStarNew = tStarForRow(i,tStar,newRow); 
            while(tStarNew<tStar){
                bstRow = i;
                tStarLast = tStar;
                tStar = tStarNew;
                clog<<"        better tStar using (2.4) "<<setprecision(11)<<tStar<<endl;
                    //<<"for slack="<<rows[i][::n]-scalprod(newRow,xbase)<<"; ";
                tStarNew = tStarForRow(i,tStar,newRow); 
            }
            if(tStar<=EPSILON){
                tStar = 0;
                break;
            }
        }

    //3. CHECK BOUNDS AND STOP IF OPTIMAL
    checkBoundsAndPrint(xbase,d,tStar);

    if(tStar==1){                  //x is feasible
        clog<<"        Found tStar=1, the outer solution is feasible"<<endl;
        return 1;
    }
    
    //4. DETERMINE A CUT TO BE RETURNED
    //This is not used to find the next interior point. The actual cut to return 
    //only needs to separate x+d.  //You can not find such a cut by separating x+tStar d 
    //because x+tStar d is not separable and can lead to a cut not separating x+d
    //tStarLast below can eliminate numerical problems associated to tStarLast=~=tStar
    tStarLast = tStar*0.5+0.5;
    tStarForRow(bstRow,tStarLast,newRow);
    rHand = rows[bstRow][::n];

    #ifndef NDEBUG
    clog<<"rHnd"<<rHand<<"new*x"<<scalprod(newRow,x)<<" "<<"new*base"<<scalprod(newRow,x)<<endl<<endl;
    assert(scalprod(newRow,x)>rHand);
    #endif

    //Code to print solutions and investigate the bang bang effects. Also put 1
    //instead of 0.1 below.
    //static int printed=0;
    //    if(printed%10==0){
    //    for(int i=0;i<15;i++)
    //        cout<<setw(5)<<left<<setprecision((xbase[i]>1000)?4:(xbase[i]>1?4:3))<<xbase[i]<<" ";
    //    cout<<"TT"<<endl;
    //}
    //printed++; if(printed==100) exit(1);

    //5. UPDATE XBASE (INTERIOR POINT)
    for(int i=0;i<n;i++){
        xbase[i] = xbase[i] + 0.1*tStar*d[i];
        if(xbase[i]<lb[i])
            xbase[i] = lb[i];
    }

    clog<<"        I return tStar="<<tStar<<endl;
    return -1;
}

/*-----------------------+------------------------------------+--------------------------
                         |  Routines for standard separation  |
                         +-----------------------------------*/
void calcAllSlacks(double* sol)
{
    if(slack==NULL)
        slack = new double[::m];
    for(int i=0;i<(::m);i++)
        switch((int)rows[i][n+1]){
            case EQUALITY:                    //for equalities, do 
                slack[i] = -1;                //not care about slacks
                break;
            case LESS_THAN_EQ:
                slack[i] = rows[i][n] - scalprod(sol,rows[i]);
                if((slack[i]<=0)&&(slack[i]>=-EPSILON))
                    slack[i] = 0;
                if(slack[i]<0){
                    cerr<<setprecision(12)<<"Warning: the nominal constraint i="<<i
                        <<"quite far from satisfied. slack="<<slack[i]<<",rhand="<<rows[i][n]<<endl;
                    if(slack[i]>=-EPSILON*50)
                        slack[i] = 0;
                }
                assert(slack[i]>=0);
                break;
            default:
                cerr<<"I only handle == or <=. Change code below for more\n";
                cerr<<"Attention: I mean a^top x <=c which is written as "
                      " c >= a^\top x in the input file\n";
                exit(EXIT_FAILURE);
        }
}

void calcMaxDeviations(double *sol)
{
    if(maxDev == NULL){
        maxDev    = new double[::m];
        profits = new double[::n];
        order   = new int[::n];
    }
    for(int i=0;i<(::m);i++)
        if(rows[i][::n+1]==LESS_THAN_EQ){
            for(int j=0;j<(::n);j++)
               profits[j] = absVal(rows[i][j]*sol[j]); 
            sortProfits(order);
            maxDev[i] = 0;
            for(int j=0;j<(::gamma);j++)
                maxDev[i] += profits[order[j]];
            maxDev[i] = maxDev[i]*DEV_NOM;
        }
}

double separation (double*x, double * newRow,double&rHand)
{
    static int iter = 0;
    iter++;
    calcAllSlacks(x);
    calcMaxDeviations(x);

    int bstRow = 0;
    while(slack[bstRow]<0){   //negative slack means equality constraint
        bstRow++;             //so that slack existence does not make sense
        assert(bstRow<(::m));
    }
    for(int i=bstRow+1;i<(::m);i++)
      if(rows[i][::n+1]==LESS_THAN_EQ)
        if(maxDev[i]-slack[i]>=maxDev[bstRow]-slack[bstRow])
            bstRow = i;
    
    //print bounds, useful for plotting the running profile
    //clog<<"("<<setw(4)<<iter<<","<<setw(6)<<setprecision(6)              
    //               <<"current opt val/nominal val times 1000:"
    //               <<1000.0*(::lowerBound-nominalObj)/absVal(nominalObj)+1000<<")"<<endl;


    if(maxDev[bstRow]<=slack[bstRow])
        return 1;                                         //opt sol

    for(int j=0;j<(::n);j++)
        profits[j] = absVal(rows[bstRow][j]*x[j]); 
    sortProfits(order);

    if(rows[bstRow][::n+1]==LESS_THAN_EQ){               //adding <= ineq
        rHand = rows[bstRow][::n];
        for(int j=0,jj=order[j];j<(::gamma);++j,jj=order[j])
            if(rows[bstRow][jj]*x[jj]>=0)
                newRow[jj] = rows[bstRow][jj] * (1+DEV_NOM);
            else
                newRow[jj] = rows[bstRow][jj] * (1-DEV_NOM);
        for(int j=::gamma,jj=order[j];j<(::n);++j,jj=order[j])
            newRow[jj] = rows[bstRow][jj];
    }
    //clog<<"max violation:"<<setw(10)<<slack[bstRow] - maxDev[bstRow];//<<"/profmax"<<setw(9)<<profits[order[0]];
    //for(int j=0,jj=order[j];j<(::gamma);++j,jj=order[j])
    //    clog<<setw(5)<<jj;
    return slack[bstRow] - maxDev[bstRow];
}

/*-----------------------+------------------------------------+--------------------------
                         |  Routines for mulit-cuts per round |
                         +-----------------------------------*/
int compare_descending (const void * a, const void * b)
{
    double cmp = - ( *(double*)a - *(double*)b );
    if(cmp<0)
        return -1;
    if(cmp>0)
        return 1;
    return 0;
}

double separation_multi (double*x, double * newRow,double&rHand,
                         double**newRows, double*newRHands, int& newMore,
                         bool multi_cuts_limited)
{
    static int iter = 0;
    iter++;
    clog<<"Multi-cut separation iteration "<<iter<<":";
    calcAllSlacks(x);
    calcMaxDeviations(x);

    int bstRow = 0;
    while(slack[bstRow]<0){   //negative slack means equality constraint
        bstRow++;             //so that slack existence does not make sense
        assert(bstRow<(::m));
    }
    for(int i=bstRow+1;i<(::m);i++)
        if(rows[i][::n+1]==LESS_THAN_EQ)
            if(maxDev[i]-slack[i]>=maxDev[bstRow]-slack[bstRow])
                bstRow = i;
    
    //print bounds, useful for plotting the running profile
    //clog<<"("<<setw(4)<<iter<<","<<setw(6)<<setprecision(6)              
    //               <<"current opt val/nominal val times 1000:"
    //               <<1000.0*(::lowerBound-nominalObj)/absVal(nominalObj)+1000<<")"<<endl;


    if(maxDev[bstRow]<=slack[bstRow])
        return 1;                                         //opt sol

    for(int j=0;j<(::n);j++)
        profits[j] = absVal(rows[bstRow][j]*x[j]); 
    sortProfits(order);

    if(rows[bstRow][::n+1]==LESS_THAN_EQ){               //adding <= ineq
        rHand = rows[bstRow][::n];
        for(int j=0,jj=order[j];j<(::gamma);++j,jj=order[j])
            if(rows[bstRow][jj]*x[jj]>=0)
                newRow[jj] = rows[bstRow][jj] * (1+DEV_NOM);
            else
                newRow[jj] = rows[bstRow][jj] * (1-DEV_NOM);
        for(int j=::gamma,jj=order[j];j<(::n);++j,jj=order[j])
            newRow[jj] = rows[bstRow][jj];
    }
    //clog<<"max violation:"<<setw(10)<<slack[bstRow] - maxDev[bstRow];//<<"/profmax"<<setw(9)<<profits[order[0]];
    //for(int j=0,jj=order[j];j<(::gamma);++j,jj=order[j])
    //    clog<<setw(5)<<jj;


    //Now I add more multiple cuts
    double min_excess = 0;
    if(multi_cuts_limited){
        double* excess = new double[::m];
        for(int ii=0;ii<(::m);ii++){
            if(rows[ii][::n+1]==LESS_THAN_EQ)
                excess[ii] = maxDev[ii]-slack[ii];
            else
                excess[ii] = INT_MIN;
        }
        qsort(excess, ::m, sizeof(double), compare_descending);
        min_excess = excess[5];
        if(min_excess<0)
            min_excess=0;
        //cout<<"excess[0]="<<excess[0]<<","<<"excess[10]="<<excess[10]<<"="<<min_excess<<endl;
        delete[] excess;
        //cout<<excess[0]<<","<<excess[1]<<"...|"<<excess[10]<<"|,"<<excess[n-1]<<endl;
        //cout<<"min_excess="<<min_excess<<endl;
        //if(!(excess[0]>=excess[1])){
        //    cout<<"!(excess[0]>=excess[1])"<<endl;
        //    exit(1);
        //}
        //if(!excess[1]>=excess[10]){
        //    cout<<"!excess[1]>=excess[10]"<<endl;
        //    exit(1);
        //}
    }

    newMore = 0;           //newMore is also the current cut 
    int row = 0;
    while(slack[row]<0)    //negative slack means equality constraint
        row++;             //so that slack existence does not make sense
    
    for(;row<(::m);row++)
      if(row!=bstRow)
        if(rows[row][::n+1]==LESS_THAN_EQ)
            if(maxDev[row]-slack[row]>min_excess){
                for(int j=0;j<(::n);j++)
                    profits[j] = absVal(rows[row][j]*x[j]); 
                sortProfits(order);
                newRHands[newMore] = rows[row][::n];
                for(int j=0,jj=order[j];j<(::gamma);++j,jj=order[j])
                    if(rows[bstRow][jj]*x[jj]>=0)
                        newRows[newMore][jj] = rows[row][jj] * (1+DEV_NOM);
                    else
                        newRows[newMore][jj] = rows[row][jj] * (1-DEV_NOM);
                for(int j=::gamma,jj=order[j];j<(::n);++j,jj=order[j])
                    newRows[newMore][jj] = rows[row][jj];
                if(total_multi_cuts+(newMore+1)<10000)
                    newMore++;
            }
    clog<<"        multiple cuts added: "<<(newMore+1)<<endl;
    return slack[bstRow] - maxDev[bstRow];
}

//return -1 if violated cut found; the cut that will be added is: newRow^T x <= rHand
double projection_multi (double*x, double * newRow,double&rHand,
                         double**newRows, double*newRHands, int& newMore,
                         bool multi_cuts_limited)
{

    //1. VARIABLES AND ITERATION/TIME COUNTERS
    static double startCpuTime  = getCPUTime();
    iter++;
    clog<<"Multi-cut projection at iteration "<<iter<<":"<<endl;
    for(int i=0;i<(::n);i++)
        d[i] = x[i]-xbase[i];
    double tStar  = 1;
    double tStarLast = 1;
    int    bstRow = 0;

    //if(iter%100==0){
    //    char* outfile = new char[200];
    //    sprintf(outfile, "stocfor3.txt.%d.start",iter);
    //    ofstream out(outfile);
    //    for(int i=0;i<n;i++)
    //        out<<setprecision(21)<<xbase[i]<<" ";
    //    out.close();
    //}
    //for(int i=0;i<(::m);i++)
    //    if(rows[i][n+1]==LESS_THAN_EQ)
    //    if(tStarForRow(i,tStar,newRow, newRowProd_d)<0.00005)
    //        clog<<"rrows["<<i<<"][n]-=tt;"<<endl;
    //exit(1);

    //2. FIND tStar ANALYZING ROW BY ROW
    newMore = 0;
    for(int i=0;i<(::m);i++)
        if(rows[i][n+1]==LESS_THAN_EQ){
            double tStarNew = tStarForRow(i,tStar,newRows[newMore]); 
            newRHands[newMore] = rows[i][::n];
            if((tStarNew<tStar) || (!multi_cuts_limited) )
              if(total_multi_cuts+(newMore+1)<10000)
                newMore++;
            //if(newMore>5) newMore=5;
            while(tStarNew<tStar){
                bstRow = i;
                tStarLast = tStar;
                tStar = tStarNew;
                clog<<"        better tStar using eq. (2.4): "<<setprecision(8)<<tStar<<endl;
                    //<< "at i="<<i;
                    //<<"for slack="<<rows[i][::n]-scalprod(newRow,xbase)<<"; ";
                //tStarNew = tStarForRow(i,tStar,newRow); 
                tStarNew = tStarForRow(i,tStar,newRows[newMore]); 
                newRHands[newMore] = rows[i][::n];
                if((tStarNew<tStar) || (!multi_cuts_limited) )
                    if(total_multi_cuts+(newMore+1)<10000)
                        newMore++;
                //if(newMore>5) newMore=5;
            }
            if(tStar<=EPSILON){
                tStar = 0;
                break;
            }
        }

    //3. CHECK BOUNDS AND STOP IF OPTIMAL
    checkBoundsAndPrint(xbase,d,tStar);

    if(tStar==1){                  //x is feasible
        clog<<"        found tStar=1, i.e., outer solution is feasible"<<endl;
        return 1;
    }
    
    //4. DETERMINE A CUT TO BE RETURNED
    //This is not used to find the next interior point. The actual cut to return 
    //only needs to separate x+d.  //You can not find such a cut by separating x+tStar d 
    //because x+tStar d is not separable and can lead to a cut not separating x+d
    //tStarLast below can eliminate numerical problems associated to tStarLast=~=tStar

    //This cut will actually be returned twice, because it is already newRows
    tStarLast = tStar*0.5+0.5;
    tStarForRow(bstRow,tStarLast,newRow);
    rHand = rows[bstRow][::n];

    #ifndef NDEBUG
    clog<<"rHnd"<<rHand<<"new*x"<<scalprod(newRow,x)<<" "<<"new*base"<<scalprod(newRow,x)<<endl<<endl;
    assert(scalprod(newRow,x)>rHand);
    #endif

    //Code to print solutions and investigate the bang bang effects. Also put 1
    //instead of 0.1 below.
    //static int printed=0;
    //    if(printed%10==0){
    //    for(int i=0;i<15;i++)
    //        cout<<setw(5)<<left<<setprecision((xbase[i]>1000)?4:(xbase[i]>1?4:3))<<xbase[i]<<" ";
    //    cout<<"TT"<<endl;
    //}
    //printed++; if(printed==100) exit(1);

    //5. UPDATE XBASE (INTERIOR POINT)
    for(int i=0;i<n;i++){
        xbase[i] = xbase[i] + 0.1*tStar*d[i];
        if(xbase[i]<lb[i])
            xbase[i] = lb[i];
    }
    clog<<"        I return tStar="<<tStar<<" Multiple cuts added: "<<(newMore+1)<<endl;
    //total_multi_cuts+=(newMore+1);
    return -1;
}




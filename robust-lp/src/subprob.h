/*See accompanying .cpp file for licence info     */
#ifndef SUBPROB_H
#define SUBPROB_H

#define EQUALITY     0
#define LESS_THAN_EQ 1

extern int      n;               //number of variables
extern int      m;               //number of initial rows
extern double*  obj;
extern double*  lb;
extern double*  xbase;
extern double*  d;         //used to solve proj_subprob(xbase->d)
extern double** rows;            //rows[...][n] = rHand, rows[...][n+1]=0/-1/1 ~ ==/>=/<=
extern int      gamma;
extern double   nominalObj;
extern int      iterLowGap;//iteration when ub<=bstLowerBound*1.2
extern double   tmLowGap  ;//tm for above

double absVal(double);
double scalprod(double *x, double* y);
double separation (double*x, double * newRow,double&newRHand);
double projection (double*x, double * newRow,double&newRHand);

double separation_multi (double*x, double * newRow,double&rHand,
                         double**newRows, double*newRHands, int& newMore,
                         bool multi_cuts_limited//if true, add maximum 10 cuts per iter
                         );
double projection_multi (double*x, double * newRow,double&newRHand,
                         double**newRows, double*newRHands, int& newMore,
                         bool multi_cuts_limited//if true, add only cuts that decrease tStar
                         );

#endif

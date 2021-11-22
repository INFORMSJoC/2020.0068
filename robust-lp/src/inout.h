/* see accompanying .cpp file for licence info  */
#ifndef INOUT_H
#define INOUT_H
extern int n;       //number of variables
extern int m;       //number of initial rows
extern double*  obj;
extern double*  lb;
extern double*  ub;
extern double** rows;
void readInstance(char* filename);
#endif

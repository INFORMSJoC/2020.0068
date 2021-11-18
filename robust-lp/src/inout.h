/*---------------+---------------------------------------------------------------+--------------+
                 | Author: Daniel Porumbel 2018       daniel.porumbel@cnam.fr    |
                 |License: Creative Commons Attribution 3.0 Unported License     |
                 |         http://creativecommons.org/licenses/by/3.0/           |             
                 +--------------------------------------------------------------*/
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

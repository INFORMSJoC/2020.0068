/*----------------------------------------------------------------------------------------------+
|               Input and output functions, they can read the n^th instance from file           |
+--------+-------------------------------------------------------------------------+------------+
         | See file LICENSE at the root of the git project for licence information |
         +------------------------------------------------------------------------*/

#ifndef INOUT_H_INCLUDED
#define INOUT_H_INCLUDED

extern double  C;                    //capacity
extern int     n;                    //the number of dual variables
extern int*    b;                    //demands
extern int*    w;                    //weights
extern int     currInst;             //instance nr (zero indexed) to solve

//Returns the instance number zeroIndexedNumber from file. Exits if error.
void readInstNrFromFile(int zeroIndexedNumber, char* file);

//Returns an integer parameter from the config file, e.g., the value of kMax
//Default return value: INT_MAX (when requestedParam is not found)
int getParamFromConfigFile(char* requestedParam);
#endif

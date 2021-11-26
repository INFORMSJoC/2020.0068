/*-------+-------------------------------------------------------------------------+------------+
         | See file LICENSE at the root of the git project for licence information |
         +------------------------------------------------------------------------*/

#include "inout.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
using namespace std;
void readInstance(char* filename)
{
    ifstream inFile (filename);
    if(!inFile.good()){
        cerr<<"I can not read instance from "<<filename<<endl;
        exit(EXIT_FAILURE);
    }
    string   line ;
    int      lineNo = 0 ;
    string   dummy;
    int      var_id;
    double   coef ;
    while(getline(inFile,line)){
        istringstream in (line);
        switch(lineNo){
            case(0):
                in>>::n;
                break;
            case(1):
                in>>::m;
                ::obj         = new double[::n];
                ::lb          = new double[::n];
                ::ub          = new double[::n];
                ::rows = new double*[::m];
                for(int i=0;i<(::m);i++)
                    ::rows[i] = new double[::n+2];
                for(int i=0;i<(::n);i++){
                    ::obj[i] = 0;
                    for(int j=0;j<(::m);j++)
                        ::rows[j][i] = 0;
                }
                break;
            case 2:
                in>>dummy;
                if(dummy=="LBVARS")
                    in>>lb[0];
                else
                    istringstream(dummy)>>lb[0];
                for(int i=1;i<(::n);i++)
                    in >> lb[i];
                break;
            case 3:
                in>>dummy;
                if(dummy=="UBVARS")
                    in>>ub[0];
                else
                    istringstream(dummy)>>ub[0];
                for(int i=1;i<(::n);i++)
                    in >> ub[i];
                break;
            case(4):
                in>>dummy;
                while(in>>var_id){
                    in>>coef;
                    obj[var_id] = coef;
                    assert(var_id<(::n));
                }
                break;
            default:
                in>>dummy;
                in>>rows[lineNo-5][::n];          //1st in line: right hand side
                in>>dummy;
                rows[lineNo-5][::n+1] = 0;        //0  means ...==rHand
                if((dummy=="G")||(dummy=="<="))
                    rows[lineNo-5][::n+1] = -1;   //-1  means ...>=rhand
                if((dummy=="L")||(dummy==">="))
                    rows[lineNo-5][::n+1] = 1;    //1 means ...<=rhand
                
                while(in>>var_id){
                    in>>coef;
                    rows[lineNo-5][var_id] = coef;
                }
                break;
        }
        lineNo++;
    }
}

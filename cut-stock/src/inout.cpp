/*----------------------------------------------------------------------------------------------+
|               Input and output functions, they can read the n^th instance from file           |
+ -------+-------------------------------------------------------------------------+------------+
         | See file LICENSE at the root of the git project for licence information |
         +------------------------------------------------------------------------*/

#include "inout.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <climits>
#include <cstdlib>
using namespace std;

int     currInst;                       //the current instance

int     instancesInFile;                //the total number of instances
int     instanceType;
fstream currStream;                     //input stream

//1 = .1bp files
//2 = classical format with no extension
#define BP1_FILE 1
#define CLASSICAL 2
#define TXT_FILE 3

string instNameFromFile;   

//This function returns the file format of the cspFile: BP1_FILE or CLASSICAL
int inputFileFormat(char*cspFile)
{
    int pointPos = strlen(cspFile)-4;
    if(cspFile[pointPos]!='.'||//no extension
      (cspFile[pointPos+1]=='d'&&cspFile[pointPos+2]=='a'&&cspFile[pointPos+3]=='t')//.dat
    ){
            clog<<"Reading a classical file (no 3-character extension):"<<cspFile<<".\n";
            return CLASSICAL;
    }
    if(cspFile[pointPos+1]=='1'&&cspFile[pointPos+2]=='b'&&cspFile[pointPos+3]=='p'){//.1bp
            clog<<"Reading a 1bp file:"<<cspFile<<".\n";
            return BP1_FILE;
    }
    if(cspFile[pointPos+1]=='t'&&cspFile[pointPos+2]=='x'&&cspFile[pointPos+3]=='t'){//.1bp
            clog<<"Reading a txt file:"<<cspFile<<".\n";
            return TXT_FILE;
    }
    cerr<<"Unknown input file format (extension '"<<cspFile+pointPos+1<<"') not recognized\n";
    exit(1);
        
}
void readNextInst()
{
    int dump;
    if(instanceType==BP1_FILE){
        currStream >> ::C;
        currStream >> dump;
        currStream >> ::n;
    }else{
        currStream >> ::n;
        currStream >> ::C;
    }
    clog<<"Reading instance nr "<<::currInst;
    clog<<" with n="<<n<<" and tot. cap="<<C<<"and";

    w = new int[::n];
    b = new int[::n];

    double sum=0;
    for(int i=0;i<n;i++){ 
        currStream >> w[i];
        if(instanceType==TXT_FILE){
            b[i] = 1;
        }else{
            currStream >> b[i];
        }
        sum+=((double)w[i]) * b[i];
    }
    clog<<"sum wi*bi="<<sum<<endl;
    instNameFromFile = "";           //signals the fact that a new instance can be read   
}


//returns 0 when there are no more instances
//the instance name is stored in instNameFromFile
int readInstNameFromFile()
{
    while(instNameFromFile.size()==0&&currStream.good()){
        getline(currStream,instNameFromFile);
    } 
    return currStream.good();
}
int moreInst()
{
    if( (instanceType==BP1_FILE) || (instanceType==TXT_FILE) )
        return ::currInst<=instancesInFile;
    else{
        return readInstNameFromFile();
    }
}


void readInstNrFromFile(int insNrRequest, char* filename)
{
    ::currInst = 0;                     //first instance to be loaded 
    currStream.open(filename, fstream::in);  
    instanceType = inputFileFormat(filename);
    if( (instanceType==BP1_FILE) || (instanceType==TXT_FILE) )
        currStream >> instancesInFile;  //number of instances
    else{
        instancesInFile = -1;           //unknown number of instances, but the read stops when
//        readInstNameFromFile();         //there are no other instances names in the file
    }

    while ( moreInst() && (::currInst<=insNrRequest) ){
        readNextInst();
        ::currInst++;
        if(::currInst<insNrRequest){
            delete[] ::w;
            delete[] ::b;
        }
    }     
    if(::currInst-1<insNrRequest){
        cerr<<"Last read instance nr="<<::currInst-1<< "while you asked inst nr"<<insNrRequest<<endl;
        exit(EXIT_FAILURE);
    }
}
//Returns an integer parameter from the config file, e.g., the value of kMax
//Default return value: INT_MAX (when requestedParam is not found)
int getParamFromConfigFile(char* requestedParam)
{
    if(rindex(requestedParam,'/')!=NULL)
        requestedParam = rindex(requestedParam,'/')+1;
    ifstream fConfig;
    char oneWord[100];
    char configFile[100] = "init.cfg";
    fConfig.open(configFile);
    string s;
    if(fConfig.good()) {
        while(!fConfig.eof()) {
            fConfig >> oneWord;
            if(!strcmp(requestedParam,oneWord)) {
                fConfig >> oneWord;
                int toRet;
                sscanf(oneWord,"%d",&toRet);
                fConfig.close();
                cout<<requestedParam<<" set to "<<toRet<<" using the value from "<<configFile<<endl;
                return toRet;
            }       
        }
        fConfig.close();        
    }   
    cerr<<"You can insert in (in "<<configFile<<") a line with the form: '"<<requestedParam<<" <VALUE>' in case you want to specify a value for this parameter."<<endl;
    return INT_MAX;
    
}


#   +---------------------------------------------------------------+
#   | Author: @ Daniel Porumbel 2021                                |
#   |License: Any person obtaining a copy of this code is free to   |
#   |         use it in any manner, subject to two conditions:      |            
#   |           1) no profit may ever be made from using this code  |
#   |           2) these 5 lines of text shall be included          |
#   +--------------------------------------------------------------*/


#Code below searches for cplex install dir. If you know it, you can erase all
#lines below (up to first comment) and just put a line as follows.
#CPLEX=your/cplex/folder/
#Otherwise let it check for several folders and it will set the first existing
#one. If all fail, the compilation will halt before trying to build main
ifneq ($(wildcard /data/cplex126/),)  #check if this folder exists
    CPLEX   = /data/cplex126/
    CPLEXVER="cplex126"
endif 
ifneq ($(wildcard /opt/ibm/ILOG/),)            
    #Cplex is often installed in a subdir of above folder.
    #take last subdir inside (ex, CPLEX_Studio...) as it can be the last version
    CPLEXVER= $(shell ls /opt/ibm/ILOG/|tail -1)
    CPLEX   = /opt/ibm/ILOG/${CPLEXVER}
endif 

CPLEXVERNO_ALLDIGITS=$(shell echo ${CPLEXVER}|sed "s/[^0-9]//g")  #remove non digits
CPLEXVERNO=$(shell echo ${CPLEXVERNO_ALLDIGITS}|cut -c 1-3)       #take first 3 letters (eg. 126)

#--------------------------------------------------------------------------------------------------------------------------------------------
#The rest depends on the cplex install dire determined above

CPLEXDIR      = $(CPLEX)/cplex/
CONCERTDIR    = $(CPLEX)/concert/

#--------------------------------------------------------------------------------------------------------------------------------------------
# Compiler selection, code optimization, debug and warning options
# ---------------------------------------------------------------------
#use DNDEBUG to remove asserts and turn off some clog debug messages
CCC = g++ 
CCFLAGS = -g -DNDEBUG -O3 -flto -m64 -Wextra -fomit-frame-pointer -funroll-loops 

#--------------------------------------------------------------------------------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

LIBFORMAT  = static_pic

#cplex related
#use SYSTEM=x86-64_sles10_4.1 and remove -ldl for cplex versions <12.6
#SYSTEM     = x86-64_linux

SYSTEM=$(shell ls $(CPLEXDIR)/lib/|grep x86)

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CCLNFLAGS =" " 
CCLNFLAGSCPLEX = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lrt -lpthread -ldl
               #-m32 -pthread -lm (add more when needed)

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGSCPLEX = $(CCFLAGS) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)  -DIL_STD #DIL_STD: CPLEX specific macro

#--------------------------------------------------------------------------------------------------------------------------------------------
# Main compilation and linking commands
# ---------------------------------------------------------------------

all:cplexdirexists main
.SILENT: cplexdirexists
cplexdirexists:
	if [ -d "$(CPLEX)" ]; then echo -e "Cplex dir $(CPLEX) found. I'll compile ./main if not already compiled.\n"; else echo "\n\nCplex dir not found!!!!! \nPlease set add a line like below on \nCPLEX=/path/to/your/cplex/install/folder/\nat the beginning of this Makefile"; exit 1; fi
main:main2
main1:CuttingPlanesEngine.o general.o inout.o frontpareto1.o subprob1.o src/main.cpp 
	$(CCC) $(CCFLAGS) *.o src/main.cpp -o main $(CCLNFLAGSCPLEX) 
main2:CuttingPlanesEngine.o general.o inout.o frontpareto2.o subprob2.o src/main.cpp 
	$(CCC) $(CCFLAGS) *.o src/main.cpp -o main $(CCLNFLAGSCPLEX) 
CuttingPlanesEngine.o : ../src_shared/CuttingPlanesEngine.cpp ../src_shared/CuttingPlanesEngine.h
	$(CCC) -c $(CCFLAGSCPLEX) -Wno-ignored-attributes -DCPLEXVER=${CPLEXVERNO} ../src_shared/CuttingPlanesEngine.cpp -o CuttingPlanesEngine.o 
inout.o : src/inout.cpp src/inout.h
	$(CCC) -c $(CCFLAGS) src/inout.cpp -o inout.o 
general.o : ../src_shared/general.cpp ../src_shared/general.h
	$(CCC) -c $(CCFLAGS) ../src_shared/general.cpp -o general.o 

frontpareto1.o : src/frontpareto1.cpp src/frontpareto1.h
	rm -f frontpareto2.o
	$(CCC) -c $(CCFLAGS) src/frontpareto1.cpp -o frontpareto1.o 
frontpareto2.o : src/frontpareto2.cpp src/frontpareto2.h
	rm -f frontpareto1.o
	$(CCC) -c $(CCFLAGS) src/frontpareto2.cpp -o frontpareto2.o 
subprob1.o : src/subprob.cpp src/subprob.h src/frontpareto1.h 
	rm -f subprob2.o
	$(CCC) -DFPARETO1 -c $(CCFLAGS) -Winline src/subprob.cpp -o subprob1.o 
subprob2.o : src/subprob.cpp src/subprob.h src/frontpareto2.h 
	rm -f subprob1.o
	$(CCC) -DFPARETO2 -c $(CCFLAGS) -Winline src/subprob.cpp -o subprob2.o 


#erase below in on line version
statManager.o : src/statManager.cpp src/statManager.h
	$(CCC) -c $(CCFLAGS) src/statManager.cpp -o statManager.o 

clean: 
	rm -f main main1 main2 *.lp *.log *.err *.cout *.o

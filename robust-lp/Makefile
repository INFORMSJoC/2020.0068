#| Author: Daniel Porumbel 2018       daniel.porumbel@cnam.fr    |
#|License: Creative Commons Attribution 3.0 Unported License     |
#|         http://creativecommons.org/licenses/by/3.0/           |             
#-----------------------------------------------------------------
#Code below searches for cplex installation dir. If you know it, you can erase all
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
#If possible I always take version 12.6
ifneq ($(wildcard /opt/ibm/ILOG/CPLEX_Studio126/),) 
    CPLEX   = /opt/ibm/ILOG/CPLEX_Studio126/
    CPLEXVER="CPLEX_Studio126"
endif 
CPLEXVERNO_ALLDIGITS=$(shell echo ${CPLEXVER}|sed "s/[^0-9]//g")  #remove non digits
CPLEXVERNO=$(shell echo ${CPLEXVERNO_ALLDIGITS}|cut -c 1-3)       #take first 3 letters (eg. 126)

#----------------------------------------------------------------------------
#The rest depends on the cplex install dire defined above

CPLEXDIR      = $(CPLEX)/cplex/
CONCERTDIR    = $(CPLEX)/concert/

# ---------------------------------------------------------------------
# Compiler selection, code optimization, debug and warning options
# ---------------------------------------------------------------------
#use DNDEBUG to remove asserts and turn off some clog debug messages
CCC = g++ 
CCFLAGS = -O3 -DNDEBUG -flto -m64 -Wno-ignored-attributes -Wextra
#maybe add -Wextra and -Wpedantic
#you could add -flto for link-time optimizaton (e.g., in-lining function 
#from a .o to another .o). However, the speed-up can be about 1 or 2% or less
#maybe -Ofast that enables -ffast-math. 
#For me, the following can be good enough instead of -Ofast. However, it
#assumes things like floating point multiplication is associative, untrue.
#-funsafe-math-optimizations -ffinite-math-only


# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------
LIBFORMAT  = static_pic
#cplex related
#use SYSTEM=x86-64_sles10_4.1 for cplex versions <12.6
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

# ---------------------------------------------------------------------
# Main compilation and linking commands
# ---------------------------------------------------------------------
all:cplexdirexists main
.SILENT: cplexdirexists
cplexdirexists:
	if [ -d "$(CPLEX)" ]; then printf "\nCplex $(CPLEX) dir found! I'll compile main if not already compiled.\n"; else printf "\nCplex dir not found!!!!! \nPlease set add a line like below on \nCPLEX=/path/to/your/cplex/install/folder/\nat the beginning of this Makefile\n"; exit 1; fi
#main:CuttingPlanesEngine.o statManager.o general.o iographs.o bendersSubprob.o multiRunsLog.o src/main.cpp 
main:CuttingPlanesEngine.o subprob.o general.o inout.o src/main.cpp 
	$(CCC) $(CCFLAGS) *.o src/main.cpp -o main $(CCLNFLAGSCPLEX) 
subprob.o : src/subprob.cpp src/subprob.h
	$(CCC) -c $(CCFLAGS) src/subprob.cpp -o subprob.o 
inout.o : src/inout.cpp src/inout.h
	$(CCC) -c $(CCFLAGS) src/inout.cpp -o inout.o 
general.o : src/general.cpp src/general.h
	$(CCC) -c $(CCFLAGS) src/general.cpp -o general.o 
CuttingPlanesEngine.o : src/CuttingPlanesEngine.cpp src/CuttingPlanesEngine.h
	$(CCC) -c $(CCFLAGSCPLEX) -DCPLEXVER=${CPLEXVERNO} src/CuttingPlanesEngine.cpp -o CuttingPlanesEngine.o 

clean: cleanobj
	rm -f main *.lp *.log *.err *.cout
cleanobj:
	rm -f *.o

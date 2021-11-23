/*---------------------------------------------------------------------------------------+
|              General very useful routines, most importantly getCPUTime()               |
|        =======================================================================         |
|        - getCPUTime() is far better than unreliable clock()/CLOCKS_PER_SEC             |
|          that suffers from the "wrap around issue" described in man 3 clock:           |
|                     "On a 32-bit system where CLOCKS_PER_SEC                           |
|                      equals 1000000 this function will return                          |
|                      the same value approximately every 72 minutes"                    |
---------+-------------------------------------------------------------------------+-----+
         | See file LICENSE at the root of the git project for licence information |
         +------------------------------------------------------------------------*/

#ifndef GENERAL_H_INCLUDED
#define GENERAL_H_INCLUDED
#include <iomanip>
#include <locale>
#include <sstream>

extern "C" {
  #if defined(_WIN32)
  #include <Windows.h>
  
  #elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
  #include <unistd.h>
  #include <sys/resource.h>
  #include <sys/times.h>
  #include <time.h>
  
  #else
  #error "Unable to define getCPUTime( ) for an unknown OS."
  #endif
  
  /**
   * Returns the amount of CPU time (including parallel usage) used by the current process,
   * in *fractional* seconds, or -1.0 if an error occurred.
   */
  double getCPUTime( );
}

std::string toString(int number);

#endif

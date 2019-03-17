#ifndef _parallel_defs
#define _parallel_defs
    #if defined(_MSC_VER) && __cplusplus >= 201500
        #pragma message "parallel active , need -fopenmp flag"
        #include <execution>
        #define par_exe_msvc
    #elif defined(__GNUG__) && __cplusplus >= 201500
        #pragma message "parallel active , need -fopenmp flag"
        #define par_exe_gcc
        #include <parallel/algorithm>
        #include <parallel/numeric>
    #else
        #pragma message "parallel deactive , need -fopenmp flag and c++17 for parallel execution"
        #include <algorithm>
        #include <numeric>
    #endif
//#define show_moves
//#define debug
//#define timeit

#endif

#ifndef __lin_alg__
    #define __lin_alg__

//#if defined(_MSC_VER) && __cplusplus >= 201500
//    #pragma message "parallel active , need -fopenmp flag"
//    #include <execution>
//    #include "omp.h"
//    #define par_exe_msvc
//#elif defined(__GNUG__) && __cplusplus >= 201500
//    #pragma message "parallel active , need -fopenmp flag"
//    #define par_exe_gcc
//    #include "omp.h"
//    #include <parallel/algorithm>
//    #include <parallel/numeric>
//#else
//    #pragma message "parallel deactive , need -fopenmp flag and c++17 for parallel execution"
//    #include <algorithm>
//    #include <numeric>
//#endif

    #include <blaze/Math.h>
    #include "matrix.hpp"
    #include "vec.hpp"
    #include "solver.hpp"
    #include "gauss.hpp"
    #include "chol.hpp"
    #include "cg.hpp"
    #include "sparse_solvers.hpp"
    #include "gs.hpp"
    #include "jacobi.hpp"
    #include "blaze_solvers.hpp"
    #include "mg_helper.hpp"
    template <typename MT1  // Type of the left-hand side dense matrix operand
             ,bool SO1      // Storage order of the left-hand side dense matrix operand
             ,typename MT2  // Type of the right-hand side dense matrix operand
             ,bool SO2      // Storage order of the right-hand side dense matrix operand
             ,typename MT3  // Type of the target dense matrix
             ,bool SO3>     // Storage order of the target dense matrix
    void kronecker( const blaze::Matrix<MT1,SO1>& left
                  , const blaze::Matrix<MT2,SO2>& right
                  , blaze::Matrix<MT3,SO3>& out)
    {
      using namespace blaze;
    
      using Left  = If_t< IsMatMatMultExpr<MT1>::value, const typename MT1::ResultType, const MT1& >;
      using Right = If_t< IsExpression<MT2>::value, const typename MT2::ResultType, const MT2& >;
    
      Left  A( ~left  );
      Right B( ~right );
    
      for (std::size_t i = 0; i < rows(A); ++i) {
        for (std::size_t j = 0; j < columns(A); ++j) {
          const std::size_t row_lower = i * rows(B);
          const std::size_t col_lower = j * columns(B);
          auto block = submatrix(out, row_lower, col_lower,
                                 rows(B), columns(B));
          (~block) = A(i, j) * B;
        }
      }
    }
    
    #include <numeric>
    #include <iostream>
    #include <cassert>
    #include <limits>
    #include <cmath>
    #include <cstring>
    #include <random>
#endif

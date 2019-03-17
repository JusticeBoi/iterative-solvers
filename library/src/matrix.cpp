//#include "matrix.hpp"
#include "linalg.hpp"

namespace linalg
{
     matrix<double> generateSPDMatrix(std::size_t size)
     {
        double lower_bound = 10.0;
        double upper_bound = 100.0;
        std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
        matrix<double> tmp(size);
     
        std::default_random_engine re;
         //#if defined(par_exe_gcc) || defined(par_exe_msvc)
         //#pragma omp parallel for schedule(static) num_threads(8)  
         //#endif
         for ( size_t row = 0; row < size; ++row)
         {
             for (size_t col = 0; col < size ; ++col)
             {
                tmp(row, col ) = unif(re);
                if ( row == col ) tmp (row,col) += unif(re);
             }
         }
         auto mat(tmp);

         tmp.transposeInPlace();

         return mat * tmp;
         //matrix<double> result(size);
         //
         //#if defined(par_exe_gcc) || defined(par_exe_msvc)
         //#pragma omp parallel for schedule(static) num_threads(8)  
         //#endif
         //for ( size_t row = 0; row < size; ++row)
         //{
         //    for (size_t col = 0; col < size ; ++col)
         //    {
         //       result(row, col ) =  (0.5) * ( tmp( row, col ) + tmp( col , row ));
         //       if ( row == col ) result( row, col ) += static_cast<double>(size);
         //    }
         //}
         
     }
}


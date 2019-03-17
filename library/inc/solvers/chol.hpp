#pragma once
#include "solver.hpp"
namespace linalg
{
        template<class vec,class mat>
        class Chol : public absSolver<vec,mat>
        {
            public:
            void printDuration(){ absSolver<vec,mat>::printDuration();};
            vec solve(mat& A, vec& b, double omega = 1.0 ) override
            {
                   int i,j,k;
                   if ( A.dim() != b.size() ) throw("cholesky cant be done, sizes dont match ");
                   auto start = std::chrono::steady_clock::now(); 
                   double epsilon = std::numeric_limits<double>::epsilon();
                   int dimension = static_cast<int>(A.dim());
                   //matrix<T> tmp(A);
                   //vec<T> tmp_vec(b);
    
                       for(k = 0; k < dimension ; ++k)
                       {
    
                           if ( std::abs( A(k,k) ) <  epsilon) throw("pivot is in choleskyDecomp, we need SPD zero");  
    
                           #if defined(par_exe_gcc) || defined(par_exe_msvc)
                           #pragma omp simd 
                           #endif
                           for ( j = 0; j <= k - 1 ; ++j )
                           {
                               A( k, k ) -= A( k, j ) * A ( k, j );
                           }
    
                           A ( k, k ) = std::sqrt(std::abs(A ( k, k )));
                           //A ( k, k ) = std::sqrt(A ( k, k ));
    
                          // #if defined(par_exe_gcc) || defined(par_exe_msvc)
                          // #pragma omp parallel for schedule(static) num_threads(8) 
                          // #endif
                           #if defined(par_exe_gcc) || defined(par_exe_msvc)
                           #pragma omp parallel for schedule(static) num_threads(8) private(i)  
                           #endif
                           for(i = k + 1 ; i < dimension; ++i)
                           {
                               //#if defined(par_exe_gcc) || defined(par_exe_msvc)
                               //#pragma omp simd  
                               //#endif
                               for ( j = 0 ; j <= k - 1;  ++j)
                               {
                                   A ( i, k ) -= A ( i, j ) * A( k, j );
                               }

                               A ( i, k ) /= A ( k, k );
                           }
                       }
                       /*forward substituion */
                       for(k = 0 ; k < dimension ; ++k )
                       {
                           for ( i = 0 ; i <= k - 1 ; ++i )
                           {
                                  b[k] -= A ( k, i ) * b[i]; 
                           }
                           b[k] /= A ( k, k );
                       }
    
                       /*backward substituion */
                       for ( k = dimension - 1 ; k >= 0 ; --k)
                       {
                           for ( i = k + 1 ; i < dimension ; ++i )
                           {
                               b[k] -= A ( i, k ) * b[i];
    
                           }
                           b[k] /= A ( k, k );
                       }
    
                       
                       auto end = std::chrono::steady_clock::now(); 

                        
                       this->duration_ = std::chrono::duration_cast<std::chrono::milliseconds>(end - start) ;
                       return b;


            }
        };


}

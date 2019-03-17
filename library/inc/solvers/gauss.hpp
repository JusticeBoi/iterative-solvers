
#pragma once
#include "solver.hpp"
namespace linalg
{
    template<class vec,class mat>
    class Gauss : public absSolver<vec,mat>
    {
        public:
        void printDuration(){ absSolver<vec,mat>::printDuration();};
        vec solve(mat& A, vec& b, double omega = 1.0 ) override
        {
                try 
                {
                    auto start = std::chrono::steady_clock::now(); 
                    int k,i,j;
                    if ( A.dim() != b.size() ) throw("gauss cant be done, sizes dont match ");
                    vec  tmpVec(b);
    
                    int dimension = static_cast<int>(A.dim());
                    double epsilon = std::numeric_limits<double>::epsilon();
    
                    /* LU decomp*/
                    for (  k = 0; k < static_cast<int>(A.dim())- 1 ; ++k)
                    {
                        if ( std::abs( A(k,k) ) <  epsilon) throw("pivot is zero");  
                        #if defined(par_exe_gcc) || defined(par_exe_msvc)
                        #pragma omp parallel for schedule(static) num_threads(8) shared(A)   
                        //#pragma omp parallel shared(A) num_threads(8) 
                        //{
                        //    #pragma omp for nowait schedule(static)
                        #endif
                        for (  i = k + 1 ; i < dimension ; ++i )
                        {
                            A(i,k) /= A( k, k); 
                            #if defined(par_exe_gcc) || defined(par_exe_msvc)
                            #pragma omp  simd 
                            #endif
                            for (  j = k + 1 ; j < dimension ; ++j)
                            {
                                A(i, j ) -= ( A(i , k) * A(k, j)); 
                            }
                        }
                        //#if defined(par_exe_gcc) || defined(par_exe_msvc)
                        //}
                        //#endif
                    }
                    std::cout << A <<'\n';
    
                    /* forward substitution */
                    for (  k = 1; k < dimension ; ++k )
                    {
                        for(  i = 0 ; i < k ; ++i )
                        {
                            tmpVec[k] -= ( A( k, i ) *  tmpVec[i]);
                        }
                    }
                    /* backward substituion */
                    for ( int  k_ = dimension - 1 ; k_ >= 0 ; --k_)
                    {
                        for ( i = k_ + 1 ; i < dimension ; ++i)
                        {
                            tmpVec[k_] -= A( k_, i) * tmpVec[i]; 
                        }
                        tmpVec[k_] /= A( k_, k_ );
                    }
                    auto end = std::chrono::steady_clock::now(); 
    
                    
                    this->duration_ = std::chrono::duration_cast<std::chrono::milliseconds>(end - start) ;
    
                    return tmpVec;
                }
    
                catch(const char * exception)
                {
                    std::cerr << "Error " <<exception << '\n';
                    return vec();
                }
    
    
    
        }
    };
}

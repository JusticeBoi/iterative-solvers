
#pragma once
#include "solver.hpp"
namespace linalg
{
        template<class vec ,class mat>
        class GaussSeidel : public absSolver<vec,mat>
        {
            public:
            void printDuration(){ absSolver<vec,mat>::printDuration();};
            vec solve(mat& A, vec& b, double omega = 1.0 ) override
            {
                try
                {
    
                          #ifdef timeit
                          auto start_2 = std::chrono::steady_clock::now(); 
                          #endif

                          if (A.dim() != b.size() ) throw ("sizes dont fit");
                         vec D = A.Diag();
                        matrix L = A.Lower();
                         vec x(b.size()); 
                         //vec<T1> answer(rhs.size()); 
                         double error = 10;
                         size_t iter_no = 0;
                         /*Precompute N*b */

                         while ( error > 1e-10)
                          {
                              ++iter_no;
                              
                              //#ifdef par_exe_gcc 
                              //#pragma omp parallel for
                              //#endif
                              for ( size_t i = 0; i < b.size(); ++i )
                              {
                              #ifdef par_exe_gcc 
                                  double correction = (b[i] - __gnu_parallel::inner_product(x.begin(), x.end(), A.getRow(i).begin(),0.0)) / D[i];
                              #else
                                  double correction = (b[i] - std::inner_product(x.begin(), x.end(), A.getRow(i).begin(),0.0)) / D[i];
                              #endif
                              
                                  x[i] +=correction;
                              }
                              if( iter_no == 1) std::cout << x << '\n';
                              error =(b - A*x).normL2();
                          } 
                         #ifdef timeit
                         auto end_2 = std::chrono::steady_clock::now(); 
                         auto diff_2 = end_2 - start_2;
                         std::cout <<"duration of GS :  "<< std::chrono::duration <double, std::milli> (diff_2).count() << " ms" << std::endl; 
                         #endif
                         std::cout <<"GS finished with " << iter_no <<" iterations!" <<'\n';
                         return x;
                    } 
                catch(const char* exception)
                {

                        std::cerr << "Error : " << exception << '\n';
                        return vec();
                }

                        return vec();
            }


        };
}

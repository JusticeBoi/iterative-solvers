#pragma once
#include "solver.hpp"
namespace linalg
{
        template<class vec,class mat>
        class Jacobi : public absSolver<vec,mat>
        {
            public:
            void printDuration(){ absSolver<vec,mat>::printDuration();};
            vec solve(mat& A, vec& b, double omega = 1.0 ) override
            {
                try
                {
                        double relax = 2.0/3.0 ;
                        //double relax = 0.9 ;
                        //double relax = 1.0 ;
                        std::cout <<"JACOBI " <<'\n';
                          #ifdef timeit
                          auto start_2 = std::chrono::steady_clock::now(); 
                          #endif
                          if (A.dim() != b.size()) throw ("sizes dont fit");
                          vec D = A.Diag();
                          matrix L = A.Lower();
                         vec x (A.dim()); 
                         vec answer (A.dim()); 
                         double error = 10;
                         size_t iter_no = 0;
                         vec Nb(b.size());
                         while ( error > 1e-10)
                          {
                              ++iter_no;

                              //#ifdef par_exe_gcc 
                              //#pragma omp parallel for
                              //#endif
                              for ( size_t i = 0; i < b.size(); ++i )
                              {
                                #ifdef par_exe_gcc 
                                  Nb[i] =__gnu_parallel::inner_product(x.begin(), x.end(), A.getRow(i).begin(),0.0);
                                    #else
                                  Nb[i] =std::inner_product(x.begin(), x.end(), A.getRow(i).begin(),0.0);
                                #endif

                              }

                              #ifdef par_exe_gcc 
                              //#pragma omp parallel for
                              #pragma omp simd 
                              #endif
                              for ( size_t i = 0; i < b.size(); ++i )
                              {
                                  //T correction = (b[i] - Nb[i]) / D[i];
                                  answer[i] = x[i] + relax* (b[i] - Nb[i]) / D[i] ;
                              }

                              x = answer;
                              //x.print();
                              error = (b - A*x).normL2();
                          }  
                         std::cout <<"jacobi finished with " << iter_no <<" iterations!" <<'\n';
                         #ifdef timeit
                         auto end_2 = std::chrono::steady_clock::now();
                         auto diff_2 = end_2 - start_2;
                         std::cout <<"duration of Jacobi :  "<< std::chrono::duration <double, std::milli> (diff_2).count() << " ms" << std::endl; 
                         #endif
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

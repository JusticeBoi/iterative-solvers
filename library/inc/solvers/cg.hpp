#pragma once
#include "solver.hpp"
namespace linalg
{
   template<class vec,class mat>
   class CG : public absSolver<vec,mat>
   {
       public:
       void printDuration(){ absSolver<vec,mat>::printDuration();};
       vec solve(mat& A, vec& b, double omega = 1.0 ) override
       {
               try
               {
                  auto start = std::chrono::steady_clock::now(); 

                   if ( A.dim() != b.size() ) throw("CG cant be done, sizes dont match ");
   
                   int dimension = static_cast<int>(A.dim());
                   int i,j;
                   vec r(dimension); 
                   vec x(dimension);

                   /* r = b - A * x */
                   #if defined(par_exe_gcc) || defined(par_exe_msvc)
                   #pragma omp parallel for schedule(static) num_threads(8) 
                   #endif
                   for ( i = 0; i < dimension ; ++i  )
                   {
                       //r[i] = b[i] - std::inner_product(x.begin(), x.end(), A.getRow(i).begin(),0.0);
                       for ( j = 0; j < dimension ; ++j )
                       {
                           r[i] += A(i,j) * x[j]; 
                       }
                       r[i] = b[i] - r[i] ;
                   }
                   vec p(r);

   
                   #ifdef par_exe_gcc
                   double rsold = __gnu_parallel::inner_product(r.begin(), r.end(), r.begin(), 0.0);
                   #else
                   double rsold = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
                   #endif
                   double rsnew = rsold;

                   size_t iterCounter = 0; 
                   vec Ap(dimension,0.0);
                   while (std::sqrt(rsnew) > 1e-8 ) 
                       {
                           ++iterCounter;
                           //#define debug
                           #ifdef debug
                           std::cout <<"iter :" <<iterCounter<<" res : " << rsnew <<'\n';
                           #endif
                           #if defined(par_exe_gcc) || defined(par_exe_msvc)
                           #pragma omp parallel for schedule(static) num_threads(8) 
                           #endif
                            for ( j = 0; j < dimension ; ++j )
                            {
                           //     for ( k = 0; k < dimension ; ++k )
                           //     {
                           //         Ap[j] += A(j,k) * p[k]; 
                           //     }
                               Ap[j] = std::inner_product(p.begin(), p.end(), A.getRow(j).begin(),0.0);
                            }
                           //Ap = A * p ; 
   
                           //std::cout <<"duration of A*p for each iterations:  "<< std::chrono::duration <double, std::milli> (diff_iter).count() << " ms" << std::endl; 

                            /* alpha is the distance we want to go in the direction of p */
                           #ifdef par_exe_gcc
                           double alpha = rsold / (__gnu_parallel::inner_product(p.begin(), p.end(), Ap.begin(), 0.0));
                           #else
                           double alpha = rsold / (std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0));
                           #endif
                           
                           x += ( p * alpha );
                           r -= (Ap * alpha);
   
                           #ifdef par_exe_gcc
                           rsnew = __gnu_parallel::inner_product(r.begin(), r.end(), r.begin(), 0.0);
                           #else
                           rsnew = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
                           #endif

                           p = r + (p * ( rsnew / rsold )); 
                           //p = p * ( rsnew / rsold ); 
                           //p+=r;
                           rsold = rsnew; 
                       }

                  auto end = std::chrono::steady_clock::now(); 
                   
                  this->duration_ = std::chrono::duration_cast<std::chrono::milliseconds>(end - start) ;
   
                   //std::cout <<"Standard CG converged in " << iterCounter <<" iterations. " <<'\n';
                   return x;
               }
               catch (const char* exception)
               {
   
                   std::cerr << "Error : " << exception << '\n';
                   return vec();
               }
   
           }


       };

   template<class vec,class mat>
   class preCg : public absSolver<vec,mat>
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

                   if ( A.dim() != b.size() ) throw("CG cant be done, sizes dont match ");
   
                   int dimension = static_cast<int>(A.dim());
                   int i,j;
                   vec r(dimension); 
                   vec x(dimension); 

                   vec C = A.inverseDiagElements();
                   //vec<T1> C(dimension, 0.001 );
                   //C.print();
                   
                   /* r = b - A * x */
                   #if defined(par_exe_gcc) || defined(par_exe_msvc)
                   #pragma omp parallel for schedule(static) num_threads(8) 
                   #endif
                   for ( i = 0; i < dimension ; ++i  )
                   {
                       r[i] = b[i] - std::inner_product(x.begin(), x.end(), A.getRow(i).begin(),0.0);
                   }

                   /* zo = Cinv * r0 */

                   vec p(dimension);
                   #if defined(par_exe_gcc) || defined(par_exe_msvc)
                   __gnu_parallel::transform(r.begin(), r.end(), C.begin(), p.begin(), std::multiplies<double>());
                   #else
                   std::transform(r.begin(), r.end(), C.begin(), p.begin(), std::multiplies<double>());
                   #endif


                   #if defined(par_exe_gcc) || defined(par_exe_msvc)
                   double rsold = __gnu_parallel::inner_product(r.begin(), r.end(), p.begin(), 0.0);
                   #else
                   double rsold = std::inner_product(r.begin(), r.end(), p.begin(), 0.0);
                   #endif


                   double rsnew = rsold;

                   #ifdef timeit
                   auto end_2 = std::chrono::steady_clock::now();
                   auto diff_2 = end_2 - start_2;
   
                   std::cout <<"duration of before iterations:  "<< std::chrono::duration <double, std::milli> (diff_2).count() << " ms" << std::endl; 
   
                   start_2 = std::chrono::steady_clock::now(); 
                   #endif
                   size_t iterCounter = 0; 
                   vec Ap(dimension);
                   vec z(dimension);
                   for( int iter = 0; iter < dimension ; ++iter)
                      {
                           ++iterCounter;
                           //#define debug
                           #ifdef debug
                           std::cout <<"iter :" <<iterCounter<<" res : " << rsnew <<'\n';
                           #endif
                           #if defined(par_exe_gcc) || defined(par_exe_msvc)
                           #pragma omp parallel for schedule(static) num_threads(8) 
                           #endif
                            for ( j = 0; j < dimension ; ++j )
                            {
                               Ap[j] = std::inner_product(p.begin(), p.end(), A.getRow(j).begin(),0.0);
                            }
                           //Ap = A*p;
                           #ifdef par_exe_gcc
                           double alpha = rsold / (__gnu_parallel::inner_product(p.begin(), p.end(), Ap.begin(), 0.0));
                           #else
                           double alpha = rsold / (std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0));
                           #endif

                           
                           x += ( p * alpha );
                           r -= (Ap * alpha);

                          /* application of jacobi preconditioning */
                           /* z i the z k+1 */
                           #ifdef par_exe_gcc 
                           __gnu_parallel::transform(r.begin(), r.end(), C.begin(), z.begin(), std::multiplies<double>());
                           #else
                           std::transform(r.begin(), r.end(), C.begin(), z.begin(), std::multiplies<double>());
                           #endif
                           
                           /* r' * r */
                           #ifdef par_exe_gcc 
                           rsnew = __gnu_parallel::inner_product(r.begin(), r.end(), z.begin(), 0.0);
                           #else
                           rsnew = std::inner_product(r.begin(), r.end(), z.begin(), 0.0);
                           #endif
                           p = z + (p * ( rsnew / rsold )); 
                           rsold = rsnew; 
                       }

                   #ifdef timeit
                   end_2 = std::chrono::steady_clock::now();
                   diff_2 = end_2 - start_2;
   
                   std::cout <<"duration of all iterations:  "<< std::chrono::duration <double, std::milli> (diff_2).count() << " ms" << std::endl; 

                   start_2 = std::chrono::steady_clock::now(); 

                   #endif
   
                   std::cout <<"Preconditioned CG converged in " << iterCounter <<" iterations. " <<'\n';
                   return x;
               }
               catch (const char* exception)
               {
   
                   std::cerr << "Error : " << exception << '\n';
                   return vec();
               }


    }
   };
}

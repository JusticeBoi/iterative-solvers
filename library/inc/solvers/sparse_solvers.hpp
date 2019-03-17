#pragma once
#include "solver.hpp"
namespace linalg
{
        template<class vec,class mat>
        class sparseCG : public absSolver<vec,mat>
        {
            public:
            void printDuration(){ absSolver<vec,mat>::printDuration();};
            vec solve(mat& A, vec& b_vec, double omega = 1.0 ) override
            {
                    try
                    {
                       auto start = std::chrono::steady_clock::now(); 
                        std::vector<double> b = b_vec.toStd();
                        if ( A.getRowCount() != b.size() ) throw("CG cant be done, sizes dont match ");
    
                        int dimension = static_cast<int>(b.size());
                        std::vector<double> r(dimension); 
                        std::vector<double> x(dimension,0.0);

                        /* r = b - A * x */
                        r = A * x ; 
                        #ifdef par_exe_gcc
                        __gnu_parallel::transform(b.begin(), b.end(), r.begin(), r.begin(),std::minus<double>());
                        #else
                        std::transform(b.begin(), b.end(), r.begin(), r.begin(),std::minus<double>());
                        #endif
                        //r[i] = b[i] - r[i] ;
                        //}
                        std::vector<double> p(r);

    
                        #ifdef par_exe_gcc
                        double rsold = __gnu_parallel::inner_product(r.begin(), r.end(), r.begin(), 0.0);
                        #else
                        double rsold = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
                        #endif
                        double rsnew = rsold;

                        size_t iterCounter = 0; 
                        std::vector<double> Ap(dimension,0.0);
                        std::vector<double> tmp(dimension); 
                        while (std::sqrt(rsnew) > 1e-10 ) 
                            {
                                ++iterCounter;
                                //#define debug
                                #ifdef debug
                                std::cout <<"iter :" <<iterCounter<<" res : " << rsnew <<'\n';
                                #endif
                                //#if defined(par_exe_gcc) || defined(par_exe_msvc)
                                //#pragma omp parallel for schedule(static) num_threads(8) 
                                //#endif
                                // for ( j = 0; j < dimension ; ++j )
                                // {
                                ////     for ( k = 0; k < dimension ; ++k )
                                ////     {
                                ////         Ap[j] += A(j,k) * p[k]; 
                                ////     }
                                //    Ap[j] = std::inner_product(p.begin(), p.end(), A.getRow(j).begin(),0.0);
                                // }
                                Ap = A * p ; 
    
                                //std::cout <<"duration of A*p for each iterations:  "<< std::chrono::duration <double, std::milli> (diff_iter).count() << " ms" << std::endl; 

                                 /* alpha is the distance we want to go in the direction of p */
                                #ifdef par_exe_gcc
                                double alpha = rsold / (__gnu_parallel::inner_product(p.begin(), p.end(), Ap.begin(), 0.0));
                                #else
                                double alpha = rsold / (std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0));
                                #endif



                                /*x += ( p * alpha );*/
                                #ifdef par_exe_gcc
                                __gnu_parallel::transform(p.begin(), p.end(), tmp.begin(),
                                     std::bind(std::multiplies<double>(), std::placeholders::_1, alpha)); 
                                __gnu_parallel::transform(x.begin(), x.end(), tmp.begin(), x.begin(),std::plus<double>());
                                #else
                                std::transform(p.begin(), p.end(), tmp.begin(),
                                     std::bind(std::multiplies<double>(), std::placeholders::_1, alpha)); 
                                std::transform(x.begin(), x.end(), tmp.begin(), x.begin(),std::plus<double>());
                                #endif 

                                /*r -= (Ap * alpha);*/
                                #ifdef par_exe_gcc
                                __gnu_parallel::transform(Ap.begin(), Ap.end(), tmp.begin(),
                                     std::bind(std::multiplies<double>(), std::placeholders::_1, alpha)); 
                                __gnu_parallel::transform(r.begin(), r.end(), tmp.begin(), r.begin(),std::minus<double>());
                                #else
                                std::transform(Ap.begin(), Ap.end(), tmp.begin(),
                                     std::bind(std::multiplies<double>(), std::placeholders::_1, alpha)); 
                                std::transform(r.begin(), r.end(), tmp.begin(), r.begin(),std::minus<double>());
                                #endif 


                                #ifdef par_exe_gcc
                                rsnew = __gnu_parallel::inner_product(r.begin(), r.end(), r.begin(), 0.0);
                                #else
                                rsnew = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
                                #endif
                                

                                /*p = r + (p * ( rsnew / rsold ));*/ 
                                double r_ratio = rsnew / rsold;
                                #ifdef par_exe_gcc
                                __gnu_parallel::transform(p.begin(), p.end(), tmp.begin(),
                                     std::bind(std::multiplies<double>(), std::placeholders::_1, r_ratio)); 
                                __gnu_parallel::transform(r.begin(), r.end(), tmp.begin(), p.begin(),std::plus<double>());
                                #else
                                std::transform(p.begin(), p.end(), tmp.begin(),
                                     std::bind(std::multiplies<double>(), std::placeholders::_1, r_ratio)); 
                                std::transform(r.begin(), r.end(), tmp.begin(), p.begin(),std::plus<double>());
                                #endif
                                
                                //p = p * ( rsnew / rsold ); 
                                //p+=r;
                                rsold = rsnew; 
                            }

                       auto end = std::chrono::steady_clock::now(); 
                        
                       this->duration_ = std::chrono::duration_cast<std::chrono::milliseconds>(end - start) ;
    
                        std::cout <<"Standard CG converged in " << iterCounter <<" iterations with res "<<rsnew <<'\n';

                       double* begin = &(*x.begin());
                        return vec(begin, x.size());
                    }
                    catch (const char* exception)
                    {
    
                        std::cerr << "Error : " << exception << '\n';
                        return vec();
                    }
    
                }


            };

        template<class vec,class mat>
        class sparsepreCg : public absSolver<vec,mat>
        {
            public:
            void printDuration(){ absSolver<vec,mat>::printDuration();};
            vec solve(mat& A, vec& b, double omega = 1.0 ) override
            {
                    try
                    {
                        if ( A.getRowCount() != static_cast<int>(b.size()) ) throw("CG cant be done, sizes dont match ");
    
                        int dimension = static_cast<int>(b.size());
                        //int i,j;
                        std::vector<double> r(dimension,0.0); 
                        std::vector<double> x(dimension,0.0); 

                        auto C = A.inverseDiagElements();
                        
                        
                        /* r = b - A * x */
                        //#if defined(par_exe_gcc) || defined(par_exe_msvc)
                        //#pragma omp parallel for schedule(static) num_threads(8) 
                        //#endif
                        //for ( i = 0; i < dimension ; ++i  )
                        //{
                        //    r[i] = b[i] - std::inner_product(x.begin(), x.end(), A.getRow(i).begin(),0.0);
                        //}

                        /* zo = Cinv * r0 */

                        r = A * x ; 
                        #ifdef par_exe_gcc
                        __gnu_parallel::transform(b.begin(), b.end(), r.begin(), r.begin(),std::minus<double>());
                        #else
                        std::transform(b.begin(), b.end(), r.begin(), r.begin(),std::minus<double>());
                        #endif

                        std::vector<double> p(dimension,0.0);
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

                        size_t iterCounter = 0; 
                        std::vector<double> Ap(dimension);
                        std::vector<double> z(dimension);
                        std::vector<double> tmp(dimension,double()); 
                            //while(rsnew > 1e-15)
                        for( int iter = 0; iter < dimension ; ++iter)
                           {
                                ++iterCounter;
                                //#define debug
                                #ifdef debug
                                std::cout <<"iter :" <<iterCounter<<" res : " << rsnew <<'\n';
                                #endif
                                //#if defined(par_exe_gcc) || defined(par_exe_msvc)
                                //#pragma omp parallel for schedule(static) num_threads(8) 
                                //#endif
                                // for ( j = 0; j < dimension ; ++j )
                                // {
                                //    Ap[j] = std::inner_product(p.begin(), p.end(), A.getRow(j).begin(),0.0);
                                // }
                                //Ap = A*p;
                                Ap = A * p; 
                                #ifdef par_exe_gcc
                                double alpha = rsold / (__gnu_parallel::inner_product(p.begin(), p.end(), Ap.begin(), 0.0));
                                #else
                                double alpha = rsold / (std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0));
                                #endif

                                
                                //x += ( p * alpha );
                                //r -= (Ap * alpha);


                                /*x += ( p * alpha );*/
                                #ifdef par_exe_gcc
                                __gnu_parallel::transform(p.begin(), p.end(), tmp.begin(),
                                     std::bind(std::multiplies<double>(), std::placeholders::_1, alpha)); 
                                __gnu_parallel::transform(x.begin(), x.end(), tmp.begin(), x.begin(),std::plus<double>());
                                #else
                                std::transform(p.begin(), p.end(), tmp.begin(),
                                     std::bind(std::multiplies<double>(), std::placeholders::_1, alpha)); 
                                std::transform(x.begin(), x.end(), tmp.begin(), x.begin(),std::plus<double>());
                                #endif 
                                /*r -= (Ap * alpha);*/
                                #ifdef par_exe_gcc
                                __gnu_parallel::transform(Ap.begin(), Ap.end(), tmp.begin(),
                                     std::bind(std::multiplies<double>(), std::placeholders::_1, alpha)); 
                                __gnu_parallel::transform(r.begin(), r.end(), tmp.begin(), r.begin(),std::minus<double>());
                                #else
                                std::transform(Ap.begin(), Ap.end(), tmp.begin(),
                                     std::bind(std::multiplies<double>(), std::placeholders::_1, alpha)); 
                                std::transform(r.begin(), r.end(), tmp.begin(), r.begin(),std::minus<double>());
                                #endif 

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

                                //p = z + (p * ( rsnew / rsold )); 
                                double r_ratio = rsnew / rsold;
                                #ifdef par_exe_gcc
                                __gnu_parallel::transform(p.begin(), p.end(), tmp.begin(),
                                     std::bind(std::multiplies<double>(), std::placeholders::_1, r_ratio)); 
                                __gnu_parallel::transform(z.begin(), z.end(), tmp.begin(), p.begin(),std::plus<double>());
                                #else
                                std::transform(p.begin(), p.end(), tmp.begin(),
                                     std::bind(std::multiplies<double>(), std::placeholders::_1, r_ratio)); 
                                std::transform(z.begin(), z.end(), tmp.begin(), p.begin(),std::plus<double>());
                                #endif

                                rsold = rsnew; 
                            }

                        std::cout <<"Preconditioned CG converged in " << iterCounter <<" iterations with res "<<rsnew <<'\n';
                       double* begin = &(*x.begin());
                        return vec(begin, x.size());
                    }
                    catch (const char* exception)
                    {
    
                        std::cerr << "Error : " << exception << '\n';
                        return vec();
                    }


         }
        };
}

#pragma once
#include "solver.hpp"
using blaze::asDiagonal;
#include <blaze/math/Band.h>
namespace linalg
{
        template<class vec,class mat>
        class sparseCGBlaze : public absSolver<vec,mat>
        {
            void printDuration(){ absSolver<vec,mat>::printDuration();};
            //using spMat = blaze::CompressedMatrix<T> ; 
            //using blazeVec = blaze::DynamicVector<T> ;
            public:
            vec solve(mat& A, vec& b, double omega = 1.0 ) override
            {
                    try
                    {
                        if ( A.rows() != b.size() ) throw("CG cant be done, sizes dont match ");
                        double alpha,delta;    
                        double beta =1;    
                        omega = 2/3;
                        int dimension = static_cast<int>(b.size());
                        vec r(dimension); 
                        vec x(dimension,0.0); 

                        /* r = b - A * x */
                        r = b - A*x ;
                        vec p(r);

                        delta = (r,r);
                        size_t iterCounter = 0; 
                        vec Ap(dimension);
                        while (std::sqrt(beta) > 1e-8 ) 
                        //while (iterCounter < b.size() ) 
                            {
                                ++iterCounter;
                                //#define debug
                                #ifdef debug
                                std::cout <<"iter :" <<iterCounter<<" res : " << beta <<'\n';
                                #endif
                                Ap = A * p ; 

                                    /* alpha is the distance we want to go in the direction of p */
                                alpha = delta / (p,Ap) ;
                                x += alpha * p;
                                r -= alpha * Ap;

                                /* beta = rsnew */
                                beta = (r,r);
                                //std::cout << p <<'\n';
                                p = r + (beta / delta) * p ;
                                delta = beta;

                            }
                        
    
                        std::cout <<"Standard CG converged in " << iterCounter <<" iterations with res "<<beta <<'\n';
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
        class jacobiBlaze : public absSolver<vec,mat>
        {

            void printDuration(){ absSolver<vec,mat>::printDuration();};
            //using spMat = blaze::CompressedMatrix<T> ; 
            //using blazeVec = blaze::DynamicVector<T> ;
            public:
            vec solve(mat& A, vec& b, double omega = 1.0 ) override
            {
                std::cout <<"jacobi blaze omega " << omega << '\n';
                //DiagonalType1 D = blaze::diagonal(A);
                //mat D = blaze::diagonal(A);
                mat D(A.rows(),A.columns(),0.0);
                for(int i = 0; i < static_cast<int>(A.rows()) ; ++i) 
                {
                   D(i,i) = A(i,i); 
                }
                int dimension = static_cast<int>(b.size());
                vec x(dimension,0.0);
                blaze::IdentityMatrix<double> identity( A.rows() );
                vec N = omega * blaze::inv(D) * b;
                //mat M = blaze::inv(D) * ( D - A ) ;
                mat M = identity - (omega * blaze::inv(D) * A) ;
                size_t iterCounter = 0;
                while(blaze::l2Norm(b-A*x)> 1e-10)
                {
                    ++iterCounter;
                    x = M*x + N;
                };
                std::cout <<"jacobiBlaze converged in "<<iterCounter <<" iterations " <<'\n';
                return x;
            }
        };
        inline blaze::CompressedMatrix<double> MGrestrictOp_mixed(size_t size)
        {


            blaze::CompressedMatrix<double> IRestrict( size / 2, size);
            IRestrict.reserve( ((size/2)-2) * 5 + 6 );

            IRestrict.append( 0, 0, 1.0);             
                              
            IRestrict.append( 0, 1, 2.0);             
                              
            IRestrict.append( 0, 2, 1.0);             
            IRestrict.finalize( 0 );

            double one_sixteenth = 1.0 / 16.0;
            double thirty_sixteenth = 30.0 / 16.0;
            for ( size_t i = 1 ; i < (size / 2) - 1 ; ++i)
            {
                IRestrict.append( i, 2 * i - 1, one_sixteenth);             

                IRestrict.append( i , 2 * i , 1.0);             

                IRestrict.append( i , 2 * i + 1 , thirty_sixteenth);             

                IRestrict.append( i , 2 * i + 2 , 1.0);             

                IRestrict.append( i , 2 * i + 3 , one_sixteenth);             

                IRestrict.finalize( i );
            }

            IRestrict.append( (size/2) - 1, (size) - 3, 1.0);             
            IRestrict.append( (size/2) - 1, (size) - 2, 2.0);             
            IRestrict.append( (size/2) - 1, (size) - 1, 1.0);             

            IRestrict.finalize( (size/2) - 1 );
            return IRestrict ;
        }


        template<class vec>
        vec MGrestrict_mixed(const vec& v )
        {
            size_t size = v.size();

            return 0.25 * MGrestrictOp_mixed(size) * v ;
        }

        template<class vec>
        vec MGProlong_mixed(const vec& v )
        {
            size_t size = v.size();
            blaze::CompressedMatrix<double> IProlong = blaze::trans(MGrestrictOp_mixed(size*2 + 1));

            return 0.5 * IProlong * v;
        }

        inline blaze::CompressedMatrix<double> MGrestrictOp(size_t size)
        {

            blaze::CompressedMatrix<double> IRestrict( size / 2, size);
            IRestrict.reserve( 3 * (size / 2) );
            for ( size_t i = 0 ; i < size / 2 ; ++i)
            {
                IRestrict.append( i, 2 * i, 1.0);             

                IRestrict.append( i , 2 * i + 1 , 2.0);             

                IRestrict.append( i , 2 * i + 2 , 1.0);             

                IRestrict.finalize( i );
            }

            return IRestrict ;
        }

        template<class vec>
        vec MGrestrict(const vec& v )
        {
            size_t size = v.size();
            return 0.25 * MGrestrictOp(size) * v ;
        }

        template<class vec>
        vec MGprolong(const vec& v )
        {
            size_t size = v.size();

            blaze::CompressedMatrix<double> IProlong = blaze::trans(MGrestrictOp(size*2 + 1));
            return 0.5 * IProlong * v ;
        }
        template<class mat>
        mat MGoperator_mixed(const mat& m )
        {
            size_t size = m.rows();

            //blaze::CompressedMatrix<double> IRestrict( size / 2, size);
            blaze::CompressedMatrix<double> IRestrict = MGrestrictOp_mixed(size);

            //jblaze::CompressedMatrix<double,blaze::columnMajor> IProlong( size , size / 2 );
            blaze::CompressedMatrix<double,blaze::columnMajor> IProlong = blaze::trans(IRestrict);
            //IProlong.reserve( 3 * (size / 2) );

            //for ( size_t i = 0 ; i < size / 2 ; ++i)
            //{
            //    IProlong.append( 2 * i, i, 1.0);             

            //    IProlong.append( 2 * i + 1 , i, 2.0);             

            //    IProlong.append( 2 * i + 2 , i, 1.0);             

            //    IProlong.finalize( i );
            //}

            return 0.125 * IRestrict * m * IProlong;
        }

        template<class mat>
        mat MGoperator(const mat& m )
        {
            size_t size = m.rows();

            //blaze::CompressedMatrix<double> IRestrict( size / 2, size);
            blaze::CompressedMatrix<double> IRestrict = MGrestrictOp(size);

            //jblaze::CompressedMatrix<double,blaze::columnMajor> IProlong( size , size / 2 );
            blaze::CompressedMatrix<double,blaze::columnMajor> IProlong = blaze::trans(IRestrict);
            //IProlong.reserve( 3 * (size / 2) );

            //for ( size_t i = 0 ; i < size / 2 ; ++i)
            //{
            //    IProlong.append( 2 * i, i, 1.0);             

            //    IProlong.append( 2 * i + 1 , i, 2.0);             

            //    IProlong.append( 2 * i + 2 , i, 1.0);             

            //    IProlong.finalize( i );
            //}

            return 0.125 * IRestrict * m * IProlong;
        }

        template<class vec,class mat>
        class twoGridJsmooth : public absSolver<vec,mat>
        {
            private:
                bool matrix_initilized = false;
                std::vector<mat> my_lower_matrices;
                std::vector<blaze::DynamicMatrix<double>> my_D_inv_matrices;
                size_t my_max_depth = 0;
                double my_omega = 2.0/3.0;


            void printDuration(){ absSolver<vec,mat>::printDuration();};
            public:

            void initMatrices(const mat& initial_A,
                            std::vector<mat>& lowermatrices,
                            size_t max_depth )
            {
                lowermatrices.reserve(max_depth);
                my_D_inv_matrices.reserve(max_depth);
                mat A = initial_A;


                blaze::DynamicMatrix<double> D(initial_A.columns(),initial_A.columns());

                #pragma omp parallel for
                for ( size_t j = 0 ; j < A.columns() ; ++j) D(j,j) = 1.0 / A(j,j);
                my_D_inv_matrices.push_back(D);

                for(size_t i = 0; i < max_depth ; ++i)
                {
                    A = MGoperator(A);
                    lowermatrices.push_back(A);
                    D = blaze::DynamicMatrix<double>(A.columns(), A.columns());
                    for ( size_t j = 0 ; j < A.columns() ; ++j) D(j,j) = 1.0 / A(j,j);
                    my_D_inv_matrices.push_back(D);
                }
            }

            void initMatrices_mixed(const mat& initial_A,
                            std::vector<mat>& lowermatrices,
                            size_t max_depth )
            {
                lowermatrices.reserve(max_depth);
                my_D_inv_matrices.reserve(max_depth);
                mat A = initial_A;


                blaze::DynamicMatrix<double> D(initial_A.columns(),initial_A.columns());

                #pragma omp parallel for
                for ( size_t j = 0 ; j < A.columns() ; ++j) D(j,j) = 1.0 / A(j,j);
                my_D_inv_matrices.push_back(D);

                for(size_t i = 0; i < max_depth ; ++i)
                {
                    A = MGoperator_mixed(A);
                    lowermatrices.push_back(A);
                    D = blaze::DynamicMatrix<double>(A.columns(), A.columns());
                    for ( size_t j = 0 ; j < A.columns() ; ++j) D(j,j) = 1.0 / A(j,j);
                    my_D_inv_matrices.push_back(D);
                }
            }

            size_t findMaxDept(const mat& A)
            {
                size_t size = A.columns();
                if (! size % 2 ) std::runtime_error("size can't be divisible by  2 " );

                size_t max_depth = 1;
                while ( size != 1 )
                {
                    ++max_depth;
                    if ( size == 0 ) break;
                    (--size) /= 2;
                    std::cout <<"size : " << size <<'\n';
                }

                return max_depth - 1;
            }

            void jacobiIterations(const mat& A,
                                  const vec& b,
                                  vec& x,
                                  size_t max_iteration,
                                  size_t level_index_ascending)
            {

                
                //mat D_inv(b.size(), b.size());
                //for ( size_t i = 0 ; i < b.size() ; ++i)
                //{
                //    D_inv(i,i) = 1.0 / A(i,i);
                //}

                //auto D_inv_vector = blaze::diagonal( D_inv);
                //std::cout <<"size of computed d_inv " << D_inv_vector.size() <<" size of precomputed one " << my_D_inv_matrices[level_index_ascending].columns() <<'\n';
                auto& D_inv = my_D_inv_matrices[level_index_ascending];
                for( size_t iter = 0; iter < max_iteration ; ++iter)
                {
                    x += decldiag(D_inv) * ( b - A * x ) * my_omega;
                    //auto q = x + intermediate;
                    //x = q;
                }

            } 


            void performVcycle(const mat& A, const vec& rhs, vec& u , size_t level_index_descending )
            {
                size_t level_index_ascending = my_max_depth - level_index_descending - 1;
                if ( !matrix_initilized)
                {
                    auto start = std::chrono::steady_clock::now(); 
                    initMatrices(A,my_lower_matrices, level_index_descending ); 
                    auto end = std::chrono::steady_clock::now(); 

                    auto diff = end - start;
   
                    //std::cout <<"duration of init matrices"<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
                    //std::cout <<"size of initilized matrices : " << my_lower_matrices.size()<<'\n';
                    matrix_initilized = true;
                }
                if ( level_index_descending != 0 )
                {
                    jacobiIterations(A, rhs, u, 3,level_index_ascending);

                    vec residual = rhs - A * u;
                    mat& restricted_A = my_lower_matrices[level_index_ascending];



                    vec restricted_rhs = MGrestrict(residual);

                    vec solution_lower_level(my_lower_matrices[level_index_ascending].columns(), 0.0 );

                    //std::cout <<" going from : " << level_index_descending<< " to : " << level_index_descending - 1 <<'\n';
                    performVcycle(restricted_A, restricted_rhs, solution_lower_level, level_index_descending -1 );

                    u = u + MGprolong(solution_lower_level);

                    jacobiIterations(A, rhs, u, 3,level_index_ascending);
                }

                if ( level_index_descending == 0 )
                {
                    //std::cout <<" deepest level, size of A is : " << A.columns() << '\n';

                    auto start = std::chrono::steady_clock::now(); 
                    //u = blaze::inv(blaze::DynamicMatrix<double>(A)) * rhs ;
                    solver<sparseCGBlaze<vec,mat>,vec,mat> b(const_cast<mat&>(A),const_cast<vec&>(rhs));
                    u = b.solve();
                    auto end = std::chrono::steady_clock::now(); 

                    auto diff = end - start;
   
                    std::cout <<"duration of coarsest level solve "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
                }

            }
            void performVcycle_mixed(const mat& A, const vec& rhs, vec& u , size_t level_index_descending )
            {
                size_t level_index_ascending = my_max_depth - level_index_descending - 1;
                if ( !matrix_initilized)
                {
                    auto start = std::chrono::steady_clock::now(); 
                    initMatrices_mixed(A,my_lower_matrices, level_index_descending ); 
                    auto end = std::chrono::steady_clock::now(); 

                    auto diff = end - start;
   
                    //std::cout <<"duration of init matrices"<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
                    //std::cout <<"size of initilized matrices : " << my_lower_matrices.size()<<'\n';
                    matrix_initilized = true;
                }
                if ( level_index_descending != 0 )
                {
                    jacobiIterations(A, rhs, u, 3,level_index_ascending);

                    vec residual = rhs - A * u;
                    mat& restricted_A = my_lower_matrices[level_index_ascending];



                    vec restricted_rhs = MGrestrict_mixed(residual);

                    vec solution_lower_level(my_lower_matrices[level_index_ascending].columns(), 0.0 );

                    //std::cout <<" going from : " << level_index_descending<< " to : " << level_index_descending - 1 <<'\n';
                    performVcycle_mixed(restricted_A, restricted_rhs, solution_lower_level, level_index_descending -1 );

                    u = u + MGProlong_mixed(solution_lower_level);

                    jacobiIterations(A, rhs, u, 3,level_index_ascending);
                }

                if ( level_index_descending == 0 )
                {
                    //std::cout <<" deepest level, size of A is : " << A.columns() << '\n';

                    auto start = std::chrono::steady_clock::now(); 
                    //u = blaze::inv(blaze::DynamicMatrix<double>(A)) * rhs ;
                    solver<sparseCGBlaze<vec,mat>,vec,mat> b(const_cast<mat&>(A),const_cast<vec&>(rhs));
                    u = b.solve();
                    auto end = std::chrono::steady_clock::now(); 

                    auto diff = end - start;
   
                    std::cout <<"duration of coarsest level solve "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
                }

            }


            vec solve_new(mat&A, vec&b, size_t max_v_cycles, double tolerance )
            {
                //my_max_depth = findMaxDept(A);
                //std::cout << findMaxDept(A) <<'\n';
                my_max_depth = 1;
                std::cout <<"max depth : " << my_max_depth <<'\n';
                vec x(b.size(), 0.0);
                vec res = b;
                double initial_res_norm = blaze::norm(res);
                double residual_norm = initial_res_norm;
                double relative_residual_norm = 1.0;
                std::cout << "iter " << 0 << " |r| " << residual_norm << " |r|/|r_0| " << relative_residual_norm << '\n';

                for ( size_t iter = 0; iter < max_v_cycles ; ++iter)
                {
                    std::cout << "MG iter : " << iter <<'\n';
                    performVcycle(A, b, x, my_max_depth - 1 );

                    residual_norm = blaze::norm(b - A * x);
                    relative_residual_norm = residual_norm / initial_res_norm ;

                    std::cout << "iter " << iter + 1 << " |r| " << residual_norm << " |r|/|r_0| " << relative_residual_norm << '\n';

                    if ( residual_norm < tolerance || relative_residual_norm < tolerance )
                    {

                        std::cout << "multigrid converged after : " << iter <<" with matrix size " << b.size() << '\n';
                        std::cout <<" |r| " << residual_norm << " |r|/|r_0| " << relative_residual_norm <<'\n'; 
                        break;
                    }
                }
                return x;

            }

            vec solve_new_mixed(mat&A, vec&b, size_t max_v_cycles, double tolerance )
            {
                //my_max_depth = findMaxDept(A);
                //std::cout << findMaxDept(A) <<'\n';
                my_max_depth = 1;
                std::cout <<"max depth : " << my_max_depth <<'\n';
                vec x(b.size(), 0.0);
                vec res = b;
                double initial_res_norm = blaze::norm(res);
                double residual_norm = initial_res_norm;
                double relative_residual_norm = 1.0;
                std::cout << "iter " << 0 << " |r| " << residual_norm << " |r|/|r_0| " << relative_residual_norm << '\n';

                for ( size_t iter = 0; iter < max_v_cycles ; ++iter)
                {
                    std::cout << "MG iter : " << iter <<'\n';
                    performVcycle_mixed(A, b, x, my_max_depth - 1 );

                    residual_norm = blaze::norm(b - A * x);
                    relative_residual_norm = residual_norm / initial_res_norm ;

                    std::cout << "iter " << iter + 1 << " |r| " << residual_norm << " |r|/|r_0| " << relative_residual_norm << '\n';

                    if ( residual_norm < tolerance || relative_residual_norm < tolerance )
                    {

                        std::cout << "multigrid converged after : " << iter <<" with matrix size " << b.size() << '\n';
                        std::cout <<" |r| " << residual_norm << " |r|/|r_0| " << relative_residual_norm <<'\n'; 
                        break;
                    }
                }
                return x;

            }




            vec solve(mat& A, vec& b, double omega = 1.0 ) 
            {
                my_omega = omega;
                my_max_depth = findMaxDept(A);
                auto start = std::chrono::steady_clock::now(); 
                mat A_c = MGoperator(A);
                mat A_cc = MGoperator(A_c);
                std::cout <<"size of finest level -1 : " << A_c.columns() << " size of finest level -2 : " << A_cc.columns() << '\n';

                //blaze::DynamicVector<double,blaze::columnVector> D_inv(A.rows());

                //blaze::DynamicVector<double,blaze::columnVector> D_inv_c(A.rows());

                blaze::DiagonalMatrix< blaze::DynamicMatrix<double> > D_inv;
                D_inv.resize(b.size(), b.size());

                blaze::DiagonalMatrix< blaze::DynamicMatrix<double> > D_inv_c;
                D_inv_c.resize(b.size() / 2, b.size() / 2);

                blaze::DiagonalMatrix< blaze::DynamicMatrix<double> > D_inv_cc;
                D_inv_cc.resize(b.size() / 4, b.size() / 4);

                //mat D(A.rows(),A.columns(),0.0);
                //mat D_c(A_c.rows(),A_c.columns(),0.0);

                for(int i = 0; i < static_cast<int>(A.rows()) ; ++i) 
                {
                   D_inv(i,i) = 1 / A(i,i);
                }

                for(int i = 0; i < static_cast<int>(A_c.rows()) ; ++i) 
                {
                   D_inv_c(i,i) = 1 / A_c(i,i);
                }

                for(int i = 0; i < static_cast<int>(A_cc.rows()) ; ++i) 
                {
                   D_inv_cc(i,i) = 1 / A_cc(i,i);
                }

                int dimension = static_cast<int>(b.size());
                vec x(dimension,0.0);

                blaze::IdentityMatrix<double> identity( A.rows() );
                blaze::IdentityMatrix<double> identity_c( A_c.rows() );
                blaze::IdentityMatrix<double> identity_cc( A_cc.rows() );

                //vec tmp(b.size());
                //std::transform(b.begin(), b.end(), D.begin(), tmp.begin(),std::multiplies<double>());
                //blaze::clear(tmp);
                vec N = omega * D_inv * b; 
                mat M = identity - (omega * D_inv * A) ;


                vec r = b - A * x;

                size_t iter_counter = 0;
                std::cout<<"before loop "<<'\n';
                while(blaze::l2Norm(r) > 1e-8)
                {
                    ++iter_counter;
                    //pre smoothing
                    x = M * x + N;
                    x = M * x + N;
                    x = M * x + N;
                    //x = M * x + N;

                    //restrict the residual 
                    vec residual_c = MGrestrict(blaze::DynamicVector<double>(b - A * x) );

                    //smoothing in course grid
                    //vec N_c = omega * blaze::inv(D_c) * residual_c;
                    //mat M_c = identity_c - (omega * blaze::inv(D_c) * A_c) ;

                    //vec N_c = omega * D_inv_c * residual_c; 
                    vec N_c = omega * D_inv_c * residual_c; 

                    mat M_c = identity_c - (omega * D_inv_c * A_c) ;

                    vec error(residual_c.size(),0.0);

                    //smooth error in course grid
                    error = M_c * error + N_c;
                    error = M_c * error + N_c;
                    error = M_c * error + N_c;

                    //error = blaze::inv(blaze::DynamicMatrix<double>(A_c)) * residual_c ;
                    
                    ////-----------------------------third level -------------------------------------------------------------------------------------
                    vec residual_cc= MGrestrict(blaze::DynamicVector<double>(residual_c - A_c * error) );

                    vec N_cc = omega * D_inv_cc * residual_cc; 

                    mat M_cc = identity_cc - (omega * D_inv_cc * A_cc) ;

                    vec error_c(residual_cc.size(),0.0);

                    //smooth error in course grid
                    //error_c = M_cc * error_c + N_cc;
                    //error_c = M_cc * error_c + N_cc;
                    error_c = blaze::inv(blaze::DynamicMatrix<double>(A_cc)) * residual_cc ;
                    
                    //error correction
                    error = error + MGprolong(error_c);

                    error = M_c * error + N_c;
                    error = M_c * error + N_c;
                    error = M_c * error + N_c;
                    //error = M_c * error + N_c;
                    //-----------------end of third level -------------------------------------------------------------------------------------------------

                    //error correction
                    x = x + MGprolong(error);

                    //post smoothing
                    x = M * x + N;
                    x = M * x + N;
                    x = M * x + N;

                    r = b - A * x;
                    std::cout <<blaze::l2Norm(r)<<'\n';
                    //if(iter_counter == 3) break;
                }
                std::cout <<"converged in " << iter_counter << " steps " <<'\n';
                auto end = std::chrono::steady_clock::now(); 

                auto diff = end - start;
   
                std::cout <<"duration of all two/three grid iterations:  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
                return x;

                //two times pre smooth


            }
        };
}

#include <iostream>
#include <blaze/Math.h>
//#include <blaze/math/DynamicMatrix.h>
//#include <blaze/math/CompressedMatrix.h>
#include "matrix_generators.hpp"
#include "linalg.hpp"
#include <iomanip>

using BlazeJacobi = linalg::jacobiBlaze<blaze::DynamicVector<double>,blaze::DynamicMatrix<double>> ;
template<typename T>
using spMat = blaze::CompressedMatrix<T> ; 


using sparse = blaze::CompressedMatrix<double>;
using vec = blaze::DynamicVector<double>;

//double BCfunc(size_t i, size_t j , double delta_h )
//{
//    return std::pow((double)(delta_h * i), 2) + std::pow((double)(delta_h * j), 2 ) ;
//}

int main()
{


 // TEST   P(X,Y) = X²+Y² with domain x -> [0,1] , y -> [0,1] with

    //dimension of the matrix is n*n x n*n AFTER ELIMINATION OF BOUNDARY CONDITIONS
    size_t n = 7;
    double h = 1. / (n + 1);
    double h2 = h * h;
    double h2_inv = 1. / h2 ;
    //std::cout <<"h : " << h <<'\n';
    //tip , n+2 with boundary conditions.
    //auto BCfunc = [](size_t i, size_t j , double delta_h )-> double
    //{
    //    return std::pow((double)(delta_h * i), 3) + std::pow((double)(delta_h * j), 3 ) ;
    //};

    auto BCfunc = [](size_t i, size_t j , double delta_h )-> double
    {
        return std::pow((double)(delta_h * i), 2) + std::pow((double)(delta_h * j), 2 ) ;
    };
    auto s_f  = [](size_t i, size_t j, double delta_h )->double{return 4.0;};
    //auto s_f  = [](size_t i, size_t j,double delta_h )->double{
    //    return -delta_h*i*6 - delta_h *j*6 ;
    //};

    //blaze::CompressedMatrix<double> K(n*n, n*n);
    //blaze::DynamicVector<double> b(n*n);
    blaze::CompressedMatrix<double> K = linalg::generateD2Q5MatrixForMG(n);
    blaze::DynamicVector<double> b = linalg::generateRHSVectorForD2Q5(n,s_f,BCfunc);
    //std::cout << "b : " << b <<'\n';
    //K.reserve();

    //K.append(0, 0, 4.0);
    //K.append(0, 1, -1.0);
    //K.append(0, n, -1.0);
    //K.finalize(0);


    //for ( size_t j = 0; j < n ; ++j)
    //{
    //    for ( size_t i = 0; i < n ; ++i)
    //    {
    //        //std::cout <<"h2 * -4.0 : " << h2 * -4. <<'\n';
    //        if ( i == 0 && j == 0 )
    //        {
    //            //K.reserve(i,3);
    //            //K.append(i, 0, 4.0);
    //            //K.append(i, 1, -1.0);
    //            //K.append(i, n, -1.0);
    //            //K.finalize(i);
    //            std::cout <<" b with i : " << i << " j : " << j  << " i + n*j : "<< i + n *j  <<'\n';
    //            b[i + n*j] = h2 *  -4.0 + BCfunc(i, j+1,h) + BCfunc(i +1 ,j , h); 
    //        }
    //        else if ( i == n - 1 && j == 0)
    //        {
    //            //K.reserve(i,3);
    //            //K.append(i,i-1,-1.0);
    //            //K.append(i,i,4.0);
    //            //K.append(i,i + n,-1.0);
    //            //K.finalize(i);
    //            // i + 1 means same x axis as me, i +2 means to the right i means to the left
    //            b[i + n*j] = h2 *  -4.0 + BCfunc(i+1, j,h) + BCfunc(i +2 ,j+1 , h); 
    //            std::cout <<" b with i : " << i << " j : " << j  << " i + n*j : "<< i + n *j  <<'\n';

    //        }
    //        else if ( i == 0 && j == (n - 1))
    //        {
    //            //K.reserve(i+ n*j,3);
    //            //K.append(i+ n*j, i+ n*j-n,-1.0);
    //            //K.append(i+ n*j, i+ n*j,4.0);
    //            //K.append(i+ n*j, i+ n*j + 1,-1.0);
    //            //K.finalize(i + n*j);

    //            b[i + n*j] = h2 *  -4.0 + BCfunc(i+1, j+2,h) + BCfunc(i ,j+1 , h); 
    //            std::cout <<" b with i : " << i << " j : " << j  << " i + n*j : "<< i + n *j  <<'\n';
    //            //std::cout <<" b ["<<i + n*j <<"] = " << b[i + n*j] <<'\n';

    //        }

    //        else if (i == n - 1 && j == n - 1) 
    //        {
    //            //K.reserve(i+ n*j,3);
    //            //K.append(i + n*j , i + (n * j) - n , -1.0);
    //            //K.append(i + n*j , i + n * j - 1, -1.0);
    //            //K.append(i + n*j , i + n*j, 4.0);

    //            //K.finalize(i + n*j);

    //            b[i + n*j] = h2 *  -4.0 + BCfunc(i+1, j+2,h) + BCfunc(i + 2 ,j+1 , h); 
    //            std::cout <<" b with i : " << i << " j : " << j  << " i + n*j : "<< i + n *j  <<'\n';
    //            //std::cout <<" b ["<<i + n*j <<"] = " << b[i + n*j] <<'\n';
    //        }
    //        else if (i == 0)
    //        {
    //            //K.reserve(i+ n*j,4);
    //            //K.append(i + n*j , i + n * j - n, -1.0);
    //            //K.append(i + n*j , i + n * j, 4.0);
    //            //K.append(i + n*j , i + n * j + 1, -1.0);
    //            //K.append(i + n*j , i + n * j + n, -1.0);
    //            //K.finalize(i + n*j);

    //            b[i + n*j] = h2 *  -4.0 + BCfunc(i, j+1,h) ;
    //            std::cout <<" b with i : " << i << " j : " << j  << " i + n*j : "<< i + n *j  <<'\n';
    //            //std::cout <<" b ["<<i + n*j <<"] = " << b[i + n*j] <<'\n';
    //        }

    //        else if (j == 0)
    //        {
    //            //K.reserve(i+ n*j,4);
    //            //K.append(i + n*j , i + n * j - 1, -1.0);
    //            //K.append(i + n*j , i + n * j, 4.0);
    //            //K.append(i + n*j , i + n * j + 1, -1.0);
    //            //K.append(i + n*j , i + n * j + n, -1.0);

    //            //K.finalize(i + n*j);

    //            b[i + n*j] = h2 *  -4.0 + BCfunc(i+1, j,h) ;
    //            std::cout <<" b with i : " << i << " j : " << j  << " i + n*j : "<< i + n *j  <<'\n';
    //            //std::cout <<" b ["<<i + n*j <<"] = " << b[i + n*j] <<'\n';
    //        }

    //        else if (i == n-1)
    //        {
    //            //K.reserve(i+ n*j,4);
    //            //K.append(i + n*j , i + n * j - n, -1.0);
    //            //K.append(i + n*j , i + n * j - 1, -1.0);
    //            //K.append(i + n*j , i + n * j, 4.0);
    //            //K.append(i + n*j , i + n * j + n, -1.0);

    //            //K.finalize(i + n*j);

    //            b[i + n*j] = h2 *  -4.0 + BCfunc(i+2, j+1,h) ;
    //            std::cout <<" b with i : " << i << " j : " << j  << " i + n*j : "<< i + n *j  <<'\n';
    //            //std::cout <<" b ["<<i + n*j <<"] = " << b[i + n*j] <<'\n';
    //        }
    //        else if (j == n-1)
    //        {
    //            //K.reserve(i+ n*j,4);
    //            //K.append(i + n*j , i + n * j - n, -1.0);
    //            //K.append(i + n*j , i + n * j - 1, -1.0);
    //            //K.append(i + n*j , i + n * j, 4.0);
    //            //K.append(i + n*j , i + n * j + 1, -1.0);

    //            //K.finalize(i + n*j);
    //            b[i + n*j] = h2 *  -4.0 + BCfunc(i+1, j+2,h) ;
    //            std::cout <<" b with i : " << i << " j : " << j  << " i + n*j : "<< i + n *j  <<'\n';
    //            //std::cout <<" b ["<<i + n*j <<"] = " << b[i + n*j] <<'\n';
    //        }
    //        else 
    //        {
    //            //K.reserve(i+ n*j,5);
    //            //K.append(i + n*j , i + n * j - n, -1.0);
    //            //K.append(i + n*j , i + n * j - 1, -1.0);
    //            //K.append(i + n*j , i + n * j, 4.0);
    //            //K.append(i + n*j , i + n * j + 1, -1.0);
    //            //K.append(i + n*j , i + n * j + n, -1.0);
    //            //K.finalize(i + n*j);

    //            b[i + n*j] = h2 *  -4.0 ;
    //            std::cout <<" b with i : " << i << " j : " << j  << " i + n*j : "<< i + n *j  <<'\n';
    //            //std::cout <<" b ["<<i + n*j <<"] = " << b[i + n*j] <<'\n';
    //        }

    //    }

    //}


    //for ( size_t i = 1; i < n - 2 ; ++i)
    //{
    //    K.append(i, i - 1, -1.0);
    //    K.append(i, i, 4.0);
    //    K.append(i, i + 1, -1.0);
    //    K.append(i, i + n, -1.0);

    //   K.finalize(i);
    //}



    //std::cout <<" K : " <<'\n'<<K<<'\n';
    //std::cout <<" b : " <<'\n'<<b<<'\n';
    std::cout<< "sol : " << blaze::inv(blaze::DynamicMatrix<double>(K)) * b <<'\n';
    std::cout <<"expected sol : " ;
    for ( size_t i = 1 ; i < n+1 ; ++i )
    {
        for ( size_t j = 1; j < n+1 ; ++j )
        {
            std::cout << BCfunc(i,j,h) << '\n';

        }
    }
    linalg::solver<linalg::twoGridJsmooth<vec,sparse>,vec, sparse>solver(K,b,2./3.);

    auto start = std::chrono::steady_clock::now(); 
    auto x2 = solver.solve_new(K,b,100, 1e-8);
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout <<"duration mgrid :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
    std::cout << "x2 : " << x2 <<'\n';
    //std::cout << K << '\n';



    /////std::cout << K <<'\n';
    //using sparse = blaze::CompressedMatrix<double>;
    //using vec = blaze::DynamicVector<double>;

    //linalg::solver<linalg::twoGridJsmooth<vec,sparse>,vec, sparse>solver(K,b,2./3.);
    //linalg::solver<linalg::sparseCGBlaze<vec,sparse>,vec, sparse>solver2(K,b,2./3.);
    ////auto x = solver.solve();
    ////std::cout <<std::setprecision(15)<< x <<'\n';

    //auto start = std::chrono::steady_clock::now(); 
    //auto x2 = solver.solve_new(K,b,100, 1e-8);
    //auto end = std::chrono::steady_clock::now();
    //auto diff = end - start;
    //std::cout <<"duration mgrid :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
    //start = std::chrono::steady_clock::now(); 

    //end = std::chrono::steady_clock::now();
    //auto x3 = solver2.solve();
    //diff = end - start;
    //std::cout <<"duration cg :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
    ////std::cout <<std::setprecision(15)<< x2 <<'\n';


    ////blaze::DynamicVector<double> v2(3, 1.0);
    //v1[0] = 64.0;
    //v1[1] = 32.0;
    //v1[2] = 16.0;
    //v1[3] = 16.0;
    //v1[4] = 8.0;
    //v1[5] = 8.0;
    //v1[6] = 4.0;
    //v1[8] = 4.0;
    //auto restricted = linalg::MGrestrict<blaze::DynamicVector<double>>(v1); 
    //auto prolonged = linalg::MGprolong<blaze::DynamicVector<double>>(restricted); 
    //auto op =linalg::MGoperator<blaze::CompressedMatrix<double>>(K);
    //std::cout << restricted <<'\n';
    //std::cout << prolonged <<'\n';
    //std::cout << op <<'\n';


    return 0;
}

#include <iostream>
#include <blaze/Math.h>
#include "matrix_generators.hpp"

using BlazeJacobi = linalg::jacobiBlaze<blaze::DynamicVector<double>,blaze::DynamicMatrix<double>> ;
template<typename T>
using spMat = blaze::CompressedMatrix<T> ; 

template<typename T>
using vec = blaze::DynamicVector<T> ;

using blaze::DynamicMatrix;

int main()
{

    //blaze::DynamicMatrix<double> A( 2UL, 2UL, 0.0 );
    //A(0,0) = 4.0;
    //A(0,1) = 3.0;
    //A(1,0) = 3.0;
    //A(1,1) = 5.0;
    //blaze::DynamicVector<double> b{1.77, 2.18};


    //linalg::solver<BlazeJacobi,blaze::DynamicVector<double> ,blaze::DynamicMatrix<double>> solver(A,b);
    //auto x = solver.solve();
    //std::cout << x <<'\n';

DynamicMatrix<double> A{ { 1, 2, 3 },
                         { 4, 5, 6 },
                         { 7, 8, 9 }};

D = decldiag(A);
std::cout << D << '\n';
blaze::invert<asDiagonal>(A);
std::cout << A <<'\n';
vec<double> b { 1, 2, 3 };


std::cout << decldiag(A) * b <<'\n';
//auto A = linalg::ReadMtxFormatBlaze("mtx/1138_bus.mtx");
//vec<double> b(A.rows(),1.0);
//
//linalg::solver<linalg::sparseCGBlaze<vec<double>,spMat<double>>,vec<double> ,spMat<double>> solver(A,b);
//
//auto start = std::chrono::steady_clock::now(); 
//auto ans =solver.solve();
//auto end = std::chrono::steady_clock::now(); 
//auto diff = end - start;
//std::cout <<"duration sparseblaze  :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
//std::cout << A*ans <<'\n';
return 0;
}

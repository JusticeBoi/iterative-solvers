#include <iostream>
#include "matrix_generators.hpp"
#include "linalg.hpp"
#include <iomanip>

using BlazeJacobi = linalg::jacobiBlaze<blaze::DynamicVector<double>,blaze::DynamicMatrix<double>> ;
template<typename T>
using spMat = blaze::CompressedMatrix<T> ; 

template<typename T>
using vec = blaze::DynamicVector<T> ;

int main()
{
    // DOMAIN FROM 0 to 1  //
    //std::cout << "hello world " <<'\n';

    //std::cout << blaze::getNumThreads()<< '\n';
    size_t n = 31 ; 
    size_t n_mixed = 7 ; 
    double h = 1.0 / (n+1);
    double h_mixed = 1.0 / (n_mixed+1);
    blaze::CompressedMatrix<double> K = linalg::generateD1Q3MatrixForMG(n);
    blaze::CompressedMatrix<double> K_2 = linalg::generateD1Q5MixedMatrixForMG(n_mixed);
    std::cout.precision(15);

    auto p = [](size_t i, double h)
    {
        return std::pow((double)(h * i), 3) ;
    };
    //auto p = [](size_t i, double h)
    //{
    //    return std::pow((double)(h * i), 2) ;
    //};
    auto sec_der = [](size_t i , double h )->double {return (double)(6.*i*h);} ;

    blaze::DynamicVector<double> b = linalg::generateRHSVectorForD1Q3(n, sec_der, p);
    blaze::DynamicVector<double> b_mixed = linalg::generateRHSVectorForD1Q3(n_mixed, sec_der, p);
    using sparse = blaze::CompressedMatrix<double>;
    using vec = blaze::DynamicVector<double>;

    linalg::solver<linalg::twoGridJsmooth<vec,sparse>,vec, sparse>solver(K,b,2./3.);
    linalg::solver<linalg::sparseCGBlaze<vec,sparse>,vec, sparse>solver2(K,b,2./3.);

    blaze::DynamicVector<double> v(7);
    v[0] = 1.0;
    v[1] = 2.0;
    v[2] = 3.0;
    v[3] = 4.0;
    v[4] = 5.0;
    v[5] = 6.0;
    v[6] = 7.0;
    //std::cout <<linalg::MGrestrict_mixed(v) << '\n';
    //std::cout <<linalg::MGProlong_mixed(linalg::MGrestrict_mixed(v)) << '\n';

    //auto x = solver.solve();
    //std::cout <<std::setprecision(15)<< x <<'\n';
    auto start = std::chrono::steady_clock::now(); 
    auto x2 = solver.solve_new(K,b,100, 1e-8);
    auto x2_mixed = solver.solve_new_mixed(K_2,b_mixed,100, 1e-8);
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout <<"duration mgrid :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
    std::cout << x2 << '\n';
    std::cout << x2_mixed << '\n';
    //start = std::chrono::steady_clock::now(); 

    //end = std::chrono::steady_clock::now();
    //auto x3 = solver2.solve();
    //diff = end - start;
    //std::cout <<"duration cg :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
    //std::cout <<std::setprecision(15)<< x3 <<'\n';



    std::cout <<"expected answer : " << '\n';
    for ( size_t i = 1 ; i < n+1 ; ++i )
    {
        std::cout << p(i,h) << '\n';
    }

    std::cout <<"expected answer mixed : " << '\n';
    for ( size_t i = 1 ; i < n_mixed+1 ; ++i )
    {
        std::cout << p(i,h_mixed) << '\n';
    }
    
    








    //
    //blaze::DynamicVector<double> expected(n);
    //blaze::DynamicVector<double> expected_mixed(n_mixed);
    //std::cout <<"expected answer : " << '\n';
    //std::generate(expected.begin(), expected.end(), [&h,n = 0.0]() mutable { double r = n+h; n+=h;return r*r; }); 
    //std::generate(expected_mixed.begin(), expected_mixed.end(), [&h_mixed,n = 0.0]() mutable { double r = n+h_mixed; n+=h_mixed;return r*r; }); 
    //std::cout <<expected<<'\n';
    //std::cout <<"expected answer mixed : " << '\n';
    //std::cout <<expected_mixed<<'\n';


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

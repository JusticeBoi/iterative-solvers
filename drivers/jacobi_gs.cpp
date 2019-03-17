#include <iostream>
#include "omp.h"
#include <chrono>
#include "linalg.hpp"

int main()
{
    auto f = [](const auto x)
    {
        return x*x + 1.0 ;

    };
    size_t n = 20;
    double x_min = 1.0;
    double x_max = 12.0;
    double h = (x_max - x_min ) / ( static_cast<double>(n+1));
    std::cout << "h:" <<h<<'\n';



    linalg::matrix<double> Jacobi_3(n+2);
    linalg::matrix<double> Jacobi_5(n+2);
    linalg::vec<double> b_3(n+2,-2.0*h*h);
    linalg::vec<double> b_5(n+2,-24.0*h*h);

    b_3[0] = -2.0 * h * h + f(x_min-h); 
    b_3[n+1] = -2.0 * h * h + f(x_max+h); 

    b_5[0] = ( -24.0 * h * h ) -f(x_min-(2*h)) + ( 16.0 * f(x_min-h)); 

    b_5[1] = ( -24.0 * h * h )- f(x_min-h);

    b_5[n+1] = ( -24.0 * h * h ) - f(x_max + (2*h)) + ( 16.0 * f(x_max+h)); 
    b_5[n] = ( -24.0 * h * h ) - f(x_max+h);

    Jacobi_3(0,0) = 2.0;
    Jacobi_3(0,1) = -1.0;

    Jacobi_5(0,0) = 30.0;
    Jacobi_5(0,1) = -16.0;
    Jacobi_5(0,2) = 1.0;

    Jacobi_5(1,0) = -16.0;
    Jacobi_5(1,1) = 30.0;
    Jacobi_5(1,2) = -16.0;
    Jacobi_5(1,3) = 1.0;
    
    for (size_t i = 1; i < n+1 ; ++i )
    {
        Jacobi_3(i,i-1) = -1.0;
        Jacobi_3(i,i) = 2.0;
        Jacobi_3(i,i+1) = -1.0;
        if ( i != n && i != 1)
        {
           Jacobi_5(i,i-2) = 1.0;
           Jacobi_5(i,i-1) = -16.0;
           Jacobi_5(i,i) = 30.0;
           Jacobi_5(i,i+1) = -16.0;
           Jacobi_5(i,i+2) = 1.0;

        }
    }

    Jacobi_3(n+1,n+1) = 2.0;
    Jacobi_3(n+1,n) = -1.0;


    Jacobi_5(n+1,n+1) = 30.0;
    Jacobi_5(n+1,n) = -16.0;
    Jacobi_5(n+1,n-1) = 1.0;

    Jacobi_5(n,n+1) = -16.0;
    Jacobi_5(n,n) = 30.0;
    Jacobi_5(n,n-1) = -16.0;
    Jacobi_5(n,n-2) = 1.0;
    
    auto GS_5 = Jacobi_5;
    auto b_5_gs = b_5;
    Jacobi_3.print(); 
    Jacobi_5.print(); 
    b_3.print();
    using JacobiSolver = linalg::Jacobi<double,linalg::matrix<double>>;
    using GSSolver = linalg::GaussSeidel<double,linalg::matrix<double>>;

    linalg::solver<JacobiSolver, double,linalg::matrix<double>> solver(Jacobi_3,b_3);
    linalg::solver<JacobiSolver, double,linalg::matrix<double>> solver_5(Jacobi_5,b_5);
    //linalg::solver<GSSolver, double,linalg::matrix<double>> solver_5_gs(GS_5,b_5_gs);

    //auto ans = solver.solve();
    auto ans_5 = solver_5.solve();
    //auto ans_5_gs = solver_5_gs.solve();
    //ans.print();
    ans_5.print();
    //ans_5_gs.print();
    //b_5.print();





    return 0;
}


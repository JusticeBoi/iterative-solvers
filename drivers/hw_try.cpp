#include <iostream>
#include "omp.h"
#include <chrono>
#include "linalg.hpp"


int main()
{
    using dense = linalg::matrix<double>;
    using vector = linalg::vec<double>;
    
    using GS = linalg::GaussSeidel<vector,dense>;
    using Gauss = linalg::Gauss<vector,dense>;


    linalg::matrix<double> A = {{3, 2},
                                {2, 6}};
    linalg::vec<double> b = {2.,-8.}; 
    linalg::solver<GS,vector ,dense> solver(A,b);
   
    auto x = solver.solve();



    return 0;
}

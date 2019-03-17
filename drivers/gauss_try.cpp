#include <iostream>
#include "omp.h"
#include <chrono>
#include "linalg.hpp"


int main()
{
    using dense = linalg::matrix<double>;
    using vector = linalg::vec<double>;
    
    using Gauss = linalg::Gauss<vector,dense>;

    linalg::matrix<double> A = {{1, 1, 0},
                                {0, 2, 0},
                                {0, -1, 2}};
    linalg::vec<double> b = {2.9, 2.8, 5.3};
   
    linalg::solver<Gauss,vector ,dense> solver(A,b);
    auto x = solver.solve();



    return 0;
}

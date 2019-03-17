#include <iostream>
#include "omp.h"
#include <chrono>
#include "linalg.hpp"

int main()
{

    #pragma omp parallel
    {
        int ID = omp_get_thread_num();
        printf(" hello(%d) ", ID);
        printf(" world(%d) \n", ID);

    }

    linalg::matrix<double> A;
    A.ReadMtxFormat("mtx/nd6k.mtx");
    std::cout << A.dim() <<'\n';
    linalg::vec<double> b(A.dim(), 1.0);
    //auto orig_mat(A);


    auto start = std::chrono::steady_clock::now(); 
    std::cout <<"start"<<'\n';
    linalg::solver<linalg::CG<double,linalg::matrix<double>>, double,linalg::matrix<double>> solver(A,b);
    auto ans = solver.solve();
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout <<"duration precond cg :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
    return 0;
}

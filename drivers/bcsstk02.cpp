
#include <iostream>
#include "omp.h"
#include <chrono>
#include "linalg.hpp"
#include <iomanip>

int main()
{

    #pragma omp parallel
    {
        int ID = omp_get_thread_num();
        printf(" hello(%d) ", ID);
        printf(" world(%d) \n", ID);

    }

    linalg::matrix<double> trial1;
    trial1.ReadMtxFormat("mtx/bcsstk02.mtx");
    linalg::vec<double> x(trial1.dim(), 1.0);
    linalg::vec<double> guess(trial1.dim());
    linalg::vec<double> b(trial1.dim(), 1.0);
    //
    //std::ifstream stream("out_bcsstk02.txt");
    //double a ;
    //int counter = 0;
    //while( stream >> a)
    //{
    //    b[counter] = a;
    //    ++counter;
    //}
    //std::cout <<"b : " <<'\n';
    //b.print();

    //std::ofstream of("out_bcsstk02.txt");
    //for( size_t i = 0 ; i < trial1.dim(); ++i)
    //{
    //    for( size_t j = 0 ; j < trial1.dim(); ++j)
    //    {
    //       b[i] += (trial1(i,j) * x[j]); 
    //    }
    //    of << b[i] <<'\n';
    //}
    //of.close();
    //
    auto start = std::chrono::steady_clock::now(); 
    //linalg::vec<double> ans = linalg::preconditionedCG(trial1,b,guess);
    //linalg::vec<double> ans = linalg::choleskyDecomp(trial1,b);
    //linalg::vec<double> ans = linalg::standardCG(trial1,b,guess);
    //linalg::vec<double> ans = linalg::gaussElim(trial1,b);
    //ans.print();
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout <<"duration gauss :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 


    return 0;
}

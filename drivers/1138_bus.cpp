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

    linalg::matrix<double> trial1;
    trial1.ReadMtxFormat("mtx/1138_bus.mtx");
    auto orig_mat(trial1);

    linalg::SparseMatrix<double> A_sparse(2);
    A_sparse.ReadMtxFormat("mtx/1138_bus.mtx");
    //std::cout << linalg::checkMatrixForJacobi(trial1)<<" check " <<'\n';
    linalg::vec<double> b(A_sparse.getRowCount(), 1.0);
    //std::ofstream of("out_1138_bus.txt");
    //for( size_t i = 0 ; i < trial1.dim(); ++i)
    //{
    //    for( size_t j = 0 ; j < trial1.dim(); ++j)
    //    {
    //       b[i] += (trial1(i,j) * x[j]); 
    //    }
    //    of << b[i] <<'\n';
    //}
    //of.close();

    //linalg::solver<linalg::CG<double,linalg::matrix<double>>, double,linalg::matrix<double>> solver(trial1,b);
    auto start = std::chrono::steady_clock::now(); 






    //linalg::solver<linalg::Gauss<double,linalg::matrix<double>>, double,linalg::matrix<double>> solver(trial1,b);
    //linalg::solver<linalg::SparseGauss<double,linalg::SparseMatrix<double>>, double,linalg::SparseMatrix<double>> solver(A_sparse,b);



    //auto ans = solver.solve();
    //ans.print();
    ////linalg::vec<double> ans = linalg::choleskyDecomp(trial1,b);
    ////linalg::vec<double> ans = linalg::standardCG(trial1,b,guess);
   ////linalg::vec<double> ans = linalg::gaussElim(trial1,b);
   ////linalg::vec<double> ans = linalg::JacobiMatrix(trial1,b,guess);
   //// linalg::vec<double> ans = linalg::preconditionedCG(trial1,b,guess);
   //// for (size_t i = 0; i < ans.size() ; ++i )
   //// {
   ////     std::cout <<i+1<<' '<<ans[i] <<'\n';
   //// }

    ////for ( int i = 0; i < 5 ; ++i) std::cout <<" x : "<< ans[i] << '\n';
    //ans.print(); 
    ////ans.write();
    //auto q = orig_mat * ans; 
    //q.print();
    //auto check = orig_mat * ans ;
    //check.print();
    //linalg::vec<double> check2(ans.size());
    //for( size_t row = 0; row < ans.size() ; ++row )
    //{
    //    for (size_t col = 0; col < ans.size() ; ++col )
    //    {
    //      check2[row] += (orig_mat(row,col) * ans[col] ) ;
    //    }
    //}
    //check2.print();

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout <<"duration gauss :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
    return 0;
}

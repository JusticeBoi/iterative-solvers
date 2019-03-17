#include <iostream>
#include "omp.h"
#include <chrono>
#include "helpers.hpp"
#include "SparseMatrixMock.hpp"
#include "testslib.hpp"
#include "matrix.hpp"
#include <sstream>
#include <fstream>

int main()
{
    linalg::SparseMatrix<double> A_sparse(2);
    A_sparse.ReadMtxFormat("mtx/1138_bus.mtx");
    linalg::vec<double> b_sparse(A_sparse.getRowCount(),1.0);

    //auto start = std::chrono::steady_clock::now(); 
    //linalg::solver<linalg::CG<double,linalg::matrix<double>>, double,linalg::matrix<double>> solver(A,b);
    //auto x = solver.solve();
    //auto end = std::chrono::steady_clock::now(); 
    //auto diff = end - start;
    //std::cout <<"duration of regular matrix CG :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 

    auto start_sparse = std::chrono::steady_clock::now(); 
    //linalg::solver<linalg::sparseCG<double,linalg::SparseMatrix<double>>, double,linalg::SparseMatrix<double>> solver_sparse(A_sparse,b_sparse);
    linalg::solver<linalg::sparseCG<double,linalg::SparseMatrix<double>>, double,linalg::SparseMatrix<double>> solver_sparse(A_sparse,b_sparse);
   // linalg::solver<linalg::sparseCG<double,linalg::SparseMatrix<double>>, double,linalg::SparseMatrix<double>> solver_sparse(A_sparse,b_sparse);
    auto x_sparse = solver_sparse.solve();
    auto end_sparse = std::chrono::steady_clock::now(); 
    auto diff_sparse = end_sparse - start_sparse;
    std::cout <<"duration of sparseCG :  "<< std::chrono::duration <double, std::milli> (diff_sparse).count() << " ms" << std::endl; 
    return 0;
}

  //std::cout << sparse <<'\n';
  //for(auto a : *sparse.getRows())
  //{
  //          std::cout << a << '\n';
  //}

  //auto result = sparse.multiply(ans);

        //for (auto a : result) std::cout << a <<'\n';

	//for (int N = 0; N < 5e3; N++) {
    //    std::cout << "\rvector multiplication... #" << N + 1 << std::flush;
    //    

	//	// generate random vector and matrix
	//	int rows = rand() % 16 + 1;
	//	int cols = rand() % 16 + 1;

    //    std::vector<int> vec = linalg::generateRandomVector<int>(cols);

    //    std::vector<std::vector<int> > classicMatrix = linalg::generateRandomMatrix<int>(rows, cols);
	//	linalg::SparseMatrixMock<int> sparseMatrix = linalg::SparseMatrixMock<int>::fromvectors(classicMatrix);

	//	// calculate result manually
    //    std::vector<int> manualResult = linalg::multiplyMatrixByVector(classicMatrix, vec);
	//	// method
    //    linalg::assertEquals<std::vector<int> >(manualResult, sparseMatrix.multiply(vec), "Incorrect vector multiplication");

	//	// operator
    //    linalg::assertEquals<std::vector<int> >(manualResult, sparseMatrix * vec, "Incorrect vector multiplication (operator *)");
	//}



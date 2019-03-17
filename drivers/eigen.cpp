#include <eigen3/Eigen/Dense>
#include <chrono>
#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/unsupported/Eigen/SparseExtra>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>
#include "omp.h"
int main()
{
    //#pragma omp parallel
    //{
    //    int ID = omp_get_thread_num();
    //    printf(" hello(%d) ", ID);
    //    printf(" world(%d) \n", ID);

    //}
    Eigen::SparseMatrix<double> A;
    //Eigen::loadMarket(A, "bcsstk02.mtx");
    //Eigen::loadMarket(A, "1138_bus.mtx");
    //Eigen::loadMarket(A, "shipsec8.mtx");
    Eigen::loadMarket(A, "mtx/1138_bus.mtx");
    Eigen::SparseMatrix<double> B = A;
    Eigen::SparseMatrix<double> C = A;

    Eigen::VectorXd b = Eigen::VectorXd::Ones(A.rows(),1);
    Eigen::VectorXd c = Eigen::VectorXd::Ones(A.rows(),1);
    Eigen::VectorXd d = Eigen::VectorXd::Ones(A.rows(),1);

    Eigen::VectorXd x = Eigen::VectorXd::Zero(A.rows(),1);
    Eigen::VectorXd y = Eigen::VectorXd::Zero(A.rows(),1);
    Eigen::VectorXd z = Eigen::VectorXd(A.rows());

    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    //Eigen::SparseLU<Eigen::SparseMatrix<double>> LUsolver;
    std::cout <<"reading finished" <<'\n';
    std::cout << Eigen::nbThreads() <<'\n';
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper > cg;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver_biCG;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;

    

    //x = solver.solve(b);
    auto start = std::chrono::steady_clock::now(); 
    //solver.compute(A);
    //z = solver.solve(d);
    cg.compute(A);
    z = cg.solve(d);
    //solver_biCG.compute(A);
    std::cout <<"compute finished" <<'\n';

    //LUsolver.compute(B);
    //y = LUsolver.solve(c);
    auto end = std::chrono::steady_clock::now();
    std::cout<<cg.iterations()<<'\n';
    if(cg.info()!=Eigen::Success) {
        // decomposition failed
        std::cout <<"compute failed " <<'\n';
     return -1;
    }
    std::cout << z <<'\n';
    auto diff = end - start;
    std::cout <<"duration eigen sparseLU :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
    //std::cout << C * z <<'\n';
   //std::cout << "#iterations:     " << cg.iterations() << std::endl;
   // std::cout << "estimated error: " << cg.error()      << std::endl; 
    //std::cout << z << '\n';
    //auto ans_A = A * x;

    //auto ans_B = B * y;

    //auto ans_C = C * z;

    //std::cout << " ans_A : " << ans_A <<'\n';

    //std::cout << " ans_B : " << ans_B <<'\n';

    //std::cout << " ans_C : " << ans_C <<'\n';

     //if(solver.info()!=Eigen::Success || LUsolver.info() != Eigen::Success) {
     //   std::cout <<"solve failed " <<'\n';
     //     return -1;
     //   }

    return 0;

}



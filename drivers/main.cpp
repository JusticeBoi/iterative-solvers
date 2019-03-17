#include <iostream>
#include "omp.h"
#include <chrono>
#include "linalg.hpp"
//#include <mkl.h>

int main()
{
    #pragma omp parallel
    {
        int ID = omp_get_thread_num();
        printf(" hello(%d) ", ID);
        printf(" world(%d) \n", ID);

    }


    /*trial 1 */ 
    linalg::matrix<double> trial1;
    trial1.ReadMtxFormat("bcsstk13.mtx");
    //trial1.print();
    //linalg::matrix<double> trial_chol(trial1);
    //trial_chol.print();

    linalg::vec<double> x(trial1.dim(), 1.0);
    linalg::vec<double> b(trial1.dim());
    linalg::vec<double> c(trial1.dim());

    std::ifstream stream("out.txt");
    //double a ;
    //int counter = 0;
    //while( stream >> a)
    //{
    //    b[counter] = a;
    //    ++counter;
    //}
    //std::cout <<"b : " <<'\n';
    //b.print();
    //std::iota(x.begin(), x.end(), -1000 );
    //std::generate(x.begin(), x.end(), [n = -10.0] ()mutable { return n+=0.03; });
    //x.print();
    //std::ofstream of("out.txt");
    //for( size_t i = 0 ; i < trial1.dim(); ++i)
    //{
    //    for( size_t j = 0 ; j < trial1.dim(); ++j)
    //    {
    //       c[i] += (trial1(i,j) * x[j]); 
    //    }
    //    of << c[i] <<'\n';
    //}
    //of.close();

        
   //linalg::vec<double> trial1_vec(trial1.dim(), 0.0);
   //linalg::vec<double> trial2_vec(trial1_vec);

   //linalg::vec<double> guess(trial1.dim(), 1.0);

   //auto start = std::chrono::steady_clock::now(); 
   ////linalg::vec<double> ans = linalg::preconditionedCG(trial1,b,guess);
   //linalg::vec<double> ans = linalg::choleskyDecomp(trial1,b);
   ////linalg::vec<double> ans = linalg::standardCG(trial1,b,guess);
   //ans.print();
   //auto end = std::chrono::steady_clock::now();
   //auto diff = end - start;
   //std::cout <<"duration cg :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 

   //auto start = std::chrono::steady_clock::now(); 
   //linalg::vec<double> ans_chol = linalg::choleskyDecomp(trial_chol,b);
   //auto end = std::chrono::steady_clock::now();
   //ans_chol.print();

   //auto diff = end - start;
   //std::cout <<"duration chol :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 

   // for(int i = 0; i < 10 ; ++i)
   // {
   //     std::cout <<"from cg " << ans[i] << " from chol : " << ans_chol[i] << '\n';
   // }

    //auto start = std::chrono::steady_clock::now(); 
    //linalg::matrix<double> spd = linalg::generateSPDMatrix(10);
    //spd.print();
    //auto for_chol(spd);
    //std::cout <<"generated " <<'\n';
    ////spd.print();
    //linalg::vec<double> b(10, 200000.0);
    //linalg::vec<double> x0(10, 0.0);
    //auto end = std::chrono::steady_clock::now();
    //auto diff = end - start;
    //std::cout <<"duration generateSPDMatrix :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 



    //start = std::chrono::steady_clock::now(); 
    //linalg::vec<double> ans_cg = linalg::standardCG(spd,b,x0);
    //ans_cg.print();
    //end = std::chrono::steady_clock::now();
    //diff = end - start;
    //std::cout <<"duration cg :  "<< std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl; 
    //
    //start = std::chrono::steady_clock::now(); 
    //linalg::vec<double> ans_chol = linalg::choleskyDecomp(for_chol,b);
    //ans_chol.print();
    //end = std::chrono::steady_clock::now();

    /*   */

    /*trial 2 */ 
    //linalg::matrix<double> trial2;
    //trial2.ReadMtxFormat("cfd1.mtx");
    ////trial2.ReadMtxFormat("bcsstk13.mtx");
    //std::cout <<" start " <<'\n';
    //linalg::vec<double> trial2_vec(trial2.dim(), 0.0);
    //auto start_2 = std::chrono::steady_clock::now(); 
    //linalg::vec<double> ans_cholesky = linalg::choleskyDecomp(trial2,trial2_vec);
    //auto end_2 = std::chrono::steady_clock::now();
    //
    //auto diff_2 = end_2 - start_2;
    //std::cout <<"duration cholesky :  "<< std::chrono::duration <double, std::milli> (diff_2).count() << " ms" << std::endl; 

    //start_2 = std::chrono::steady_clock::now(); 
    //linalg::vec<double> ans_gauss = linalg::gaussElim(trial2,trial2_vec);
    //end_2 = std::chrono::steady_clock::now();
    //diff_2 = end_2 - start_2;
    //std::cout <<"duration gauss :  "<< std::chrono::duration <double, std::milli> (diff_2).count() << " ms" << std::endl; 

    //for ( size_t i = 0 ; i < trial2.dim() ; ++i )
    //{
    //    if ( std::abs(ans_gauss[i] - ans_cholesky[i] ) > 0.0001 ) std::cout<< " solutions are not equal =( " <<'\n';

    //}

    return 0;

}

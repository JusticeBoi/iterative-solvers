
#include "catch.hpp"
#include "linalg.hpp"

using dense = linalg::matrix<double>;
using sparse = linalg::SparseMatrix<double>;
using vector = linalg::vec<double>;

using Gauss = linalg::Gauss<vector,dense>;
using Chol = linalg::Chol<vector,dense>;
using CG = linalg::CG<vector,dense>;
using preCg = linalg::preCg<vector,dense>;
using Jacobi = linalg::Jacobi<vector,dense> ;
using BlazeJacobi = linalg::jacobiBlaze<blaze::DynamicVector<double>,blaze::DynamicMatrix<double>> ;
using GS = linalg::GaussSeidel<vector,dense>;
using SparsepreCg = linalg::sparsepreCg<vector, sparse>;
using SparseCg = linalg::sparseCG<vector, sparse>;

//using blazeVec = 
TEST_CASE( "easy matrix solve" )
{
    std::cout <<"Inside Solver Tests... " <<'\n';
    linalg::matrix<double> A = {{1.0, 4.0},
                                {7.0, 8.0}};
    linalg::vec<double> b = {0.9, 2.3};
   
    linalg::solver<Gauss,vector ,dense> solver(A,b);
    auto x = solver.solve();
    CHECK(std::abs( x[0] - 0.1) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( x[1] - 0.2) < std::numeric_limits<double>::epsilon() );    

}

TEST_CASE( "medium matrix solve" )
{
    linalg::matrix<double> A = {{3.5, 4.5, 5.5},
                                {2.0, 4.0, 6.0},
                                {8.0, 9.0, 9.0}};
    linalg::vec<double> b = {2.9, 2.8, 5.3};
   
    linalg::solver<Gauss,vector ,dense> solver(A,b);
    auto x = solver.solve();

    CHECK(std::abs( x[0] - 0.1) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( x[1] - 0.2) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( x[2] - 0.3) < std::numeric_limits<double>::epsilon() );    

}

TEST_CASE( "easy SPD matrix chol solve" )
{

    linalg::matrix<double> A = {{4.0, 3.0},
                                {3.0, 5.0}};
   
    linalg::vec<double> b = {1.77, 2.18};


    linalg::solver<Chol,vector ,dense> solver(A,b);
    auto x = solver.solve();

    CHECK(std::abs( x[0] - 0.21) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( x[1] - 0.31) < std::numeric_limits<double>::epsilon() );    

}

TEST_CASE( "easy SPD matrix solve chol 2" )
{

    linalg::matrix<double> A = {{4.0, 1.0,1.0},
                                {1.0, 4.0,1.0},
                                {1.0, 1.0, 4.0}};
   
    linalg::vec<double> b = {6.0, 3.0,-9.0};

    linalg::solver<Chol,vector ,dense> solver(A,b);
    auto x = solver.solve();

    CHECK(std::abs( x[0] - 2.0) < 1e-12);    
    CHECK(std::abs( x[1] - 1.0) < 1e-12);    
    CHECK(std::abs( x[2] + 3.0) < 1e-12);    

}

TEST_CASE( "easy diagonal matrix solve CG" )
{

    linalg::matrix<double> A = {{7.1, 0.0},
                                {0.0, 2.3}};
   
    linalg::vec<double> b = {8.52, 2.99};

    linalg::solver<CG,vector ,dense> solver(A,b);
    auto x = solver.solve();
    CHECK(std::abs( x[0] - 1.2) < 1e-10 );    
    CHECK(std::abs( x[1] - 1.3) < 1e-10 );    

}

TEST_CASE( "easy diagonal matrix solve PCG" )
{

    linalg::matrix<double> A = {{7.1, 0.0},
                                {0.0, 2.3}};
   
    linalg::vec<double> b = {8.52, 2.99};
    linalg::vec<double> guess = {0.0, 0.0};

    linalg::solver<preCg,vector ,dense> solver(A,b);

    auto x = solver.solve();

    CHECK(std::abs( x[0] - 1.2) < 1e-10 );    
    CHECK(std::abs( x[1] - 1.3) < 1e-10 );    

}
TEST_CASE( "hard SPD matrix solve CG" )
{

    linalg::matrix<double> A = {{30.0, -16.0, 1.0, 0.0},
                                {-16.0, 30.0, -16.0, 1.0},
                                {1.0, -16.0, 30.0, -16.0},
                                {0.0, 1.0, -16.0, 30.0 }};
                               
   
    linalg::vec<double> b = {47.0, -63.0, 63.0, -47.0};

    linalg::solver<CG,vector ,dense> solver(A,b);
    auto x = solver.solve();
    CHECK(std::abs( x[0] - 1.0) < 1e-9 );    
    CHECK(std::abs( x[1] + 1.0) < 1e-9 );    
    CHECK(std::abs( x[2] - 1.0) < 1e-9 );    
    CHECK(std::abs( x[3] + 1.0) < 1e-9 );    

}
TEST_CASE( "hard SPD matrix solve CG precond" )
{

    linalg::matrix<double> A = {{30.0, -16.0, 1.0, 0.0},
                                {-16.0, 30.0, -16.0, 1.0},
                                {1.0, -16.0, 30.0, -16.0},
                                {0.0, 1.0, -16.0, 30.0 }};
                               
   
    linalg::vec<double> b = {47.0, -63.0, 63.0, -47.0};

    linalg::solver<preCg,vector ,dense> solver(A,b);
    auto x = solver.solve();
    CHECK(std::abs( x[0] - 1.0) < 1e-9 );    
    CHECK(std::abs( x[1] + 1.0) < 1e-9 );    
    CHECK(std::abs( x[2] - 1.0) < 1e-9 );    
    CHECK(std::abs( x[3] + 1.0) < 1e-9 );    

}

TEST_CASE( "hard SPD matrix solve gauss" )
{

    linalg::matrix<double> A = {{30.0, -16.0, 1.0, 0.0},
                                {-16.0, 30.0, -16.0, 1.0},
                                {1.0, -16.0, 30.0, -16.0},
                                {0.0, 1.0, -16.0, 30.0 }};
                               
   
    linalg::vec<double> b = {47.0, -63.0, 63.0, -47.0};

    linalg::solver<Gauss,vector ,dense> solver(A,b);
    auto x = solver.solve();

    CHECK(std::abs( x[0] - 1.0) < 1e-10 );    
    CHECK(std::abs( x[1] + 1.0) < 1e-10 );    
    CHECK(std::abs( x[2] - 1.0) < 1e-10 );    
    CHECK(std::abs( x[3] + 1.0) < 1e-10 );    

}

TEST_CASE( "hard SPD matrix solve chol" )
{

    linalg::matrix<double> A = {{30.0, -16.0, 1.0, 0.0},
                                {-16.0, 30.0, -16.0, 1.0},
                                {1.0, -16.0, 30.0, -16.0},
                                {0.0, 1.0, -16.0, 30.0 }};
                               
   
    linalg::vec<double> b = {47.0, -63.0, 63.0, -47.0};

    linalg::solver<Chol,vector ,dense> solver(A,b);
    auto x = solver.solve();

    CHECK(std::abs( x[0] - 1.0) < 1e-10 );    
    CHECK(std::abs( x[1] + 1.0) < 1e-10 );    
    CHECK(std::abs( x[2] - 1.0) < 1e-10 );    
    CHECK(std::abs( x[3] + 1.0) < 1e-10 );    

}
TEST_CASE( "cg solver test " )
{
    linalg::matrix<double> A = {{4.0, 3.0},
                                {3.0, 5.0}};
   
    linalg::vec<double> b = {1.77, 2.18};


    linalg::solver<CG,vector ,dense> solver(A,b);
    auto x = solver.solve();
    CHECK(std::abs( x[0] - 0.21) < 1e-10 );    
    CHECK(std::abs( x[1] - 0.31) < 1e-10 );    

}
TEST_CASE( "cg preconditioned solver test " )
{
    linalg::matrix<double> A = {{4.0, 3.0},
                                {3.0, 5.0}};
   
    linalg::vec<double> b = {1.77, 2.18};

    linalg::solver<preCg,vector ,dense> solver(A,b);
    auto x = solver.solve();

    CHECK(std::abs( x[0] - 0.21) < 1e-10 );    
    CHECK(std::abs( x[1] - 0.31) < 1e-10 );    

}

TEST_CASE( "jacobi solver test " )
{
    linalg::matrix<double> A = {{4.0, 3.0},
                                {3.0, 5.0}};
   
    linalg::vec<double> b = {1.77, 2.18};


    linalg::solver<Jacobi,vector ,dense> solver(A,b);
    auto x = solver.solve();
    CHECK(std::abs( x[0] - 0.21) < 1e-10 );    
    CHECK(std::abs( x[1] - 0.31) < 1e-10 );    
}

TEST_CASE( "jacobi blaze solver test " )
{
    blaze::DynamicMatrix<double> A( 2UL, 2UL, 0.0 );
    A(0,0) = 4.0;
    A(0,1) = 3.0;
    A(1,0) = 3.0;
    A(1,1) = 5.0;
    blaze::DynamicVector<double> b{1.77, 2.18};


    linalg::solver<BlazeJacobi,blaze::DynamicVector<double> ,blaze::DynamicMatrix<double>> solver(A,b,0.5);
    auto x = solver.solve();
    CHECK(std::abs( x[0] - 0.21) < 1e-10 );    
    CHECK(std::abs( x[1] - 0.31) < 1e-10 );    
}

TEST_CASE( "GS solver test " )
{
    linalg::matrix<double> A = {{4.0, 3.0},
                                {3.0, 5.0}};
   
    linalg::vec<double> b = {1.77, 2.18};


    linalg::solver<GS,vector ,dense> solver(A,b);
    auto x = solver.solve();

    CHECK(std::abs( x[0] - 0.21) < 1e-10 );    
    CHECK(std::abs( x[1] - 0.31) < 1e-10 );    

}
TEST_CASE( "hard SPD matrix solve GS" )
{

    linalg::matrix<double> A = {{30.0, -16.0, 1.0, 0.0},
                                {-16.0, 30.0, -16.0, 1.0},
                                {1.0, -16.0, 30.0, -16.0},
                                {0.0, 1.0, -16.0, 30.0 }};
                               
   
    linalg::vec<double> b = {47.0, -63.0, 63.0, -47.0};

    linalg::solver<GS,vector ,dense> solver(A,b);
    auto x = solver.solve();

    CHECK(std::abs( x[0] - 1.0) < 1e-10 );    
    CHECK(std::abs( x[1] + 1.0) < 1e-10 );    
    CHECK(std::abs( x[2] - 1.0) < 1e-10 );    
    CHECK(std::abs( x[3] + 1.0) < 1e-10 );    

}
TEST_CASE( "hard SPD matrix solve precg" )
{

    linalg::matrix<double> A = {{30.0, -16.0, 1.0, 0.0},
                                {-16.0, 30.0, -16.0, 1.0},
                                {1.0, -16.0, 30.0, -16.0},
                                {0.0, 1.0, -16.0, 30.0 }};
                               
   
    linalg::vec<double> b = {47.0, -63.0, 63.0, -47.0};

    linalg::solver<preCg,vector ,dense> solver(A,b);
    auto x = solver.solve();

    CHECK(std::abs( x[0] - 1.0) < 1e-10 );    
    CHECK(std::abs( x[1] + 1.0) < 1e-10 );    
    CHECK(std::abs( x[2] - 1.0) < 1e-10 );    
    CHECK(std::abs( x[3] + 1.0) < 1e-10 );    

}
TEST_CASE( "hard SPD matrix solve cg" )
{

    linalg::matrix<double> A = {{30.0, -16.0, 1.0, 0.0},
                                {-16.0, 30.0, -16.0, 1.0},
                                {1.0, -16.0, 30.0, -16.0},
                                {0.0, 1.0, -16.0, 30.0 }};
                               
   
    linalg::vec<double> b = {47.0, -63.0, 63.0, -47.0};

    linalg::solver<CG,vector ,dense> solver(A,b);
    auto x = solver.solve();

    CHECK(std::abs( x[0] - 1.0) < 1e-10 );    
    CHECK(std::abs( x[1] + 1.0) < 1e-10 );    
    CHECK(std::abs( x[2] - 1.0) < 1e-10 );    
    CHECK(std::abs( x[3] + 1.0) < 1e-10 );    

}

TEST_CASE( "cg preconditioned sparse solver test " )
{
    linalg::SparseMatrix<double> A(2); 

    A.set(4.0,1,1);
    A.set(3.0,1,2);
    A.set(3.0,2,1);
    A.set(5.0,2,2);
    linalg::vec<double> b = {1.77, 2.18};

    linalg::solver<SparsepreCg,vector ,sparse> solver(A,b);
    auto x = solver.solve();

    CHECK(std::abs( x[0] - 0.21) < 1e-10 );    
    CHECK(std::abs( x[1] - 0.31) < 1e-10 );    

}

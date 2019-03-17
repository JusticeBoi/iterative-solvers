
#include "catch.hpp"
#include "linalg.hpp"

TEST_CASE( "DoubleMatrix default constructor" )
{
    std::cout <<"Inside Matrix Tests... " <<'\n';
    linalg::matrix<double> mat;

    CHECK( mat.dim( ) == 0 );
    
}

TEST_CASE( "DoubleMatrix  constructor" )
{
    linalg::matrix<double> mat(5);

    CHECK( mat.dim( ) == 5 );
    
}

TEST_CASE( "DoubleMatrix  constructor 2" )
{
    linalg::matrix<double> mat(7,4.0);

    CHECK( mat.dim( ) == 7 );
    for( int i = 0 ; i < 7 ; ++i )
    {
        for( int j = 0 ; j < 7 ; ++j )
        {
            CHECK( std::abs(mat(i, j) - 4.0) < std::numeric_limits<double>::epsilon()) ;
        }

    }

    
}

TEST_CASE( "DoubleMatrix  constructor 3" )
{
    linalg::matrix<double> mat(7,4.0);

    linalg::matrix<double> mat2 = mat;

    CHECK( mat.dim( ) == 7 );
    for( int i = 0 ; i < 7 ; ++i )
    {
        for( int j = 0 ; j < 7 ; ++j )
        {
            CHECK( std::abs(mat(i, j) - 4.0) < std::numeric_limits<double>::epsilon()) ;
        }

    }

    
}

TEST_CASE( "DoubleMatrix  reassign an element" )
{
    linalg::matrix<double> mat(5);
    mat(0, 1) = 2.0;

    CHECK( std::abs(mat(0, 1) - 2.0) < std::numeric_limits<double>::epsilon()) ;
    
}

TEST_CASE( "DoubleMatrix  destruction" )
{
    {
        linalg::matrix<double> mat(5);
        mat(0, 1) = 2.0;
    }

    
}
TEST_CASE( "DoubleMatrix  initializer_list construction" )
{
    linalg::matrix<double> mat = {{1.0, 2.0, 3.0},
                                {4.0, 5.0, 6.0},
                                {7.0, 8.0, 9.0}};

    CHECK( std::abs(mat(0, 0) - 1.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(0, 1) - 2.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(0, 2) - 3.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(1, 0) - 4.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(1, 1) - 5.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(1, 2) - 6.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(2, 0) - 7.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(2, 1) - 8.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(2, 2) - 9.0) < std::numeric_limits<double>::epsilon()) ;
}
TEST_CASE( "DoubleMatrix  transpose" )
{
    linalg::matrix<double> mat = {{1.0, 2.0, 3.0},
                                {4.0, 5.0, 6.0},
                                {7.0, 8.0, 9.0}};
    mat.transposeInPlace();

    CHECK( std::abs(mat(0, 0) - 1.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(0, 1) - 4.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(0, 2) - 7.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(1, 0) - 2.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(1, 1) - 5.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(1, 2) - 8.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(2, 0) - 3.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(2, 1) - 6.0) < std::numeric_limits<double>::epsilon()) ;
    CHECK( std::abs(mat(2, 2) - 9.0) < std::numeric_limits<double>::epsilon()) ;

}
TEST_CASE( "DoubleMatrix  multip" )
{
    linalg::matrix<double> mat_left = {{1.0, 2.0},
                                        {4.0, 5.0}};
                                
    linalg::matrix<double> mat_right = {{3.0, 4.0},
                                {1.0, 7.0}};
    linalg::matrix<double> result = mat_left * mat_right; 
    CHECK( mat_left.dim() == 2 );
    CHECK( mat_right.dim() == 2 );

    CHECK(std::abs( result(0,0) - 5.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( result(0,1) - 18.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( result(1,0) - 17.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( result(1,1) - 51.0) < std::numeric_limits<double>::epsilon() );    
}

TEST_CASE( "DoubleMatrix  multip with vector" )
{
    linalg::matrix<double> mat = {{1.0, 2.0, 3.0},
                                {4.0, 5.0, 6.0},
                                {7.0, 8.0, 9.0}};
                                
    linalg::vec<double> vec = {2.0, 3.0, 0.1};
    linalg::vec<double> result = mat * vec; 

    CHECK(std::abs( result[0] - 8.3) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( result[1] - 23.6) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( result[2] - 38.9) < std::numeric_limits<double>::epsilon() );    
}

TEST_CASE( "DoubleMatrix  multip with vector bigger" )
{
    linalg::matrix<double> mat_left = {{1.0, 2.0},
                                        {4.0, 5.0}};
                                
    linalg::vec<double> vec = {2.0, 3.0};
    linalg::vec<double> result = mat_left * vec; 

    CHECK(std::abs( result[0] - 8.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( result[1] - 23.0) < std::numeric_limits<double>::epsilon() );    
}

TEST_CASE( "DoubleMatrix  inverseDiagElements" )
{

    linalg::matrix<double> mat = {{1.0, 2.0, 3.0},
                                {4.0, 5.0, 6.0},
                                {7.0, 8.0, 9.0}};

    linalg::vec<double> v = mat.inverseDiagElements();

    CHECK(std::abs( v[0] - 1.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( v[1] - 1.0/5.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( v[2] - 1.0/9.0) < std::numeric_limits<double>::epsilon() );    

}

TEST_CASE( "DoubleMatrix getROw  " )
{
    linalg::matrix<double> mat = {{1.0, 2.0, 3.0},
                                {4.0, 5.0, 6.0},
                                {7.0, 8.0, 9.0}};
                                
    linalg::vec<double> row_0 =  mat.getRow(0);
    linalg::vec<double> row_1 =  mat.getRow(1);
    linalg::vec<double> row_2 =  mat.getRow(2);

    CHECK(std::abs( row_0[0] - 1.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( row_0[1] - 2.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( row_0[2] - 3.0) < std::numeric_limits<double>::epsilon() );    

    CHECK(std::abs( row_1[0] - 4.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( row_1[1] - 5.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( row_1[2] - 6.0) < std::numeric_limits<double>::epsilon() );    

    CHECK(std::abs( row_2[0] - 7.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( row_2[1] - 8.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( row_2[2] - 9.0) < std::numeric_limits<double>::epsilon() );    
}

TEST_CASE( "DoubleMatrix multiply with scalar  " )
{
    linalg::matrix<double> mat = {{1.0, 2.0, 3.0},
                                {4.0, 5.0, 6.0},
                                {7.0, 8.0, 9.0}};
    auto mat1 = mat * 3.0;                                

    CHECK(std::abs( mat1(0,0) - 3.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( mat1(0,1) - 6.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( mat1(0,2) - 9.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( mat1(1,0) - 12.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( mat1(1,1) - 15.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( mat1(1,2) - 18.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( mat1(2,0) - 21.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( mat1(2,1) - 24.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( mat1(2,2) - 27.0) < std::numeric_limits<double>::epsilon() );    
}

TEST_CASE( "DoubleMatrix upper lower" )
{
    linalg::matrix<double> mat = {{1.0, 2.0, 3.0},
                                {4.0, 5.0, 6.0},
                                {7.0, 8.0, 9.0}};
    auto lower = mat.Lower();                                
    auto upper = mat.Upper();                                

    CHECK(std::abs( lower(0,0) - 1.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( lower(1,0) - 4.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( lower(1,1) - 5.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( lower(2,0) - 7.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( lower(2,1) - 8.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( lower(2,2) - 9.0) < std::numeric_limits<double>::epsilon() );    

    CHECK(std::abs( lower(0,1) - 0.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( lower(0,2) - 0.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( lower(1,2) - 0.0) < std::numeric_limits<double>::epsilon() );    
}


TEST_CASE( "sparsematrix  multip with vector bigger" )
{
    linalg::SparseMatrix<double> mat(2);         
    mat.set(1.0,1,1);
    mat.set(2.0,1,2);
    mat.set(4.0,2,1);
    mat.set(5.0,2,2);

    std::vector<double> vec = {2.0, 3.0};
    auto result = mat* vec; 

    CHECK(std::abs( result[0] - 8.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( result[1] - 23.0) < std::numeric_limits<double>::epsilon() );    
}

TEST_CASE( "Sparse  inverseDiagElements" )
{
    linalg::SparseMatrix<double> mat(2);         
    mat.set(1.0,1,1);
    mat.set(2.0,1,2);
    mat.set(4.0,2,1);
    mat.set(5.0,2,2);

    auto v = mat.inverseDiagElements();

    CHECK(std::abs( v[0] - 1.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( v[1] - 1.0/5.0) < std::numeric_limits<double>::epsilon() );    

}

TEST_CASE( "DoubleMatrix minors" )
{
    linalg::matrix<double> mat = {{1.0, 2.0, 3.0},
                                {4.0, 5.0, 6.0},
                                {7.0, 8.0, 9.0}};
    auto minor_1 = mat.minor(0);                                
    auto minor_2 = mat.minor(1);                                
    auto minor_3 = mat.minor(2);                                

    CHECK(std::abs( minor_1(0,0) - 5.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_1(0,1) - 6.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_1(1,0) - 8.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_1(1,1) - 9.0) < std::numeric_limits<double>::epsilon() );    

    CHECK(std::abs( minor_2(0,0) - 4.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_2(0,1) - 6.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_2(1,0) - 7.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_2(1,1) - 9.0) < std::numeric_limits<double>::epsilon() );    

    CHECK(std::abs( minor_3(0,0) - 4.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_3(0,1) - 5.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_3(1,0) - 7.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_3(1,1) - 8.0) < std::numeric_limits<double>::epsilon() );    
}


TEST_CASE( "DoubleMatrix minors2" )
{
    linalg::matrix<double> A = {{30.0, -16.0, 1.0, 0.0},
                                {-16.0, 30.0, -16.0, 1.0},
                                {1.0, -16.0, 30.0, -16.0},
                                {0.0, 1.0, -16.0, 30.0 }};

    auto minor_1 = A.minor(0);                                
    auto minor_3 = A.minor(2);                                

    CHECK(std::abs( minor_1(0,0) - 30.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_1(0,1) + 16.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_1(0,2) - 1.0) < std::numeric_limits<double>::epsilon() );    

    CHECK(std::abs( minor_1(1,0) + 16.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_1(1,1) - 30.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_1(1,2) + 16) < std::numeric_limits<double>::epsilon() );    

    CHECK(std::abs( minor_1(2,0) - 1.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_1(2,1) + 16.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_1(2,2) - 30.0) < std::numeric_limits<double>::epsilon() );    

    CHECK(std::abs( minor_3(0,0) + 16.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_3(0,1) - 30.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_3(0,2) - 1.0) < std::numeric_limits<double>::epsilon() );    
                          
    CHECK(std::abs( minor_3(1,0) - 1.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_3(1,1) + 16.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_3(1,2) + 16.0) < std::numeric_limits<double>::epsilon() );    
                          
    CHECK(std::abs( minor_3(2,0) - 0.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_3(2,1) - 1.0) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( minor_3(2,2) - 30.0) < std::numeric_limits<double>::epsilon() );    

}

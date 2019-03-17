

#include "catch.hpp"
#include "linalg.hpp"

TEST_CASE( "DoubleVector_defaultInitialization" )
{
std::cout <<"Inside Vector Tests... " <<'\n';
    int n = 5;
    linalg::vec<double> vector( n );

    REQUIRE( vector.size( ) == n );
    
    for( int i = 0; i < n; ++i )
    {
        CHECK(std::abs( vector[i] - 0.0) < std::numeric_limits<double>::epsilon() );    
    }
}

TEST_CASE( "DoubleVector construction" )
{
    double * t = new double[3];
    t[0] = 4.0;
    t[1] = 5.0;
    t[2] = 6.0;

    linalg::vec<double> vector( t, 3 );

    REQUIRE( vector.size( ) == 3 );
    REQUIRE( vector.capacity( ) == 3 );
    
    for( int i = 0; i < 3; ++i )
    {
        CHECK(std::abs( vector[i] - double(i+4)) < std::numeric_limits<double>::epsilon() );    
    }
}



TEST_CASE( "DoubleVector_Initialization & dot product" )
{
    int n = 3;
    linalg::vec<double> vector1( n, 2.0 );
    linalg::vec<double> vector2( n ,3.0 );
    
    REQUIRE( vector1.size( ) == n );
    REQUIRE( vector2.size( ) == n );


    for( int i = 0; i < n; ++i )
    {
        CHECK(std::abs( vector1[i] - 2.0) < std::numeric_limits<double>::epsilon() );    
        CHECK(std::abs( vector2[i] - 3.0) < std::numeric_limits<double>::epsilon() );    
    }

    CHECK(std::abs( (vector1 * vector2) - 18.0 )   < std::numeric_limits<double>::epsilon() );    

}
TEST_CASE( "DoubleVector default initalization" )
{
    
    linalg::vec<double> vector1 = { 1.0, 2.0, 3.0 };
    linalg::vec<double> vector3{ 5.0, 7.0, 2.0 };
}
TEST_CASE( "DoubleVector move ctor" )
{
    
    linalg::vec<double> vector1 = { 1.0, 2.0, 3.0 };
    linalg::vec<double> vector3(std::move(vector1));

    CHECK(vector1.size() == 0);    
    CHECK(vector1.capacity() == 0);    
    CHECK(vector1.begin() == nullptr);    

    CHECK(vector3.size() == 3);    
    CHECK(vector3.capacity() == 3);    

    CHECK(std::abs( *(vector3.end() - 1)  - 3.0 ) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( *(vector3.end() - 2)  - 2.0 ) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( vector3[0]  - 1.0 ) < std::numeric_limits<double>::epsilon() );    

}

TEST_CASE( "DoubleVector copy construction" )
{
    
    linalg::vec<double> vector1 = { 1.0, 2.0, 3.0 };
    linalg::vec<double> vector3(vector1);

    CHECK(vector3.size() == 3);    
    CHECK(vector3.capacity() == 3);    
    CHECK(std::abs( vector3[2]  - 3.0 ) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( vector3[1]  - 2.0 ) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( vector3[0]  - 1.0 ) < std::numeric_limits<double>::epsilon() );    

    CHECK(vector1.size() == 3);    
    CHECK(vector1.capacity() == 3);    
    CHECK(std::abs( vector1[2]  - 3.0 ) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( vector1[1]  - 2.0 ) < std::numeric_limits<double>::epsilon() );    
    CHECK(std::abs( vector1[0]  - 1.0 ) < std::numeric_limits<double>::epsilon() );    
}
TEST_CASE( "DoubleVector move assignment" )
{
    
    linalg::vec<double> vector1 = { 1.0, 2.0, 3.0 };
    linalg::vec<double> vector3;

    vector3 = std::move(vector1); 

    CHECK(vector1.size() == 0);    
    CHECK(vector1.capacity() == 0);    
    CHECK(vector1.begin() == nullptr);    

    CHECK(vector3.size() == 3);    
    CHECK(vector3.capacity() == 3);    
    CHECK(std::abs( *(vector3.end() - 1)  - 3.0 ) < std::numeric_limits<double>::epsilon() );    

}

TEST_CASE( "DoubleVector copy assignment" )
{

    linalg::vec<double> vector1 = { 4.0, 5.0, 6.0 };

    linalg::vec<double> vector2 = vector1;
    

    CHECK(vector2.size() == 3);    
    CHECK(vector2.capacity() == 3);    

    for( int i = 0; i < 3; ++i )
    {
        CHECK(std::abs( vector2[i] - double(i+4)) < std::numeric_limits<double>::epsilon() );    
    }


}

TEST_CASE( "DoubleVector copy assignment2" )
{

    linalg::vec<double> vector1 = { 4.0, 5.0, 6.0 };

    linalg::vec<double> vector2 = {1.0,2.0,3.0};
    

    vector2 = vector1;
    CHECK(vector2.size() == 3);    
    CHECK(vector2.capacity() == 3);    

    for( int i = 0; i < 3; ++i )
    {
        CHECK(std::abs( vector2[i] - double(i+4)) < std::numeric_limits<double>::epsilon() );    
    }


}
TEST_CASE( "DoubleVector reset" )
{

    linalg::vec<double> vector1 = { 4.0, 5.0, 6.0 };
    
    vector1.reset();
    
    CHECK(vector1.size() == 0);    
    CHECK(vector1.capacity() == 0);    
    CHECK(vector1.begin() == nullptr);    

}

TEST_CASE( "DoubleVector initializer_list assignment" )
{
    
    linalg::vec<double> vector1; 

    CHECK(vector1.size() == 0);    
    CHECK(vector1.capacity() == 0);    
    CHECK(vector1.begin() == nullptr);    
    
    vector1 = { 1.0, 2.0, 3.0};

    CHECK(vector1.size() == 3);    
    CHECK(vector1.capacity() == 3);    
    CHECK(vector1.begin() != nullptr);    

    for( int i = 0; i < 3; ++i )
    {
        CHECK(std::abs( vector1[i] - double(i+1)) < std::numeric_limits<double>::epsilon() );    
    }


}

TEST_CASE( "DoubleVector push back" )
{

    linalg::vec<double> vector1 = { 1.0, 2.0, 3.0 };
    vector1.push_back(4.0);

    CHECK(vector1.size() == 4);    
    CHECK(vector1.capacity() >= 4);    

    for( int i = 0; i < 4; ++i )
    {
        CHECK(std::abs( vector1[i] - double(i+1)) < std::numeric_limits<double>::epsilon() );    
    }
}


TEST_CASE( "DoubleVector reserve" )
{

    linalg::vec<double> vector1 = { 1.0, 2.0, 3.0, 4.0 };

    vector1.reserve(10);

    CHECK(vector1.capacity() == 10);    
    CHECK(vector1.size() == 4);    
    for( int i = 0; i < 4; ++i )
    {
        CHECK(std::abs( vector1[i] - double(i+1)) < std::numeric_limits<double>::epsilon() );    
    }
}


TEST_CASE( "DoubleVector copy ctor and ranged for loop" )
{
    int n = 3;
    linalg::vec<double> vector1( n, 3.0 );
    auto vector2(vector1);
    
    for (auto b : vector2)
    {
        CHECK(std::abs( b - 3.0 )   < std::numeric_limits<double>::epsilon() );    
    }

}

TEST_CASE( "DoubleVector sorting" )
{
    linalg::vec<double> vector1 {5.0, 1.0, 3.0, 0.0, 2.0 };

    linalg::vec<double> vector2 {0.0, 1.0, 2.0, 3.0, 5.0 };
    
    std::sort(vector1.begin(), vector1.end());

    for( int i = 0; i < 5; ++i )
    {
        CHECK(std::abs( vector1[i] -vector2[i] ) < std::numeric_limits<double>::epsilon() );    
    }

}
TEST_CASE( "DoubleVector norm" )
{
    linalg::vec<double> vector1 {5.0, 1.0, 3.0, 0.0, 2.0 };

    CHECK(vector1.size() == 5);    
    CHECK(vector1.capacity() >= 5);    

    CHECK(std::abs( vector1.normL2()- std::sqrt(39.0)) < std::numeric_limits<double>::epsilon() );    

}

TEST_CASE( "DoubleVector operator * " )
{
    linalg::vec<double> vector1 {5.0, 1.0, 3.0};
    linalg::vec<double> vector2 {2.0, 2.0, 2.0};


    CHECK(std::abs( vector1 * vector2 -  18.0 ) < std::numeric_limits<double>::epsilon() );    

}
TEST_CASE( "DoubleVector comparison" )
{
    linalg::vec<double> vector1 {5.0, 1.0, 3.0, 0.0, 2.0 };
    linalg::vec<double> vector2 {5.0, 1.0, 3.0, 0.0, 2.0 };
    linalg::vec<double> vector3(vector2); 
    linalg::vec<double> vector4;
    vector4 = vector3;

    CHECK(vector1 == vector2);    
    CHECK(vector1 == vector3);    
    CHECK(vector2 == vector3);    
    CHECK(vector1 == vector4);    
    CHECK(vector4 == vector3);    
    CHECK(vector4 == vector2);    

}

TEST_CASE( "DoubleVector += operator" )
{
    linalg::vec<double> vector1 {5.0, 1.0, 3.0, 0.0, 2.0 };
    linalg::vec<double> vector2 {2.0, 3.0, 1.0, 6.0, -2.0 };
    linalg::vec<double> vector3 {7.0, 4.0, 4.0, 6.0, 0.0 };

    CHECK((vector1 += vector2) == vector3);    

}

TEST_CASE( "DoubleVector -= operator" )
{
    linalg::vec<double> vector1 {5.0, 1.0, 3.0, 0.0, 2.0 };
    linalg::vec<double> vector2 {5.0, 3.0, 1.0, 6.0, -2.0 };
    linalg::vec<double> vector3 {0.0, -2.0, 2.0, -6.0, 4.0 };

    CHECK((vector1 -= vector2) == vector3);    

}

TEST_CASE( "DoubleVector + operator" )
{
    linalg::vec<double> vector1 {5.0, 1.0, 3.0, 0.0, 2.0 };
    linalg::vec<double> vector2 {2.0, -3.0, 1.0, 6.0, -2.0 };
    linalg::vec<double> vector3 {7.0, -2.0, 4.0, 6.0, 0.0 };

    CHECK((vector1 + vector2) == vector3);    
}

TEST_CASE( "DoubleVector - operator" )
{
    linalg::vec<double> vector1 {5.0, 1.0, 3.0, 0.0, 2.0 };
    linalg::vec<double> vector2 {2.0, -3.0, 1.0, 6.0, -2.0 };
    linalg::vec<double> vector3 {3.0, 4.0, 2.0, -6.0, 4.0 };

    CHECK((vector1 - vector2) == vector3);    
}

TEST_CASE( "DoubleVector * with operator" )
{
    linalg::vec<double> vector1 {-5.0, 1.3, 3.0, 0.0, 2.0 };
    double a = 3.0;
    linalg::vec<double> vector2 {-15.0, 3.9, 9.0, 0.0, 6.0 };
    CHECK( vector1 * a == vector2 );   

}
TEST_CASE( "DoubleVector multip with std transform" )
{
    linalg::vec<double> vector1 {5.0, 1.2, 3.0, 0.0, 2.0 };
    linalg::vec<double> vector2 {-15.0, 3.0, 9.0, 0.0, 6.0 };
    linalg::vec<double> vector3 {-75.0, 3.6, 27.0, 0.0, 12.0 };
    std::transform(vector1.begin(), vector1.end(), vector2.begin() ,vector1.begin(), std::multiplies<double>());

    CHECK(std::abs( vector1[0] - vector3[0]   ) < 1e-15);    
    CHECK(std::abs( vector1[1] - vector3[1]   ) < 1e-15);    
    CHECK(std::abs( vector1[2] - vector3[2]   ) < 1e-15);    

}

TEST_CASE( "DoubleVector tostd" )
{
    linalg::vec<double> vector1 {5.0, 1.0, 3.0, 0.0, 2.0 };
    std::vector<double> vector2 = vector1.toStd();

    CHECK(std::abs( vector1[0] - vector2[0]   ) < 1e-18);    
    CHECK(std::abs( vector1[1] - vector2[1]   ) < 1e-18);    
    CHECK(std::abs( vector1[2] - vector2[2]   ) < 1e-18);    
    CHECK(std::abs( vector1[3] - vector2[3]   ) < 1e-18);    
    CHECK(std::abs( vector1[4] - vector2[4]   ) < 1e-18);    

}

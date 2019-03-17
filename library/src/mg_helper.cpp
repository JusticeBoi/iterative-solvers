#include "../inc/mg_helper.hpp"
namespace linalg
{
    blaze::CompressedMatrix<double> generateD1Q3MatrixForMG( size_t n )
    {

        blaze::CompressedMatrix<double> K(n, n);

        K.reserve((n - 2) * 3 + 4 );

        K.append(0, 0, 2.0);
        K.append(0, 1, -1.0);

        K.finalize(0);

        for ( size_t i = 1; i < K.rows() -1 ; ++i)
        {
           K.append(i, i-1, -1.0); 
           K.append(i, i, 2.0); 
           K.append(i, i+1, -1.0); 

           K.finalize(i);
        }

        K.append(n-1, n-2, -1.0);
        K.append(n-1, n-1, 2.0);

        K.finalize(n-1);
        return K;
    }
    blaze::CompressedMatrix<double> generateD1Q5MixedMatrixForMG( size_t n )
    {
        assert( n > 4 );
        double one_twelwth = 1. / 12. ;
        blaze::CompressedMatrix<double> K(n, n);
        K.reserve((n - 2) * 5 );

        K.append(0, 0, 2.0);
        K.append(0, 1, -1.0);
        K.finalize(0);

        K.append(1, 0, -1.0);
        K.append(1, 1, 2.0);
        K.append(1, 2, -1.0);
        K.finalize(1);

        for ( size_t i = 2; i < K.rows() - 2 ; ++i)
        {
           K.append(i, i - 2, one_twelwth); 
           K.append(i, i - 1, -16. * one_twelwth); 
           K.append(i, i , 30. * one_twelwth); 
           K.append(i, i + 1, -16. * one_twelwth); 
           K.append(i, i + 2, one_twelwth); 

           K.finalize(i);

        }
        K.append(n - 2,  n - 3, -1.0);
        K.append(n - 2,  n - 2 , 2.0);
        K.append(n - 2,  n - 1, -1.0);
        K.finalize(n - 2);

        K.append(n - 1, n-2, -1.0);
        K.append(n - 1, n-1, 2.0);
        K.finalize(n - 1);

        return K;
    }
    blaze::CompressedMatrix<double> generateD2Q5MatrixForMG( size_t n )
    {
        blaze::CompressedMatrix<double> K(n*n, n*n);
        for ( size_t j = 0; j < n ; ++j)
        {
            for ( size_t i = 0; i < n ; ++i)
            {
                if ( i == 0 && j == 0 )
                {
                    K.reserve(i,3);
                    K.append(i, 0, 4.0);
                    K.append(i, 1, -1.0);
                    K.append(i, n, -1.0);
                    K.finalize(i);
                }
                else if ( i == n - 1 && j == 0)
                {
                    K.reserve(i,3);
                    K.append(i,i-1,-1.0);
                    K.append(i,i,4.0);
                    K.append(i,i + n,-1.0);
                    K.finalize(i);
                    
                    // i + 1 means same x axis as me, i +2 means to the right i means to the left

                }
                else if ( i == 0 && j == (n - 1))
                {
                    K.reserve(i+ n*j,3);
                    K.append(i+ n*j, i+ n*j-n,-1.0);
                    K.append(i+ n*j, i+ n*j,4.0);
                    K.append(i+ n*j, i+ n*j + 1,-1.0);
                    K.finalize(i + n*j);


                }

                else if (i == n - 1 && j == n - 1) 
                {
                    K.reserve(i+ n*j,3);
                    K.append(i + n*j , i + (n * j) - n , -1.0);
                    K.append(i + n*j , i + n * j - 1, -1.0);
                    K.append(i + n*j , i + n*j, 4.0);

                    K.finalize(i + n*j);

                }
                else if (i == 0)
                {
                    K.reserve(i+ n*j,4);
                    K.append(i + n*j , i + n * j - n, -1.0);
                    K.append(i + n*j , i + n * j, 4.0);
                    K.append(i + n*j , i + n * j + 1, -1.0);
                    K.append(i + n*j , i + n * j + n, -1.0);
                    K.finalize(i + n*j);

                }

                else if (j == 0)
                {
                    K.reserve(i+ n*j,4);
                    K.append(i + n*j , i + n * j - 1, -1.0);
                    K.append(i + n*j , i + n * j, 4.0);
                    K.append(i + n*j , i + n * j + 1, -1.0);
                    K.append(i + n*j , i + n * j + n, -1.0);

                    K.finalize(i + n*j);

                }

                else if (i == n-1)
                {
                    K.reserve(i+ n*j,4);
                    K.append(i + n*j , i + n * j - n, -1.0);
                    K.append(i + n*j , i + n * j - 1, -1.0);
                    K.append(i + n*j , i + n * j, 4.0);
                    K.append(i + n*j , i + n * j + n, -1.0);

                    K.finalize(i + n*j);

                }
                else if (j == n-1)
                {
                    K.reserve(i+ n*j,4);
                    K.append(i + n*j , i + n * j - n, -1.0);
                    K.append(i + n*j , i + n * j - 1, -1.0);
                    K.append(i + n*j , i + n * j, 4.0);
                    K.append(i + n*j , i + n * j + 1, -1.0);

                    K.finalize(i + n*j);
                }
                else 
                {
                    K.reserve(i+ n*j,5);
                    K.append(i + n*j , i + n * j - n, -1.0);
                    K.append(i + n*j , i + n * j - 1, -1.0);
                    K.append(i + n*j , i + n * j, 4.0);
                    K.append(i + n*j , i + n * j + 1, -1.0);
                    K.append(i + n*j , i + n * j + n, -1.0);
                    K.finalize(i + n*j);

                }

            }

        }
        return K;
    }

}

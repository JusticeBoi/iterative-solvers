
template <typename second_der_func, typename orig_func>
blaze::DynamicVector<double> generateRHSVectorForD2Q5(size_t n, second_der_func&& s_f, orig_func&& f )
{
    double h = 1. / (n + 1);
    double h2 = h * h;
    blaze::DynamicVector<double> b(n*n);
    for ( size_t j = 0; j < n ; ++j)
    {
        for ( size_t i = 0; i < n ; ++i)
        {
            if ( i == 0 && j == 0 )
            {
                b[i + n*j] = h2 *  -s_f(i+1,j+1,h) + f(i, j+1,h) + f(i +1 ,j , h); 
            }
            else if ( i == n - 1 && j == 0)
            {
                b[i + n*j] = h2 *  -s_f(i+1,j+1,h) + f(i+1, j,h) + f(i +2 ,j+1 , h); 

            }
            else if ( i == 0 && j == (n - 1))
            {
                b[i + n*j] = h2 *  -s_f(i+1,j+1,h) + f(i+1, j+2,h) + f(i ,j+1 , h); 
            }

            else if (i == n - 1 && j == n - 1) 
            {
                b[i + n*j] = h2 *  -s_f(i+1,j+1,h) + f(i+1, j+2,h) + f(i + 2 ,j+1 , h); 
            }
            else if (i == 0)
            {
                b[i + n*j] = h2 *  -s_f(i+1,j+1,h) + f(i, j+1,h) ;
            }

            else if (j == 0)
            {
                b[i + n*j] = h2 *  -s_f(i+1,j+1,h) + f(i+1, j,h) ;
            }

            else if (i == n-1)
            {
                b[i + n*j] = h2 *  -s_f(i+1,j+1,h) + f(i+2, j+1,h) ;
            }
            else if (j == n-1)
            {
                b[i + n*j] = h2 *  -s_f(i+1,j+1,h) + f(i+1, j+2,h) ;
            }
            else 
            {
                b[i + n*j] = h2 *  -s_f(i+1,j+1,h) ;
            }

        }

    }
    return b;
}
    template <typename second_der_func, typename orig_func>
    blaze::DynamicVector<double> generateRHSVectorForD1Q3(size_t n, second_der_func&& s_f, orig_func&& f )
{
    blaze::DynamicVector<double> b(n);
    double h = 1. / (n + 1);
    double h2 = h * h;

    b[0] = -s_f(1,h) * h2 + f(0,h);

    for(size_t i = 1; i <n-1 ; ++i)
    {
        b[i] = -s_f(i+1,h) * h2;
    }

    b[n-1] = -s_f(n,h) * h2 + f(n+1,h);
    return b;
}


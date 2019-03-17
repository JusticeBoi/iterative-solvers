    //#ifndef __solver__hpp
    //#define __solver__hpp
    #pragma once
    #include <chrono>
    #include "linalg.hpp"
    #include <iomanip>
    #include "parallelDefs.hpp"
    #include <blaze/math/Band.h>
    #include <blaze/math/DiagonalMatrix.h>
    #include <blaze/math/IdentityMatrix.h>
  
     #define timeit

    namespace linalg
    {
        template<typename T>
        struct Methodtraits;
        
        template<typename T>
        using methodSupportsSparse = typename Methodtraits<T>::isSparse;

        template <class Method, class vec, class mat >
        class solver : public Method
        {
            public:
                solver(mat& A, vec& b, double omega = 1.0): A_(A), b_(b),m_(Method{}), omega_(omega){}; 
                vec solve()
                {
                    return m_.solve(A_,b_,omega_);
                }
            private:
                mat& A_;
                vec& b_;
                Method m_; 
                double omega_;
        };

       
        template <class vec, class mat>
        class absSolver
        {
            public:
                virtual vec solve(mat& A, vec& b, double omega = 1.0) = 0;
                virtual ~absSolver(){};
                virtual void printDuration() = 0;

            protected:
                std::chrono::duration <double, std::milli> duration_; 

        };
        template <class vec, class mat>
        void absSolver<vec,mat>::printDuration()
          {
              std::cout <<"Duration : " << duration_.count() <<" ms. " <<'\n';
          }


    } 

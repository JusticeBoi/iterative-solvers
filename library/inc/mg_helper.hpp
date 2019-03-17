#pragma once
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicVector.h> 
#include <iostream>
namespace linalg
{
    blaze::CompressedMatrix<double> generateD1Q3MatrixForMG( size_t n );
    blaze::CompressedMatrix<double> generateD1Q5MixedMatrixForMG( size_t n );
    blaze::CompressedMatrix<double> generateD2Q5MatrixForMG( size_t n );

    template <typename second_der_func, typename orig_func>
    blaze::DynamicVector<double> generateRHSVectorForD2Q5(size_t n, second_der_func&& s_f, orig_func&& f );
    template <typename second_der_func, typename orig_func>
    blaze::DynamicVector<double> generateRHSVectorForD1Q3(size_t n, second_der_func&& s_f, orig_func&& f );

    #include "mg_helper_impl.hpp"
}


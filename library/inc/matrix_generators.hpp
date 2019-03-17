    //#ifndef __INC__MATRIX_GENERATORS_HPP
    //#define __INC__MATRIX_GENERATORS_HPP
#pragma once
#include "linalg.hpp"
#include "parallelDefs.hpp"

namespace linalg
{

    template<typename T>
    matrix<T> eye(size_t dim )
    {
        matrix<T> result(dim);    

        #ifdef par_exe_gcc 
        #pragma omp parallel for
        #endif
        for(size_t i = 0 ; i < dim ; ++i)
        {
            result(i,i) = 1.0;
        }

            return result;
    }
    template<typename T>
    using spMat = blaze::CompressedMatrix<T> ; 

    template<typename T>
    using vector = blaze::DynamicVector<T> ;

    spMat<double> ReadMtxFormatBlaze( std::string path )
    {
            std::ifstream read(path);
            if(read.is_open())
            {
                std::string line;
                bool banner = true;
                bool is_sym = false;

                while (banner )
                {
                    std::getline(read,line);
                    if ( line.substr(0,2) == "%%" && line.find("symmetric") != std::string::npos )
                    {
                            is_sym = true;
                    }

                   if ( std::isdigit(line[0]) ) break;  
                }
                int tmp_row = 0; 
                int total_no_of_entries = 0; 
                int tmp_col = 0; 
                double tmp_val = 0.0; 

                std::istringstream iss(line);
                iss >> tmp_row >> tmp_col >> total_no_of_entries;

                iss.clear();


                if ( tmp_row != tmp_col || !tmp_row || !tmp_col ) 
                {
                    std::cout <<" read something wrong, returning! " <<'\n';
                    return spMat<double>();
                }


                spMat<double> A(tmp_row,tmp_col);

                while( read >> tmp_row)
                {
                    read >> tmp_col;
                    read >> tmp_val;

                    A(tmp_row-1, tmp_col-1 ) = tmp_val;

                    if ( is_sym && tmp_col != tmp_row ) 
                    {
                        A(tmp_col-1, tmp_row-1 ) = tmp_val;
                    }
                }
                return A;


            }
            else 
            {
                std::cout <<"not open " <<'\n';
            }
            read.close();


    }
    
}


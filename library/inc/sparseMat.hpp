//#ifndef __sparse_matrix__hpp
//#define __sparse_matrix__hpp
//
//#include <random>
//#include <iostream>
//#include <cmath>
//#include <string>
//#include <fstream>
//#include <sstream>
//#include <cctype>
//#include <cstring>
//#include <functional>
//#include <limits>
//#include "parallelDefs.hpp"
//namespace linalg
//{
//
//    template<typename T>
//    class sparse
//    {
//        private:
//            std::vector<T>* vals_;
//            std::vector<T>* cols_;
//            std::vector<T>* rowptr_;
//            size_t numberOfRows_;
//            size_t numberOfCols_;
//        public:
//            /*constuctor */
//            sparse();
//            sparse(std::initializer_list<std::initializer_list<T>> init);
//            sparse(size_t dim);
//            sparse(size_t dim, T initialValue );
//            sparse(const sparse& _rhs);
//            sparse(sparse && ) noexcept;
//            ~sparse();
//            void destruct();
//            void construct(size_t row, size_t col);
//            /*assignments*/
//            sparse& operator=( sparse rhs );
//            sparse& operator=(std::initializer_list<T> reinit );
//
//            T& operator()(size_t row, size_t col);
//            const T& operator()(size_t row, size_t col) const;
//    };
//
//    
//
//    template<typename T>
//    void sparse<T>::construct(size_t rows, size_t cols)
//    {
//    	if (rows < 1 || cols < 1) {
//    		throw ("Matrix dimensions cannot be zero or negative.");
//    	}
//
//        this->numberOfRows_ = rows;
//        this->numberOfCols_ = cols;
//
//        this->vals_ = nullptr;
//        this->cols_ = nullptr;
//        this->rowptr_ = new std::vector<T>(rows+1,1);
//    }
//    template<typename T>
//    void sparse<T>::destruct()
//    {
//    	if (this->vals != nullptr) {
//    		delete this->vals_;
//    		delete this->cols_;
//    	}
//    
//    	delete this->rowptr_;
//
//    }
//    template<typename T>
//    sparse<T>::~sparse()
//    {
//        this->destruct();
//    }
//
//
//
//}
//#endif

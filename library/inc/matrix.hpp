#ifndef __matrix__hpp
#define __matrix__hpp

#include <random>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <cctype>
#include <cstring>
#include <functional>
#include <limits>
#include "parallelDefs.hpp"
#include <cassert>
#include <iomanip>

namespace linalg
{

    template<typename T >
    class vec;

    
    template<typename T>
    class matrix
    {
        private:
            T** _2darray;
            //vec<T> _mat;
            size_t _dim;

            void assignValue( T val);

            template<typename T1>
            friend void swap( matrix<T1>& lhs, matrix<T1>& rhs);
        public:

            /*constuctor */
            matrix();
            matrix(std::initializer_list<std::initializer_list<T>> init);
            matrix(size_t dim);
            matrix(size_t dim, T initialValue );
            matrix(const matrix& _rhs);
            matrix(matrix && ) noexcept;
            void allocateMemForMatrix();
            /*assignments*/
            matrix& operator=( matrix rhs );
            matrix& operator=(std::initializer_list<T> reinit );

            T& operator()(size_t row, size_t col);
            const T& operator()(size_t row, size_t col) const;
            matrix& transposeInPlace ( ); 
            matrix operator*( const matrix& mat ); 
            vec<T> operator*( vec<T>& rhs);
            matrix<T> operator*( T t );

            matrix<T> lowerInvert();

            vec<T> inverseDiagElements();
            vec<T> Diag();
            matrix<T> Lower();
            matrix<T> Upper();
            matrix<T> minor(size_t);
            double det();
            bool checkMatrixForJacobi( matrix<T>& mat);

            vec<T> getRow(size_t row ) const;
            vec<T> getRow(size_t row ) ;
            /*destruction*/
            void reset();
            ~matrix();
            void ReadMtxFormat( std::string path);
            void print() const ;
            /*getters*/
            size_t dim() const;
    };

        /* default constructor */
        template <typename T>
        matrix<T>::matrix():_2darray(nullptr),_dim(0)/*_mat(vec<T>()*/
        {
        }
        
        template <typename T>
        void matrix<T>::assignValue(T val)
        {
            //std::fill(_mat.begin(), _mat.end(), val );
            for (size_t row = 0; row < _dim ; ++row)
            {
                for (size_t col = 0; col < _dim ; ++col )
                {
                    _2darray[row][col] = val;
                }
            }
        }

        /* constructor */
        template<typename T>
        matrix<T>::matrix(size_t dim):_dim(dim)
        {
            allocateMemForMatrix();
            assignValue(0.0);
        }
        
        /* constructor */
        template<typename T>
        matrix<T>::matrix(size_t dim, T val):_dim(dim)
        {
            allocateMemForMatrix();
            assignValue(val);
        }

        /* reaching/assigning an element */
        template<typename T>
        inline T& matrix<T>::operator()(size_t row, size_t col)
        {
            return _2darray[row][col];

        }
        /*copy constructor */
        template <typename T>
        matrix<T>::matrix(const matrix& rhs):_dim(rhs._dim)
        {
            allocateMemForMatrix();
            //std::copy(rhs._mat.begin() , rhs._mat.end(), _mat.begin());
            for (size_t i = 0; i < _dim ; ++i)
            {
                std::memcpy(_2darray[i], rhs._2darray[i], _dim * sizeof(T));
            }
        }

        /* reaching/assigning an element */
        template<typename T>
        inline const T& matrix<T>::operator()(size_t row, size_t col) const
        {
            //return _mat[col*_dim + row];
            return _2darray[row][col];

        }

        /* initializer list constuctor */
        template<typename T>
        matrix<T>::matrix(std::initializer_list<std::initializer_list<T>> listlist )
        {
            if( (listlist.begin())->size() != listlist.size()) throw("onyl square matrix is allowed");
            _dim = listlist.size();
            allocateMemForMatrix();
            for ( size_t i = 0 ; i < listlist.size() ; ++i)
            {
                for ( size_t j = 0; j < listlist.size(); ++j )
                {
                   _2darray[i][j] = ((listlist.begin() + i)->begin())[j]; 
                }
            }

            
        }
        template <typename T>
        matrix<T>& matrix<T>::operator=( matrix rhs )
        {
            using std::swap;

            swap(_2darray,rhs._2darray);
            swap(_dim,rhs._dim);

            return *this;
        }

        /*allocate memory for the matrix, with known dim already */ 
        template <typename T>
        inline void matrix<T>::allocateMemForMatrix()
        {
            //_mat = vec<T>(_dim * _dim );
            if ( _dim != 0)
            {
                _2darray = new T*[_dim];
                for (size_t row = 0; row < _dim ; ++row) 
                {
                    _2darray[row] = new T[_dim]; 
                }

            }
            else 
            {
                throw("memory for the matrix could not be allocated, you used this function before assiging _dim ");
            }

        }

        /* reset */
        template <typename T>
        void matrix<T>::reset()
        {
            //_mat.reset();
            for ( size_t i = 0; i < _dim ; ++i)
            {
                delete[] _2darray[i];
            }

            delete[] _2darray;

            _dim = 0;
            
        }

        /*destructor */
        template <typename T>
        matrix<T>::~matrix()
        {
            reset();
        }


        template <typename T>
        inline size_t matrix<T>::dim() const
        { return _dim;}
        
        template <typename T>
        void matrix<T>::print() const
        {
            if(_dim != 0)
            {
                for (size_t row = 0; row < _dim ; ++row)
                {
                    std::cout << "[ ";
                    for (size_t col = 0; col < _dim ; ++col)
                    {
                        std::cout << _2darray[row][col] <<' ';
                    }
                    std::cout <<"]"<<'\n';
                }
                std::cout <<'\n';
            }

        }

        template <typename T>
        void matrix<T>::ReadMtxFormat(std::string path)
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
                    return;
                }
                _dim = tmp_row;

                allocateMemForMatrix(); 

                while( read >> tmp_row)
                {
                    read >> tmp_col;
                    read >> tmp_val;
                    _2darray[tmp_row-1][tmp_col-1] = tmp_val;

                    if (is_sym && tmp_col != tmp_row ) 
                    {
                        _2darray[tmp_col-1][tmp_row-1] = tmp_val;
                    }
                }


            }
            else 
            {
                std::cout <<"not open " <<'\n';
            }
            read.close();
        }

        template <typename T>
        matrix<T>& matrix<T>:: transposeInPlace ()
        {
            for ( size_t i = 0; i < this->dim() - 1 ; ++i )
            {
                for ( size_t k = i + 1 ; k < this->dim() ; ++k)
                {
                    std::swap(_2darray[i][k], _2darray[k][i]);
                }
            }

            return *this;

        }
        template <typename T>
        matrix<T> matrix<T>::operator*( const matrix<T>& mat )
        {
            if ( mat._dim != this->_dim ) throw("matrix multiplication, dimensions does not match " );

            matrix<T> result(mat._dim);
            #if defined(par_exe_gcc) || defined(par_exe_msvc)
            #pragma omp parallel for schedule(static) num_threads(8)  
            #endif
            for ( size_t row = 0 ; row < mat._dim ; ++row)
            {
                for( size_t col = 0 ; col < mat._dim ; ++col)
                {
                    #if defined(par_exe_gcc) || defined(par_exe_msvc)
                    #pragma omp simd
                    #endif
                        for( size_t it = 0; it < mat._dim ; ++it)
                        {
                            result(row, col ) += ( this->_2darray[row][it] * mat._2darray[it][col] );
                        }
                }
            }
            return result;
        }

        template <typename T>
        vec<T> matrix<T>::inverseDiagElements( )
        {
            vec<T> result(this->_dim );

        #ifdef par_exe_gcc 
        #pragma omp parallel for
        for (size_t i = 0; i < _dim ; ++i)
        {
            result[i] = 1.0/(this->_2darray[i][i]);
        }

        #else
        std::for_each(result.begin(), result.end(), [this, i = 0]  (T & t) mutable {
                     t = 1.0/(this->_2darray[i][i]);
                    ++i;
                });
        #endif

        return result;


        }

        template<typename T>
        vec<T> matrix<T>::operator*( vec<T>& rhs)
        {
            if ( this->_dim != rhs.size() ) throw("dimensions dont fit broh");
            vec<T> result(rhs.size());
        #ifdef par_exe_gcc 
        #pragma omp parallel for
            for ( size_t row = 0; row < rhs.size() ; ++row )
            {
                result[row] = std::inner_product(rhs.begin(), rhs.end(), this->getRow(row).begin(),0.0);
            }
            return result;
        #else
            for ( size_t row = 0; row < rhs.size() ; ++row )
            {
                result[row] = std::inner_product(rhs.begin(), rhs.end(), this->getRow(row).begin(),0.0);
            }
            return result;
        #endif

        }


        template<typename T>
        inline vec<T> matrix<T>::getRow(size_t row ) const
        {
            return vec<T>(_2darray[row], _dim );
        }

        template<typename T>
        inline vec<T> matrix<T>::getRow(size_t row ) 
        {
            return vec<T>(_2darray[row], _dim );
        }

        template<typename T>
        matrix<T> matrix<T>::operator*( T t )
        {
            matrix<T> result(_dim);
            #ifdef par_exe_gcc 
            #pragma omp parallel for collapse(2)
            #endif
            for ( size_t row = 0; row < _dim ; ++row )
            {
                for ( size_t col = 0; col < _dim ; ++col )
                {
                    result(row,col) = t * this->_2darray[row][col];
                }
            }
            return result;
        }
        template<typename T>
        vec<T> matrix<T>::Diag()
        {
            vec<T> result(_dim);
            #ifdef par_exe_gcc 
            #pragma omp parallel for
            #endif
            for ( size_t i = 0 ; i < _dim ; ++i)
            {
                result[i] = this->_2darray[i][i];
            }
            return result;

        }

        template<typename T>
        matrix<T> matrix<T>::Lower()
        {
            matrix<T> result(_dim);
            #ifdef par_exe_gcc 
            #pragma omp parallel for
            #endif
            for (size_t i  = 0; i < _dim ; ++i)
            {
                for ( size_t j = 0 ; j <= i ; ++j)
                {
                    result(i,j) = this->_2darray[i][j];
                    
                }
            }
            return result;

        }

        template<typename T>
        matrix<T> matrix<T>::Upper()
        {
            matrix<T> result(_dim);
            #ifdef par_exe_gcc 
            #pragma omp parallel for
            #endif
            for (size_t i  = 0; i < _dim ; ++i)
            {
                for ( size_t j = 0 ; j >= i && j < _dim ; ++j)
                {
                    result(i,j) = this->_2darray[i][j];
                }
            }
            return result;

        }

        //TODO
     template<typename T>
     bool checkMatrixForJacobi( matrix<T>& mat)
     {
         std::vector<T> row_sum;
         row_sum.reserve(mat.dim());
         double tmp = 0.0;
         for(size_t i = 0; i < mat.dim() ; ++i)
         {
             for(size_t j = 0; j < mat.dim() ; ++j)
             {
                 if ( i != j) tmp += std::abs(mat(i,j)); 
             }
                double criteria = tmp / static_cast<double>(mat(i,i));       
            if ( criteria > 1.1)
            {
                std::cout <<"criteria = " << criteria << " is bigger than 1 " <<'\n';
                return false;
            }
            //else if(criteria < 1.000000000001 &&  criteria > 0.99999999999)
            //{
            //    std::cout <<" may or may not converge I guess " <<'\n';
            //    return false;
            //}
             
         }
         return true;

     }

    template<typename T>
    matrix<T> matrix<T>::minor(size_t s)
    {
        assert( this->_dim > 2 && "Dimension is smaller than two, no minor");
        assert( this->_dim > s && "Wrong index for minor ");
        matrix<T> result(_dim -1);
        #ifdef par_exe_gcc 
        #pragma omp parallel for
        #endif
        for( size_t i = 1 ; i < _dim ; ++i)
        {
            for ( size_t j = 0; j < _dim ;++j)
            {
                if ( j == s) continue;
                else if ( j < s) result(i - 1, j ) = this->_2darray[i][j];
                else if ( j > s) result(i - 1, j-1 ) = this->_2darray[i][j];
            }
        }
        return result;

    }


     matrix<double> generateSPDMatrix(size_t size);

    template <typename T>
    std::ostream& operator << ( std::ostream& os, const matrix<T>& mat )
    {
        for ( size_t i = 0; i < mat.dim() ; ++i)
        {
            for ( size_t j = 0 ; j < mat.dim(); ++j)
            {
                os <<std::setprecision(8)<< std::fixed << mat(i,j) << ' ';
            }
            os << '\n';
        }
        os << '\n';
            return os;
    };
    

	template<typename T>
	class SparseMatrix
	{

		public:

			// === CREATION ==============================================

			SparseMatrix(int n); // square matrix n×n
			SparseMatrix(int rows, int columns); // general matrix

			SparseMatrix(const SparseMatrix<T> & m); // copy constructor
			SparseMatrix<T> & operator = (const SparseMatrix<T> & m);

			~SparseMatrix(void);


			// === GETTERS / SETTERS ==============================================

			int getRowCount(void) const;
			int getColumnCount(void) const;
            std::vector<int>* getRows() const;
            std::vector<int>* getCols() const;
            std::vector<double>* getVals() const;

            std::vector<T> inverseDiagElements();

			// === READERS ============================================
            void ReadMtxFormat(std::string path);


			// === VALUES ==============================================

			T get(int row, int col) const;
			SparseMatrix & set(T val, int row, int col);


			// === OPERATIONS ==============================================

			std::vector<T> multiply(const std::vector<T> & x) const;
			std::vector<T> operator * (const std::vector<T> & x) const;

			SparseMatrix<T> multiply(const SparseMatrix<T> & m) const;
			SparseMatrix<T> operator * (const SparseMatrix<T> & m) const;

			SparseMatrix<T> add(const SparseMatrix<T> & m) const;
			SparseMatrix<T> operator + (const SparseMatrix<T> & m) const;

			SparseMatrix<T> subtract(const SparseMatrix<T> & m) const;
			SparseMatrix<T> operator - (const SparseMatrix<T> & m) const;


			// === FRIEND FUNCTIONS =========================================

			template<typename X>
			friend bool operator == (const SparseMatrix<X> & a, const SparseMatrix<X> & b);

			template<typename X>
			friend bool operator != (const SparseMatrix<X> & a, const SparseMatrix<X> & b);

			template<typename X>
			friend std::ostream & operator << (std::ostream & os, const SparseMatrix<X> & matrix);


		protected:

			int m_, n_;

			std::vector<T> * vals_;
            std::vector<int> * rows_, * cols_;


			// === HELPERS / VALIDATORS ==============================================

			void construct(int m, int n);
			void destruct(void);
			void deepCopy(const SparseMatrix<T> & m);
			void validateCoordinates(int row, int col) const;
			void insert(int index, int row, int col, T val);
			void remove(int index, int row);

	};

    template<typename T>
    SparseMatrix<T>::SparseMatrix(int n)
    {
    	this->construct(n, n);
    }
    
    
    template<typename T>
    SparseMatrix<T>::SparseMatrix(int rows, int columns)
    {
    	this->construct(rows, columns);
    }
    
    
    template<typename T>
    SparseMatrix<T>::SparseMatrix(const SparseMatrix<T> & matrix)
    {
    	this->deepCopy(matrix);
    }
    
    
    template<typename T>
    SparseMatrix<T> & SparseMatrix<T>::operator = (const SparseMatrix<T> & matrix)
    {
    	if (&matrix != this) {
    		this->destruct();
    		this->deepCopy(matrix);
    	}
    
    	return *this;
    }
    
    
    template<typename T>
    void SparseMatrix<T>::deepCopy(const SparseMatrix<T> & matrix)
    {
    	this->m_ = matrix.m_;
    	this->n_ = matrix.n_;
    	this->rows_ = new std::vector<int>(*(matrix.rows_));
    
    	if (matrix.vals_ != NULL) {
    		this->cols_ = new std::vector<int>(*(matrix.cols_));
    		this->vals_ = new std::vector<T>(*(matrix.vals_));
    	}
    }
    
    
    template<typename T>
    SparseMatrix<T>::~SparseMatrix(void)
    {
    	this->destruct();
    }
    
    
    template<typename T>
    void SparseMatrix<T>::construct(int rows, int columns)
    {
    	if (rows < 1 || columns < 1) {
    		throw ("Matrix dimensions cannot be zero or negative.");
    	}
    
    	this->m_ = rows;
    	this->n_ = columns;
    
    	this->vals_ = NULL;
    	this->cols_ = NULL;
    	this->rows_ = new std::vector<int>(rows + 1, 1);
    }
    
    
    template<typename T>
    void SparseMatrix<T>::destruct(void)
    {
    	if (this->vals_ != NULL) {
    		delete this->vals_;
    		delete this->cols_;
    	}
    
    	delete this->rows_;
    }
    
    
    // === GETTERS / SETTERS ==============================================
    
    template<typename T>
    inline int SparseMatrix<T>::getRowCount(void) const
    {
    	return this->m_;
    }
    
    
    template<typename T>
    inline int SparseMatrix<T>::getColumnCount(void) const
    {
    	return this->n_;
    }

    template<typename T>
    std::vector<int>* SparseMatrix<T>::getRows() const
    {
        return this->rows_;
    }

    template<typename T>
    std::vector<int>* SparseMatrix<T>::getCols() const
    {
        return this->cols_;
    }

    template<typename T>
    std::vector<double>* SparseMatrix<T>::getVals() const
    {
        return this->vals_;
    }

    
    
    // === VALUES ==============================================
    
    template<typename T>
    inline T SparseMatrix<T>::get(int row, int col) const
    {
    	this->validateCoordinates(row, col);
    
    	int currCol;
    
    	for (int pos = (*(this->rows_))[row - 1] - 1; pos < (*(this->rows_))[row] - 1; ++pos) {
    		currCol = (*(this->cols_))[pos];
    
    		if (currCol == col) {
    			return (*(this->vals_))[pos];
    
    		} else if (currCol > col) {
    			break;
    		}
    	}
    	return T();
    }
    
    
    template<typename T>
    SparseMatrix<T> & SparseMatrix<T>::set(T val, int row, int col)
    {
    	this->validateCoordinates(row, col);
    
    	int pos = (*(this->rows_))[row - 1] - 1;
    	int currCol = 0;
    
    	for (; pos < (*(this->rows_))[row] - 1; pos++) {
    		currCol = (*(this->cols_))[pos];
    
    		if (currCol >= col) {
    			break;
    		}
    	}
    
    	if (currCol != col) {
    		if (!(val == T())) {
    			this->insert(pos, row, col, val);
    		}
    
    	} else if (val == T()) {
    		this->remove(pos, row);
    
    	} else {
    		(*(this->vals_))[pos] = val;
    	}
    
    	return *this;
    }
    
    
    // === OPERATIONS ==============================================
    
    template<typename T>
    std::vector<T> SparseMatrix<T>::multiply(const std::vector<T> & x) const
    {
    	if (this->n_ != (int) x.size()) {
    		throw ("Cannot multiply: Matrix column count and vec size don't match.");
    	}
    
    	std::vector<T> result(this->m_, T());
    
    	if (this->vals_ != NULL) { // only if any value set
            #pragma omp parallel for  
    		for (int i = 0; i < this->m_; i++) {
    			T sum = T();
                #pragma omp simd 
    			for (int j = (*(this->rows_))[i]; j < (*(this->rows_))[i + 1]; j++) {
    				sum = sum + (*(this->vals_))[j - 1] * x[(*(this->cols_))[j - 1] - 1];
    			}
    
    			result[i] = sum;
    		}
    	}
    
    	return result;
    }
    
    
    template<typename T>
    std::vector<T> SparseMatrix<T>::operator * (const std::vector<T> & x) const
    {
    	return this->multiply(x);
    }
    
    
    template<typename T>
    SparseMatrix<T> SparseMatrix<T>::multiply(const SparseMatrix<T> & m) const
    {
    	if (this->n_ != m.m) {
    		throw ("Cannot multiply: Left matrix column count and right matrix row count don't match.");
    	}
    
    	SparseMatrix<T> result(this->m_, m.n);
    
    	T a;
    
    	// TODO: more efficient?
    	// @see http://www.math.tamu.edu/~srobertp/Courses/Math639_2014_Sp/CRSDescription/CRSStuff.pdf
    
    	for (int i = 1; i <= this->m_; i++) {
    		for (int j = 1; j <= m.n; j++) {
    			a = T();
    
    			for (int k = 1; k <= this->n_; k++) {
    				a = a + this->get(i, k) * m.get(k, j);
    			}
    
    			result.set(a, i, j);
    		}
    	}
    
    	return result;
    }
    
    
    template<typename T>
    SparseMatrix<T> SparseMatrix<T>::operator * (const SparseMatrix<T> & m) const
    {
    	return this->multiply(m);
    }
    
    
    template<typename T>
    SparseMatrix<T> SparseMatrix<T>::add(const SparseMatrix<T> & m) const
    {
    	if (this->m_ != m.m || this->n_ != m.n) {
    		throw ("Cannot add: matrices dimensions don't match.");
    	}
    
    	SparseMatrix<T> result(this->m_, this->n_);
    
    	// TODO: more efficient?
    	// @see http://www.math.tamu.edu/~srobertp/Courses/Math639_2014_Sp/CRSDescription/CRSStuff.pdf
    
    	for (int i = 1; i <= this->m_; i++) {
    		for (int j = 1; j <= this->n_; j++) {
    			result.set(this->get(i, j) + m.get(i, j), i, j);
    		}
    	}
    
    	return result;
    }
    
    
    template<typename T>
    SparseMatrix<T> SparseMatrix<T>::operator + (const SparseMatrix<T> & m) const
    {
    	return this->add(m);
    }
    
    
    template<typename T>
    SparseMatrix<T> SparseMatrix<T>::subtract(const SparseMatrix<T> & m) const
    {
    	if (this->m_ != m.m || this->n_ != m.n) {
    		throw ("Cannot subtract: matrices dimensions don't match.");
    	}
    
    	SparseMatrix<T> result(this->m_, this->n_);
    
    	// TODO: more efficient?
    	// @see http://www.math.tamu.edu/~srobertp/Courses/Math639_2014_Sp/CRSDescription/CRSStuff.pdf
    
    	for (int i = 1; i <= this->m_; i++) {
    		for (int j = 1; j <= this->n_; j++) {
    			result.set(this->get(i, j) - m.get(i, j), i, j);
    		}
    	}
    
    	return result;
    }
    
    
    template<typename T>
    SparseMatrix<T> SparseMatrix<T>::operator - (const SparseMatrix<T> & m) const
    {
    	return this->subtract(m);
    }
    
    
    // === HELPERS / VALIDATORS ==============================================
    
    template<typename T>
    inline void SparseMatrix<T>::validateCoordinates(int row, int col) const
    {
    	if (row < 1 || col < 1 || row > this->m_ || col > this->n_) {
    		throw ("Coordinates out of range.");
    	}
    }
    
    
    template<typename T>
    void SparseMatrix<T>::insert(int index, int row, int col, T val)
    {
    	if (this->vals_ == NULL) {
    		this->vals_ = new std::vector<T>(1, val);
    		this->cols_ = new std::vector<int>(1, col);
    
    	} else {
    		this->vals_->insert(this->vals_->begin() + index, val);
    		this->cols_->insert(this->cols_->begin() + index, col);
    	}
    
    	for (int i = row; i <= this->m_; i++) {
    		(*(this->rows_))[i] += 1;
    	}
    }
    
    
    template<typename T>
    void SparseMatrix<T>::remove(int index, int row)
    {
    	this->vals_->erase(this->vals_->begin() + index);
    	this->cols_->erase(this->cols_->begin() + index);
    
    	for (int i = row; i <= this->m_; i++) {
    		(*(this->rows_))[i] -= 1;
    	}
    }

    template<typename T>
    void SparseMatrix<T>::ReadMtxFormat(std::string path)
    {
            this->destruct();
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
                    return;
                }

                this->n_ = tmp_row;

                this->construct(this->n_,this->n_);

                while( read >> tmp_row)
                {
                    read >> tmp_col;
                    read >> tmp_val;

                    this->validateCoordinates(tmp_row, tmp_col);

                    this->set(tmp_val, tmp_row, tmp_col );

                    if ( is_sym && tmp_col != tmp_row ) 
                    {
                        this->validateCoordinates(tmp_col, tmp_row);
                        this->set(tmp_val, tmp_col, tmp_row );
                    }
                }


            }
            else 
            {
                std::cout <<"not open " <<'\n';
            }
            read.close();

    }

    template<typename T>
    std::vector<T> SparseMatrix<T>::inverseDiagElements()
    {
        std::vector<T> result(this->n_);
        for(int i = 0; i < this->n_ ; ++i)
        {
            result[i] = 1.0/(this->get(i+1,i+1));
        }
        return result;

        
    }

    
    
    // === FRIEND FUNCTIONS =========================================
    
    template<typename T>
    bool operator == (const SparseMatrix<T> & a, const SparseMatrix<T> & b)
    {
    	return ((a.vals_ == NULL && b.vals_ == NULL)
    				|| (a.vals_ != NULL && b.vals_ != NULL && *(a.vals_) == *(b.vals_)))
    			&& ((a.cols_ == NULL && b.cols_ == NULL)
    				|| (a.cols_ != NULL && b.cols_ != NULL && *(a.cols_) == *(b.cols_)))
    			&& *(a.rows_) == *(b.rows_);
    }
    
    
    template<typename T>
    bool operator != (const SparseMatrix<T> & a, const SparseMatrix<T> & b)
    {
    	return !(a == b);
    }
    
    
    template<typename T>
    std::ostream & operator << (std::ostream & os, const SparseMatrix<T> & matrix)
    {
    	for (int i = 1; i <= matrix.m_; i++) {
    		for (int j = 1; j <= matrix.n_; j++) {
    			if (j != 1) {
    				os << " ";
    			}
    
    			os << matrix.get(i, j);
    		}
    
    		if (i < matrix.m_) {
    			os << '\n';
    		}
    	}
    
    	return os;
    }

        
}
#endif

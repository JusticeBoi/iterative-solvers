#ifndef __vector__hpp
#define __vector__hpp

#include <iostream>
#include <cstring>
#include <cmath>
#include "parallelDefs.hpp"
#include <functional>
#include <fstream>
#include <iomanip>
//#if defined(_MSC_VER) && __cplusplus >= 201500
//    #pragma message "parallel active , need -fopenmp flag"
//    #include <execution>
//    #define par_exe_msvc
//#elif defined(__GNUG__) && __cplusplus >= 201500
//    #pragma message "parallel active , need -fopenmp flag"
//    #define par_exe_gcc
//    #include <parallel/algorithm>
//    #include <parallel/numeric>
//#else
//    #pragma message "parallel deactive , need -fopenmp flag and c++17 for parallel execution"
//    #include <algorithm>
//    #include <numeric>
//#endif


//#define par_exe_gcc
namespace linalg
{
    template <typename T>
    class matrix;

    template<typename T>
    class vec
    {
        private:

            T* _array;
            size_t _size;
            size_t _capacity;

            template<typename T1>
            friend void swap( vec<T1>& lhs, vec<T1>& rhs);

            template<typename T1>
            friend vec<T1> operator*( vec<T1>& rhs);

            template<typename T1>
            friend bool operator== ( const vec<T1>& lhs, const vec<T1>& rhs );

            template<typename T1>
            friend std::ostream& operator<<( std::ostream& os, vec<T1>& v );

        public:
            /*constuctor */
            vec();
            vec(std::initializer_list<T> init);
            vec(T* arr, size_t size);
            vec(size_t size);
            vec(size_t size, T initialValue );
            vec(const vec& _rhs);
            vec(vec && ) noexcept;

            /*assignments*/
            vec& operator=( vec rhs );
            vec& operator=(std::initializer_list<T> reinit );
            
            /*operators*/
            T& operator[](size_t loc);
            const T& operator[](size_t loc) const;
            
            void assignValue( T val);
            void print() const;
            /*destruction*/
            void reset();
            ~vec();

            void reserve(size_t );
            void push_back(const T& pb);

            /*iterators*/
            T* begin();
            const T* begin() const;
            T* end();
            const T* end() const ;

            T at ( size_t index );

            void write() const;
            std::vector<T> toStd();

            /*operations */
            T operator*(const vec&);
            vec<T> operator*(const T&);
            vec<T>& operator+= ( const vec<T>& rhs );
            vec<T>& operator-= ( const vec<T>& rhs );
            vec<T> operator+ ( const vec<T>& rhs );
            vec<T> operator- ( const vec<T>& rhs );
            double normL2()const ; 

            /* getters */
            size_t size() const ;
            size_t capacity() const ;
    };

    /*default constructor */
    template<typename T>
    vec<T>::vec():_array(nullptr),_size(0), _capacity(0)
    {
    }


    /* reset */
    template< typename T>
    void vec<T>::reset()
    {
        delete[] _array;
        _size = 0;
        _capacity = 0;
        _array = nullptr;
    }

    /*constructor */
    template<typename T>
    vec<T>::vec(size_t size):_size(size), _capacity(size)
    {
        _array = new T[size];         
        assignValue(0);
    }
    
    /* constructo from array */
    template <typename T>
    vec<T>::vec(T* arr, size_t size):_size(size),_capacity(size)
    {
        _array = new T[_size];
        std::memcpy(_array, arr, _size * sizeof(T));

    }

    /*constructor */
    template<typename T>
    vec<T>::vec(size_t size, T initialValue ):_size(size),_capacity(size)
    {
       _array = new T[size];
        assignValue(initialValue); 
    }
   

    /*copy constructor */
   template <typename T>
    vec<T>::vec(const vec& rhs): _size(rhs._size), _capacity(rhs._capacity)
    {
        _array = new T[_size];
        //std::memcpy(_array, rhs._array, _size * sizeof(T));
        std::copy(rhs.begin(), rhs.end(), this->begin() ); 
        //for (size_t i = 0; i < _size ; ++i) _array[i] = rhs._array[i]; 
    }
   
    /*Destructor */
    template <typename T>
    vec<T>::~vec()
    {
        reset();
    }

    /*move constructor */
    template<typename T>
    vec<T>::vec(vec && rhs) noexcept:_size(rhs._size), _capacity(rhs._capacity) 
    {
        #ifdef show_moves
        std::cout <<" vec move constructor " <<'\n';
        #endif
        _array = std::move(rhs._array);
        rhs._array = nullptr;
        rhs._size = 0;
        rhs.reset();
    }

    /*copy assignment */
    //here rhs passed by value, calls copy constructor
    template <typename T>
    vec<T>& vec<T>::operator=(vec rhs )
    {
        /*copy and swap */
        swap(*this, rhs);
        return *this;
    }

    /*initializer list  assignment */
    template <typename T>
    vec<T>& vec<T>::operator=(std::initializer_list<T> reinit )
    {
        reset();
        _capacity = reinit.size();
        _size = reinit.size();
        _array = new T[_capacity];
        for ( size_t i = 0 ; i < _size ; ++i) _array[i] = reinit.begin()[i]; 
        //for(T i : reinit) this->push_back(i);

        return *this;
        

    }

    /*assign value */
    template <typename T >
    inline void vec<T>::assignValue(T val)
    {
        std::fill(this->begin(), this->end(), val );
    }


    /* swap */
    template <typename T>
    inline void swap(vec<T>& lhs, vec<T>& rhs)
    {
        using std::swap;
        swap(lhs._size, rhs._size);
        swap(lhs._capacity, rhs._capacity);
        swap(lhs._array, rhs._array);
    }

    template<typename T>
    void vec<T>::push_back(const T& pb)
    {
       if ( _capacity == _size ) 
       {
            #ifdef debug
           std::cout <<"inefficient stuff " << '\n';
            #endif
            T* tmp = _array;
            _array = new T[_capacity + 5];

            std::memcpy(_array, tmp, _size * sizeof(T) );            

            _capacity += 5;
            _array[_size] = pb;
            ++_size;

            delete[] tmp;
       }
       else
       {
            _array[_size] = pb;
            ++_size;
       }

    }

    /*initilizer list construction */
    template<typename T>
    vec<T>::vec(std::initializer_list<T> init ):_size(init.size()), _capacity(init.size())
    {
        _array = new T[_capacity];
        for ( size_t i = 0 ; i < _size ; ++i) _array[i] = init.begin()[i]; 
        //for(auto i : init) this->push_back(i);
    }


    template <typename T>
    inline T* vec<T>::begin()
    {
        return _array;
    }

    template <typename T>
    inline const T* vec<T>::begin() const 
    {
        return _array;
    }

    template <typename T>
    inline T* vec<T>::end()
    {
        return _array + _size;
    }

    template <typename T>
    inline const T* vec<T>::end() const
    {
        return _array + _size;
    }


    template<typename T>
    T vec<T>::at(size_t index)
    {
        if ( index > _size  || static_cast<int>(index) < 0 ) std::cout << " out of bounds, wr " <<'\n';
        return _array[index]; 
    }

    template<typename T>
    inline T& vec<T>::operator[](size_t loc)
    {
        return _array[loc];

    }

    template<typename T>
    inline const T& vec<T>::operator[](size_t loc) const
    {
        return _array[loc];

    }


    /*reserve*/
    template <typename T>
    void vec<T>::reserve(size_t new_cap )
    {
        if( _capacity < new_cap )
        {
            T* tmp = _array;

            _array = new T[new_cap];
            _capacity = new_cap;

            //copy(tmp, _array, _size );            
            std::memcpy(_array, tmp, _size * sizeof(T) );            

            delete[] tmp;

        }
        else
        {
            /*TODO*/

        }


    }

    /*size*/
    template<typename T>
    inline size_t vec<T>::size() const
    {
        return _size;
    }


    /*capacity*/
    template<typename T>
    inline size_t vec<T>::capacity() const
    {
        return _capacity;
    }

    /*dot product*/
    template <typename T>
    T vec<T>::operator*(const vec& rhs)
    {
        if ( _size != rhs._size )
        {
            throw("sizes of the vectors are not the same");
        }
        else
        {
        #ifdef par_exe_gcc 
        return __gnu_parallel::inner_product(this->begin(), this->end(), rhs.begin(), 0.0);
        #else
        return std::inner_product(this->begin(), this->end(), rhs.begin(), 0.0);
        #endif 
        }

    }

    template <typename T>
    inline vec<T> vec<T>::operator*(const T& rhs)
    {
        if ( !rhs ) return vec<T>(this->size());
        #ifdef par_exe_gcc 
        vec<T> result(this->_size);
        __gnu_parallel::transform(this->begin(), this->end(), result.begin(),
             std::bind(std::multiplies<T>(), std::placeholders::_1, rhs)); 
        return result;
        #else
        vec<T> result(this->_size);
        std::transform(this->begin(), this->end(), result.begin(),
             std::bind(std::multiplies<T>(), std::placeholders::_1, rhs)); 
        return result;
        #endif 
        //#ifdef par_exe_gcc 
        //if ( !rhs ) return vec<T>(this->size());

        //vec<T> result(*this);
        //__gnu_parallel::for_each(result.begin(), result.end(), [&rhs](auto& t)
        //        {
        //          t *=rhs;  
        //        });
        //return result;
        //#else
        //if ( !rhs ) return vec<T>(this->size());
        //vec<T> result(*this);
        //std::for_each(result.begin(), result.end(), [&rhs](auto& t)
        //        {
        //          t *= rhs;  
        //        });
        //return result;
        //#endif 

    }

    /* normL2 */
    template <typename T>
    double vec<T>::normL2() const
    {
        #ifdef par_exe_gcc 
        return std::sqrt(__gnu_parallel::inner_product(this->begin(), this->end(), this->begin(), 0.0));
        #else
        return std::sqrt(std::inner_product(this->begin(), this->end(), this->begin(), 0.0));
        #endif 
    }

    /*print */
    template <typename T>
    void vec<T>::print() const
    {
        if(_size != 0 )
        {
            std::cout <<"[ ";
            for (size_t i = 0; i < _size ; ++i) std::cout<<std::setprecision(10) <<_array[i] <<' ' ;
            std::cout <<"]\n";
            std::cout <<'\n';
        }
        else 
        {
            std::cerr <<"cant show the vector coz dimension is 0 " <<'\n';
        }


    }
    
    //friends
    template<typename T>
    inline std::ostream& operator<<( std::ostream& os, vec<T>& v )
    {

        std::for_each(v.begin(), v.end(),[&os](T& t){ os << t<<'\n'; });
        return os;
    }

    template <typename T>
    bool operator== ( const vec<T>& lhs, const vec<T>& rhs)
    {
        if ( lhs.size() != rhs.size()) return false;
        else if ( std::memcmp(lhs._array, rhs._array, sizeof(lhs._array))) return false;
        //else if ( std::copy(lhs.begin(),lhs.end(), rhs.begin())) return false;
        else return true;
    }

    template<typename T>
    inline vec<T> vec<T>::operator+ ( const vec<T>& rhs )
    {
        if ( this->size() != rhs.size())
        {
            throw("error");
            return *this;
        }
        else 
        {
            vec<T> result(this->size());
            #ifdef par_exe_gcc 
            __gnu_parallel::transform( this->begin(), this->end(), rhs.begin(), result.begin(), std::plus<T>()); 
            #else
            std::transform( this->begin(), this->end(), rhs.begin(), result.begin(), std::plus<T>()); 
            #endif
            return result;
    
            //vec<T> result(this->size());
            //#ifdef par_exe_gcc 
            //__gnu_parallel::transform( this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<T>()); 
            //#else
            //std::transform( this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<T>()); 
            //#endif
            //return *this;
        }
    }
    template<typename T>
    inline vec<T> vec<T>::operator- ( const vec<T>& rhs )
    {
        if ( this->size() != rhs.size())
        {
            throw("error");
            return *this;
        }
        else 
        {
            vec<T> result(this->size());
            #ifdef par_exe_gcc 
            __gnu_parallel::transform( this->begin(), this->end(), rhs.begin(), result.begin(), std::minus<T>()); 
            #else
            std::transform( this->begin(), this->end(), rhs.begin(), result.begin(), std::minus<T>()); 
            #endif
            return result;
    
            //vec<T> result(this->size());
            //#ifdef par_exe_gcc 
            //__gnu_parallel::transform( this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<T>()); 
            //#else
            //std::transform( this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<T>()); 
            //#endif
            //return *this;
        }
    }

    template<typename T>
    inline vec<T>& vec<T>::operator+= ( const vec<T>& rhs )
    {
        if ( this->size() != rhs.size())
        {
            return *this;
        }
        else 
        {
            #ifdef par_exe_gcc 
            __gnu_parallel::transform( this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<T>()); 
            #else
            std::transform( this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<T>()); 
            #endif
            return *this;
    
        }
    }

    template<typename T>
    inline vec<T>& vec<T>::operator-= ( const vec<T>& rhs )
    {
        if ( this->size() != rhs.size())
        {
            return *this;
        }
        else 
        {
            #ifdef par_exe_gcc 
             /*The binary operation std::minus is applied to pairs of elements from two ranges:
              * one defined by [first1, last1) and the other beginning at first2.*/

            __gnu_parallel::transform( this->begin(), this->end(), rhs.begin(), this->begin(), std::minus<T>()); 
            #else
            std::transform( this->begin(), this->end(), rhs.begin(), this->begin(), std::minus<T>()); 
            #endif
            return *this;
    
        }
    }
    template<typename T>
    void vec<T>::write() const
    {
        std::string name = "vector_of " + std::to_string(_size);
        std::ofstream of(name);
        std::for_each(this->begin(), this->end(), [&of,this](const auto val)
                {
                    of << val <<'\n';
                });
        of.close();
    }

    template<typename T>
    std::vector<T> vec<T>::toStd()
    {
        std::vector<T> result(_size, T());
        std::copy(this->begin(), this->end(), result.begin());
        return result;
    }



}
#endif

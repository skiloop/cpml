#ifndef DATA1D_H
#define DATA1D_H

#include <string>

template<class T>
class data1d {
public:
    /**
     * constructor
     * @param num array length
     * @param val initial variable
     */
    data1d(unsigned num, T val = 0);

    /**
     * default donstructor
     * set @c p NULL
     * array length zero
     */
    data1d() : p(NULL), n(0) {
    };

    /**
     * copy constructor
     * @param orig
     */
    data1d(const data1d& orig) ;

    /**
     * deconstructor
     */
    ~data1d();

    /**
     * create array with length @c num
     * @param num
     */
    void createArray(unsigned num) ;
    
    /**
     * create array with length @c num and initial array with @c val
     * @param num
     * @param val
     */
    void createArray(unsigned num, T val) ;

    /**
     * initial array to @c initval
     * @param initval
     */
    void initArray(T initval = 0) ;

    /**
     * reset array to zero
     */
    void resetArray() ;
    
    /**
     * save array with file name @c name
     * @param name
     */
    void save(const std::string name) ;
public:
    T* p;
    unsigned n;
};

#include "data1d.hpp"

#endif // DATA1D_H

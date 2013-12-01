#ifndef DATA1D_H
#define DATA1D_H

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;

template<typename T>
class data1d {
public:

    /**
     * constructor
     * @param num array length
     * @param val initial variable
     */
    data1d(unsigned num, T val = 0) : p(NULL), n(num) {
        createArray(num);
        initArray(val);
    };

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
    data1d(const data1d& orig) {
        if (orig.p != NULL && orig.n > 0) {
            p = new T[orig.n];
            for (unsigned i = 0; i < orig.n; i++)p[i] = orig.p[i];
            this->n = orig.n;
        } else if (orig.n != 0) {
            createArray(orig.n, 0);
            this->n = orig.n;
        }
    };

    /**
     * deconstructor
     */
    ~data1d() {
        if (p != NULL)delete []p;
    };

    /**
     * create array with length @c num
     * @param num
     */
    void createArray(unsigned num) {
        if (num > 0) {
            if(NULL!=p){
                delete []p;
            }
            p = new T[num];
            n = num;
        }
    };
    
    /**
     * create array with length @c num and initial array with @c val
     * @param num
     * @param val
     */
    void createArray(unsigned num, T val) {
        createArray(num);
        initArray(val);
    };

    /**
     * initial array to @c initval
     * @param initval
     */
    void initArray(T initval = 0) {
        if (p == NULL)return;
        for (unsigned i = 0; i < n; i++)p[i] = initval;
    };

    /**
     * reset array to zero
     */
    void resetArray() {
        initArray();
    };
    
    /**
     * save array with file name @c name
     * @param name
     */
    void save(const string name) {
        ofstream out;
        out.open(name.c_str());
        if (out.is_open()) {
            for (int i = 0; i < n; i++) {
                out.setf(ios::fixed);
                out << p[i] << endl;
            }
            out.close();
        }
    };
public:
    T* p;
    unsigned n;
};

#endif // DATA1D_H

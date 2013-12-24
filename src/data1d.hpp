/**
 * class template data1d implementations
 * 
 *  
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;

template<class T>
data1d<T>::data1d(unsigned num, T val): p(NULL), n(num) {
    createArray(num);
    initArray(val);
}

template<class T>
data1d<T>::data1d(const data1d<T>& orig) {
    if (orig.p != NULL && orig.n > 0) {
        p = new T[orig.n];
        for (unsigned i = 0; i < orig.n; i++)p[i] = orig.p[i];
        this->n = orig.n;
    } else if (orig.n != 0) {
        createArray(orig.n, 0);
        this->n = orig.n;
    }
}

template<class T>
data1d<T>::~data1d(){
    if (p != NULL)delete []p;
}

template<class T>
void data1d<T>::createArray(unsigned num) {
    if (num > 0) {
        if (NULL != p) {
            delete []p;
        }
        p = new T[num];
        n = num;
    }
}

template<class T>
void data1d<T>::createArray(unsigned num, T val) {
    createArray(num);
    initArray(val);
};

template<class T>
void data1d<T>::initArray(T initval) {
    if (p == NULL)return;
    for (unsigned i = 0; i < n; i++)p[i] = initval;
}

template<class T>
void data1d<T>::resetArray() {
    initArray();
}

template<class T>
void data1d<T>::save(const string name) {
    ofstream out;
    out.open(name.c_str());
    if (out.is_open()) {
        for (int i = 0; i < n; i++) {
            out.setf(ios::fixed);
            out << p[i] << endl;
        }
        out.close();
    }
}


/* 
 * File:   data2d.cpp
 * Author: skiloop
 * 
 * Created on 2013年12月17日, 下午5:08
 */


#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>

#include "data2d.h"

using namespace std;

template<class DataType> const string data2d<DataType>::OUTPUT_FILE_NAME_TAIL = ".dat";

template<class DataType> data2d<DataType>::data2d(const data2d< DataType >& obj) : p(NULL) {
    create2DArray(obj);
    try {
        unsigned i;
        for (i = 0; i < obj.nx; i++) {
            memcpy(p[i], obj.p[i], obj.ny * sizeof (DataType));
        }
    } catch (exception & e) {
        cerr << e.what() << endl;
    }
}

template<class DataType>
data2d<DataType>::~data2d() {
    freeArray();
}

template<class DataType>
int data2d<DataType>::create2DArray(unsigned nnx, unsigned nny) {
#if(DEBUG>=6)
    cout << "(" << nnx << ',' << nny << ')' << endl;
#endif
    if (p != NULL) {
        cerr << "pointer not NULL" << endl;
        return -1;
    }
    unsigned i, j;
    if (nnx == 0 || nny == 0 || nnx > MAX_ARRAY_SIZE || nny > MAX_ARRAY_SIZE) {
        cerr << "Invalid size!" << endl;
        return -1;
    }
    try {

        p = new DataType*[nnx];
        if (p == NULL) {
            cout << "Failed to create space for data2d!" << endl;
            return -1;
        }
        for (i = 0; i < nnx; i++) {
            p[i] = new DataType[nny];
        }

        nx = nnx;
        ny = nny;
    } catch (exception & e) {
        cerr << "Exception in " << __FUNCTION__ << e.what() << endl;
        exit(-1);
    }
#if(DEBUG>=9)
    cout << "create successfully" << endl;
#endif
    return 0;
}

template<class DataType>
int data2d<DataType>::create2DArray(unsigned nnx, unsigned nny, DataType initVal) {
    if (create2DArray(nnx, nny) < 0)
        return -1;
    return resetArray(initVal);
}

template<class DataType>
int data2d<DataType>::resetArray(DataType Val) {
    unsigned i, j;
    if (!checkArray())
        return -1;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            p[i][j] = Val;
        }
    }
    return 0;
}

template<class DataType>
void data2d<DataType>::printArray() {
    unsigned i, j;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            cout << p[i][j] << '\t';
        }
        cout << endl;
    }
}

template<class DataType>
bool data2d<DataType>::checkArray() {
    if (ny <= 0 || nx <= 0 || p == NULL) {
        return false;
    }
    return true;
}

template<class DataType>
void data2d<DataType>::saveArrayData(const unsigned num, unsigned leap) {
    unsigned i, j;

    stringstream ss;
    ss << mName << num << OUTPUT_FILE_NAME_TAIL;
    ofstream ofile(ss.str().c_str());

    if (!ofile.is_open()) {
        cerr << "Cannot open" << ss << endl;
        return;
    }
    //check leap
    if (leap > nx || leap > ny)leap = 0;

    leap += 1;
    for (i = 0; i < nx; i += leap) {
        for (j = 0; j < ny; j += leap) {
            ofile << p[i][j] << '\t';
        }
        ofile << endl;
    }
    ofile.close();
}

template<class DataType>
void data2d<DataType>::operator =(data2d< DataType > const &other) {
    unsigned i, j;
    if (&other == this || p == NULL || other.p == NULL ||
            other.nx != nx || other.ny != ny)
        return;
    for (i = 0; i < nx; i++) {
        memcpy(p[i][j], other.p[i][j], ny * sizeof (DataType));
    }
}

template<class DataType>
void data2d<DataType>::initArray(DataType initVal) {
    unsigned i, j;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            p[i][j] = initVal;
        }
    }
}

template<class DataType>
void data2d<DataType>::saveData(unsigned leap, unsigned step) {
    stringstream ss;
    ss << mName << "_" << step << OUTPUT_FILE_NAME_TAIL;
    string fname = ss.str();
    if (leap >= nx || leap >= ny) {
        cerr << "Invalid leap for saving p!" << endl;
        return;
    }
    ofstream out(fname.c_str(), ios_base::binary);
    if (!out.is_open()) {
        cerr << "File " << fname << "cannot be opened!" << endl;
        return;
    }
    if (leap <= 0)leap = 1;
    unsigned i, j;
    for (i = 0; i < nx; i += leap) {
        for (j = 0; j < ny; j += leap) {
            out << p[i][j] << '\t';
        }
        out << endl;
    }
    out.close();
}

/**
 *
 * @param k
 * @param leap
 */
template<class DataType>
void data2d<DataType>::save(unsigned k, unsigned leap, unsigned step, int type) {
    savePlain(k, leap, step, type);
}

template<class DataType>
void data2d<DataType>::save(int leap) {
    stringstream ss;
    ss << mName << OUTPUT_FILE_NAME_TAIL;
    string fname = ss.str();

    ofstream out(fname.c_str(), ios_base::binary);
    if (!out.is_open()) {
        cerr << "File " << fname << "cannot be opened!" << endl;
        return;
    }
    if (leap <= 0)leap = 1;
    unsigned i, j;
    for (i = 0; i < nx; i += leap) {
        for (j = 0; j < ny; j += leap) {
            out << p[i][j] << '\t';
        }
        out << endl;
    }

    out.close();
}

/**
 *
 * @param k
 * @param leap
 */
template<class DataType>
void data2d<DataType>::savePlain(unsigned k, unsigned leap, unsigned step, int type = 1) {
    switch (type) {
        case 1:
            saveXPlain(k, leap, step);
            break;
        case 2:
            saveYPlain(k, leap, step);
            break;
        default:
            saveXPlain(k, leap, step);
    }
}

/**
 *
 * @param k
 * @param leap
 * @param step
 */
template<class DataType>
void data2d<DataType>::saveYPlain(unsigned k, unsigned leap, unsigned step) {
    stringstream ss;
    ss << mName << "_y_" << step << OUTPUT_FILE_NAME_TAIL;
    if (leap >= nx || leap >= ny) {
        cerr << "Invalid leap for saving p!" << endl;
        return;
    }
    string fname = ss.str();
    ofstream out(fname.c_str(), ios_base::binary);
    if (!out.is_open()) {
        cerr << "File " << fname << "cannot be opened!" << endl;
        return;
    }
    if (leap <= 0)leap = 1;
    unsigned i;
    if (k >= ny)k = ny / 2;
    for (i = 0; i < nx; i += leap) {
        out << p[i][k] << '\t';
    }
    out.close();
}

/**
 *
 * @param k
 * @param leap
 * @param step
 */
template<class DataType>
void data2d<DataType>::saveXPlain(unsigned k, unsigned leap, unsigned step) {
    stringstream ss;
    ss << mName << "_x_" << step << OUTPUT_FILE_NAME_TAIL;
    if (leap >= nx || leap >= ny) {
        cerr << "Invalid leap for saving p!" << endl;
        return;
    }
    string fname = ss.str();
    ofstream out(fname.c_str(), ios_base::binary);
    if (!out.is_open()) {
        cerr << "File " << fname << "cannot be opened!" << endl;
        return;
    }
    if (leap <= 0)leap = 1;
    unsigned j;
    if (k >= nx)k = nx / 2;
    for (j = 0; j < ny; j += leap) {
        out << p[k][j] << '\t';
    }
    out.close();
}

template<class DataType>
int data2d<DataType>::create2DArray(const data2d< DataType > &stru) {
    return create2DArray(stru.nx, stru.ny);
}

template<class DataType>
int data2d<DataType>::create2DArray(const data2d< DataType > &stru, DataType initVal) {
    return create2DArray(stru.nx, stru.ny, initVal);
}

/**
 * check if is nan at point (i,j,k)
 *
 */
template<class DataType>
bool data2d<DataType>::isNaN(unsigned i, unsigned j) {
    if (isnan(p[i][j])) {
        cout << "nan var found for " << getName() << " at:(" << i << "," << j << ")" << endl;
        return true;
    }
    return false;
}

/**
 * check if is inf at point (i,j,k)
 *
 */
template<class DataType>
bool data2d<DataType>::isInf(unsigned i, unsigned j) {
    if (isinf(p[i][j])) {
        cout << "inf var found for " << getName() << " at:(" << i << "," << j << ")" << endl;
        return true;
    }
    return false;
}

/**
 * check if is valid var at point (i,j,k)
 */
template<class DataType>
bool data2d<DataType>::isValid(unsigned i, unsigned j) {
    return isNaN(i, j) | isInf(i, j);
}

template<class DataType>
void data2d<DataType>::whenLargerThan(unsigned i, unsigned j, MyDataF limit, void (*fun)()) {
    if (p[i][j] > limit) {
        //cout << "larger var found for " << getName() << " at:(" << i << "," << j << "," << k << ")" << endl;
        if (NULL != fun) {
            (*fun)();
        }
    }
}

template<class DataType>
void data2d<DataType>::freeArray() {
    if (p != NULL) {
        for (unsigned i = 0; i < nx; i++) {
            if (p[i] != NULL) {
                delete [] p[i];
            }
        }
        delete []p;
    }
}

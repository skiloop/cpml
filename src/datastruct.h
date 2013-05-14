/*
<one line to give the program's name and a brief idea of what it does.>
Copyright (C) 2011  <copyright holder> <email>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file datastruct.h
 * @version 0.0.0
 * @author skiloop ( skiloop@126.com )
 * @date 31/08/2011 0.0.0 created, by skiloop
 */
#ifndef DATASTRUCT_H
#define DATASTRUCT_H

#ifdef MATLAB_SIMULATION

#include <engine.h>
#include <mex.h>
#ifdef printf
#undef printf
#endif

#endif// end MATLAB_SIMULATION
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;
//#include "microdef.h"
#define MAX_ARRAY_SIZE 300000
//#include "common.h"

template<typename T>
class data1d {
public:

    data1d(unsigned num, T val = 0) : p(NULL), n(num) {
        createArray(num);
        initArray(val);
    };

    data1d() : p(NULL), n(0) {
    };

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

    ~data1d() {
        if (p != NULL)delete []p;
        p = NULL;
        n = 0;
    };

    void createArray(unsigned num) {
        if (num > 0) {
            p = new T[num];
            n = num;
        }
    };

    void CreateStruct(unsigned num) {
        createArray(num);
    };

    void initArray(T initval = 0) {
        if (p == NULL)return;
        for (unsigned i = 0; i < n; i++)p[i] = initval;
    };

    void resetArray() {
        initArray();
    };

    void createArray(unsigned num, T val) {
        createArray(num);
        initArray(val);
    };

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

/**
 * nx:point count in x direction
 * ny:point count in x direction
 * p:pointer to store data;
 */
template<class DataType>
class data3d {
public:
    unsigned int nx;
    unsigned int ny;
    unsigned int nz;
    DataType*** p;

public:
    static unsigned int cnt;
    static std::string tail;
#ifdef MATLAB_SIMULATION
    static Engine *ep;
private:
    static bool isMatlabEngineStarted;
#endif
private:
    std::string name;

private:
    int Number;
#ifdef MATLAB_SIMULATION
    mxArray *num;
    mxArray *MyArray;
#endif
public:

    /**
     * if cx==0 then p = NULL
     * else p has cx pointers;
     * 
     * if cy == 0 then all cx pointers of p is NULL
     * else p[i] has cy pointers with i from 0 to cx-1;
     * 
     * when cannot create space for p and p[i],exit program;
     */
    data3d(unsigned int cx, unsigned int cy, unsigned cz)
    : nx(cx), ny(cy), nz(cz), p(NULL) {
        unsigned i, j;
        if (cx == 0 || cy == 0 || cz == 0) {
            return;
        }
        try {

            p = new DataType**[cx];
            for (i = 0; i < cx; i++) {
                p[i] = new DataType*[cy];
            }
            for (i = 0; i < cx; i++) {
                for (j = 0; j < cy; j++) {
                    p[i][j] = new DataType[cz];
                }
            }
        } catch (exception & e) {
            cerr << e.what() << endl;
            return;
        }

    };

    data3d() : nx(0), ny(0), nz(0), p(NULL) {
    };
    data3d(const data3d< DataType > &obj);
    ~data3d();
    /**
     * print p in struct data3d
     */
    void PrintData();
    /**
     * free space created for data3d @c mst
     */
    void FreeStructData();

    /**
     * Set Data to val
     */
    int ResetStructData(DataType val = 0);

    /**
     * check p of data3d @c mst is valid
     * if p is not NULL and none of its subpointers,then 
     * return true,otherwise false
     */
    bool CheckStruct();

    /**
     * Create Space for struct data3d and initialize its @c nx and @c ny
     */
    int CreateStruct(unsigned nnx, unsigned nny, unsigned nnz);
    /**
     * Create Space for struct data3d and initialize its @c nx and @c ny
     */
    int CreateStruct(unsigned nnx, unsigned nny, unsigned nnz, DataType initVal);

    /**
     * Copy all p in st to stpre
     * Dimensions of @c st and that of @c pstruct must macth,and both with valid 
     * p
     */
    int Backupdata3d(const data3d< DataType > &mstru);

    /**
     * @brief Save p of  data3d data skipping p rows and p columns 
     * 
     */
    void CaptData(const unsigned num, unsigned leap = 0);

    void operator=(data3d< DataType > const &other);
    void InitStructData(DataType initVal = 0);
    void SaveData(unsigned leap, unsigned step);
    void SaveData(unsigned k, unsigned leap, unsigned step);
    void SaveXPlain(unsigned i, unsigned leap, unsigned step);
    void SaveYPlain(unsigned j, unsigned leap, unsigned step);
    void SaveZPlain(unsigned k, unsigned leap, unsigned step);
    void SaveData(unsigned k, unsigned leap, unsigned step, int type);
    void save(unsigned k, unsigned leap, unsigned step, int type);

    /**
     * @brief Create a data3d with the same size;
     * @param stru the source data3d to be copied.
     * @return
     */
    int CreateStruct(const data3d< DataType > &stru);
    /**
     * @brief Create a data3d with the same size as @c stru and initial all var to @c initVal;
     * @param stru the source data3d to be copied.
     */
    int CreateStruct(const data3d< DataType > &stru, DataType initVal);

    //set name

    void setName(const std::string &sn) {
        name = sn;
    }

    string getName() {
        return name;
    }
public:
    void ClearSim();
    void PlotArrays();
    void InitPlot();
public:
    static int InitMatlabEngine();
    static int CloseEngine();
#if DEBUG
    void nanOperator(unsigned i, unsigned j, unsigned k);
#endif

};

#include "datastruct.cpp"

#endif // DATASTRUCT_H

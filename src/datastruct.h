/*
<one line to give the program's name and a brief idea of what it does.>
Copyright (C) 2011  skiloop <skiloop@126.com>

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

#endif // end MATLAB_SIMULATION
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;
//#include "microdef.h"
#define MAX_ARRAY_SIZE 300000
#include "common.h"

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
    void CreateStruct(unsigned num,T val) {
        createArray(num,val);
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

template<class DataType> unsigned int data3d<DataType>::cnt = 0;
template<class DataType> string data3d<DataType>::tail = ".dat";


#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>

using namespace std;

#ifdef MATLAB_SIMULATION
template<class DataType> Engine* data3d<DataType>::ep = NULL;
template<class DataType> bool data3d<DataType>::isMatlabEngineStarted = false;
#endif

template<class DataType> data3d<DataType>::data3d(const data3d< DataType >& obj) : p(NULL) {
    CreateStruct(obj);
    try {
        unsigned i, j;
        for (i = 0; i < obj.nx; i++) {
            for (j = 0; j < obj.ny; j++) {
                memcpy(p[i][j], obj.p[i][j], obj.ny * sizeof (DataType));
            }
        }
    } catch (exception & e) {
        cerr << e.what() << endl;
    }
}

template<class DataType>
data3d<DataType>::~data3d() {

    if (p == NULL)
        return;
    unsigned i, j;
    for (i = 0; i < nx; i++) {
        if (p[i] != NULL) {
            for (j = 0; j < ny; j++) {
                if (p[i][j] != NULL)
                    delete [] p[i][j];
            }
            delete [] p[i];
        }
    }
    //delete []p;
    delete []p;
}

template<class DataType>
int data3d<DataType>::CreateStruct(unsigned nnx, unsigned nny, unsigned nnz) {
#if(DEBUG>=6)
    cout << "(" << nnx << ',' << nny << ',' << nnz << ')' << endl;
#endif
    if (p != NULL) {
        cerr << "pointer not NULL" << endl;
        return -1;
    }
    unsigned i, j;
    if (nnx == 0 || nny == 0 || nnz == 0 || nnx > MAX_ARRAY_SIZE || nny > MAX_ARRAY_SIZE || nnz > MAX_ARRAY_SIZE) {
        cerr << "Invalid size!" << endl;
        return -1;
    }
    try {

        p = new DataType**[nnx];
        if (p == NULL) {
            cout << "Failed to create space for data3d!" << endl;
            return -1;
        }
        for (i = 0; i < nnx; i++) {
            p[i] = new DataType*[nny];
        }
        for (i = 0; i < nnx; i++) {
            for (j = 0; j < nny; j++) {
                p[i][j] = new DataType[nnz];
            }
        }
        nx = nnx;
        ny = nny;
        nz = nnz;
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
int data3d<DataType>::CreateStruct(unsigned nnx, unsigned nny, unsigned nnz, DataType initVal) {
    if (CreateStruct(nnx, nny, nnz) < 0)
        return -1;
    return ResetStructData(initVal);
}

template<class DataType>
int data3d<DataType>::ResetStructData(DataType Val) {
    unsigned i, j, k;
    if (!CheckStruct())
        return -1;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            for (k = 0; k < nz; k++) {
                p[i][j][k] = Val;
            }
        }
    }
    return 0;
}

template<class DataType>
void data3d<DataType>::PrintData() {
    unsigned i, j, k;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            for (k = 0; k < nz; k++) {
                cout << p[i][j][k] << '\t';
            }
            cout << endl;
        }
        cout << endl;
    }
}

template<class DataType>	
bool data3d<DataType>::CheckStruct() {
    if (ny <= 0 || nx <= 0 || nz <= 0 || p == NULL) {
        return false;
    }
    return true;
}

template<class DataType>
void data3d<DataType>::CaptData(const unsigned num, unsigned leap) {
    unsigned i, j, k;

    stringstream ss;
    ss << name << num << tail;
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
            for (k = 0; k < ny; k += leap) {
                ofile << p[i][j][k] << '\t';
            }
            ofile << endl;
        }
        ofile << endl;
    }
    ofile.close();
}

template<class DataType>
void data3d<DataType>::operator =(data3d< DataType > const &other) {
    unsigned i, j;
    if (&other == this || p == NULL || other.p == NULL ||
            other.nx != nx || other.ny != ny || other.nz != nz)
        return;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            memcpy(p[i][j], other.p[i][j], nz * sizeof (DataType));
        }
    }
}

template<class DataType>
void data3d<DataType>::InitStructData(DataType initVal) {
    unsigned i, j, k;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            for (k = 0; k < nz; k++) {
                p[i][j][k] = initVal;
            }
        }
    }
}

template<class DataType>
void data3d<DataType>::SaveData(unsigned leap, unsigned step) {
    stringstream ss;
    ss << name << "_" << step << tail;
    string fname = ss.str();
    if (leap >= nx || leap >= ny || leap >= nz) {
        cerr << "Invalid leap for saving p!" << endl;
        return;
    }
    ofstream out(fname.c_str(), ios_base::binary);
    if (!out.is_open()) {
        cerr << "File " << fname << "cannot be opened!" << endl;
        return;
    }
    if (leap <= 0)leap = 1;
    unsigned i, j, k;
    for (i = 0; i < nx; i += leap) {
        for (j = 0; j < ny; j += leap) {
            for (k = 0; k < ny; k += leap) {
                out << p[i][j][k] << '\t';
            }
            out << endl;
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
void data3d<DataType>::save(unsigned k, unsigned leap, unsigned step, int type) {
    SaveData(k, leap, step, type);
}

/**
 * 
 * @param k
 * @param leap
 */
template<class DataType>
void data3d<DataType>::SaveData(unsigned k, unsigned leap, unsigned step, int type) {
    switch (type) {
        case 1:
            SaveXPlain(k, leap, step);
            break;
        case 2:
            SaveYPlain(k, leap, step);
            break;
        case 3:
            SaveZPlain(k, leap, step);
            break;
        default:
            SaveZPlain(k, leap, step);
    }
}

/**
 * 
 * @param k
 * @param leap
 */
template<class DataType>
void data3d<DataType>::SaveData(unsigned k, unsigned leap, unsigned step) {
    SaveZPlain(k, leap, step);
}

/**
 * 
 * @param k
 * @param leap
 * @param step
 */
template<class DataType>
void data3d<DataType>::SaveZPlain(unsigned k, unsigned leap, unsigned step) {
    stringstream ss;
    ss << name << "_z_" << step << tail;
    if (leap >= nx || leap >= ny || leap >= nz) {
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
    unsigned i, j;
    if (k >= nz)k = nz / 2;
    for (i = 0; i < nx; i += leap) {
        for (j = 0; j < ny; j += leap) {
            out << p[i][j][k] << '\t';
        }
        out << endl;
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
void data3d<DataType>::SaveYPlain(unsigned k, unsigned leap, unsigned step) {
    stringstream ss;
    ss << name << "_y_" << step << tail;
    if (leap >= nx || leap >= ny || leap >= nz) {
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
    unsigned i, j;
    if (k >= nz)k = nz / 2;
    for (i = 0; i < nx; i += leap) {
        for (j = 0; j < nz; j += leap) {
            out << p[i][k][j] << '\t';
        }
        out << endl;
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
void data3d<DataType>::SaveXPlain(unsigned k, unsigned leap, unsigned step) {
    stringstream ss;
    ss << name << "_x_" << step << tail;
    if (leap >= nx || leap >= ny || leap >= nz) {
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
    unsigned i, j;
    if (k >= nx)k = nx / 2;
    for (j = 0; j < ny; j += leap) {
        for (i = 0; i < nz; i += leap) {
            out << p[k][j][i] << '\t';
        }
        out << endl;
    }
    out.close();
}

template<class DataType>
int data3d<DataType>::InitMatlabEngine() {

#ifdef MATLAB_SIMULATION
    if (isMatlabEngineStarted) {
        return 0;
    }
    if (ep != NULL) {
        return -2;
    }
    if ((ep = engOpen(NULL)) == NULL) {
        cerr << "Can't start matlab engine!" << endl;
        return -1;
    }
    isMatlabEngineStarted = true;
#endif
    return 0;

}

template<class DataType>
int data3d<DataType>::CloseEngine() {
#ifdef MATLAB_SIMULATION
    if (isMatlabEngineStarted) {
        engEvalString(ep, "close all;clear;");
        engClose(ep);
    }
#endif
    return 0;
}

template<class DataType>
int data3d<DataType>::CreateStruct(const data3d< DataType > &stru) {
    return CreateStruct(stru.nx, stru.ny, stru.nz);
}

template<class DataType>
int data3d<DataType>::CreateStruct(const data3d< DataType > &stru, DataType initVal) {
    return CreateStruct(stru.nx, stru.ny, stru.nz, initVal);
}

template<class DataType>
void data3d<DataType>::ClearSim() {
#ifdef MATLAB_SIMULATION
    if (!isMatlabEngineStarted)return;
    mxDestroyArray(MyArray);
    mxDestroyArray(num);
#endif
}

template<class DataType>
void data3d<DataType>::PlotArrays() {
#ifdef MATLAB_SIMULATION
    if (!isMatlabEngineStarted)return;

    DataType *pData = (DataType*) malloc(nx * ny * sizeof (DataType));
    for (unsigned i = 0; i < nx; i++)
        for (unsigned j = 0; j < ny; j++)
            pData[i * ny + j] = p[i][j][nz / 2];
    engPutVariable(ep, "ind", num);
    engEvalString(ep, "ind=int32(ind);");
    memcpy(mxGetPr(MyArray), pData, nx * ny * sizeof (DataType));
    engPutVariable(ep, "array", MyArray);
    engEvalString(ep, "obj(ind).array=array;clear array;");
    engEvalString(ep, "set(obj(ind).img,'CData',obj(ind).array);drawnow;");
    free(pData);
#endif
}

template<class DataType>
void data3d<DataType>::InitPlot() {
#ifdef MATLAB_SIMULATION
    if (!isMatlabEngineStarted)return;
#endif
    cnt++;
    string filename = name + tail;
    Number = cnt;
#ifdef MATLAB_SIMULATION

    mxArray *mxStr = mxCreateString(filename.c_str());
    DataType *pData = (DataType*) malloc(nx * ny * sizeof (DataType));
    if (pData == NULL)return;
    MyArray = mxCreateDoubleMatrix(ny, nx, mxREAL);
    num = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    memcpy(mxGetPr(num), &Number, sizeof (unsigned));
    engPutVariable(ep, "ind", num);
    engEvalString(ep, "ind=int32(ind);");
    engEvalString(ep, "obj(ind).fig=figure('NumberTitle','OFF');");

    engPutVariable(ep, "name", mxStr);
    engEvalString(ep, "obj(ind).name=name;");

    for (unsigned i = 0; i < nx; i++)
        for (unsigned j = 0; j < ny; j++)
            pData[i * ny + j] = p[i][j][nz / 2];
    memcpy(mxGetPr(MyArray), pData, nx * ny * sizeof (DataType));
    engPutVariable(ep, "array", MyArray);

    engEvalString(ep, "obj(ind).array=array;clear array;");
    engEvalString(ep, "obj(ind).img=imagesc(obj(ind).array);obj(ind).ax=gca;title(obj(ind).ax,obj(ind).name);drawnow;");
    engEvalString(ep, "set(gca,'YDir','Normal');colorbar;");
    free(pData);
    mxDestroyArray(mxStr);
#endif
}

#if DEBUG
template<class DataType>
void data3d<DataType>::nanOperator(unsigned i,unsigned j,unsigned k){
    if(isnan(p[i][j][k])){
        cout<<"nan var found for "<< getName()<<" at:("<<i<<","<<j<<","<<k<<")"<<endl;
    }
}
#endif

//#include "datastruct.cpp"

#endif // DATASTRUCT_H

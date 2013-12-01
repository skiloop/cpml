#ifndef DATA3D_H
#define DATA3D_H


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
#include "Point.h"

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
    std::string name; // name for this to save file

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

    /**
     * default constructor set
     *  @c nx = 0
     *  @c ny = 0
     *  @c nz = 0
     *  @c p=NULL
     */
    data3d() : nx(0), ny(0), nz(0), p(NULL) {
    };

    /**
     * copy constructor
     * @param obj
     */
    data3d(const data3d< DataType > &obj);

    /**
     * deconstructor
     */
    ~data3d();

    /**
     * print p in struct data3d
     */
    void printArray();

    /**
     * free space created for data3d @c mst
     */
    void FreeStructData();

    /**
     * Set Data to val
     */
    int resetArray(DataType val = 0);

    /**
     * check p of data3d @c mst is valid
     * if p is not NULL and none of its subpointers,then
     * return true,otherwise false
     */
    bool checkArray();

    /**
     * Create Space for struct data3d and initialize its @c nx and @c ny
     */
    int create3DArray(unsigned nnx, unsigned nny, unsigned nnz);

    /**
     * Create Space for struct data3d and initialize its @c nx and @c ny
     */
    int create3DArray(unsigned nnx, unsigned nny, unsigned nnz, DataType initVal);

    /**
     * Copy all p in st to stpre
     * Dimensions of @c st and that of @c pstruct must macth,and both with valid
     * p
     */
    int backup3DArray(const data3d< DataType > &mstru);

    /**
     * @brief Save p of  data3d data skipping p rows and p columns
     *
     */
    void saveArrayData(const unsigned num, unsigned leap = 0);

    /**
     *
     * @param other
     */
    void operator=(data3d< DataType > const &other);

    /**
     *  return this.p[index.x][index.y][index.z]
     *
     * @param index
     * @return
     */
    DataType operator[](const Point index) const;

    /**
     * initial array to @c initVal
     * @param initVal
     */
    void initArray(DataType initVal = 0);

    /**
     * save array
     * @param leap
     * @param step
     */
    void saveData(unsigned leap, unsigned step);

    /**
     *
     * @param k
     * @param leap
     * @param step
     */
    void saveData(unsigned k, unsigned leap, unsigned step);

    /**
     *
     * @param i
     * @param leap
     * @param step
     */
    void saveXPlain(unsigned i, unsigned leap, unsigned step);

    /**
     *
     * @param j
     * @param leap
     * @param step
     */
    void saveYPlain(unsigned j, unsigned leap, unsigned step);

    /**
     *
     * @param k
     * @param leap
     * @param step
     */
    void saveZPlain(unsigned k, unsigned leap, unsigned step);

    /**
     *
     * @param k
     * @param leap
     * @param step
     * @param type
     */
    void savePlain(unsigned k, unsigned leap, unsigned step, int type);

    /**
     *
     * @param k
     * @param leap
     * @param step
     * @param type
     */
    void save(unsigned k, unsigned leap, unsigned step, int type);

    /**
     * @brief Create a data3d with the same size;
     * @param stru the source data3d to be copied.
     * @return
     */
    int create3DArray(const data3d< DataType > &stru);

    /**
     * @brief Create a data3d with the same size as @c stru and initial all var to @c initVal;
     * @param stru the source data3d to be copied.
     */
    int create3DArray(const data3d< DataType > &stru, DataType initVal);

    /**
     * set @name to @sn
     * @param sn
     */
    void setName(const std::string &sn) {
        name = sn;
    }

    /**
     * get name
     * @return @c name
     */
    string getName() {
        return name;
    }

public:

    void clearMatlabEngineArray();
    void plotArrays();
    void preparePlotting();
public:
    static int initMatlabEngine();
    static int closeMatlabEngine();
    bool isNaN(unsigned i, unsigned j, unsigned k);
    bool isInf(unsigned i, unsigned j, unsigned k);
    bool isValid(unsigned i, unsigned j, unsigned k);

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
    create3DArray(obj);
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
int data3d<DataType>::create3DArray(unsigned nnx, unsigned nny, unsigned nnz) {
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
int data3d<DataType>::create3DArray(unsigned nnx, unsigned nny, unsigned nnz, DataType initVal) {
    if (create3DArray(nnx, nny, nnz) < 0)
        return -1;
    return resetArray(initVal);
}

template<class DataType>
int data3d<DataType>::resetArray(DataType Val) {
    unsigned i, j, k;
    if (!checkArray())
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
void data3d<DataType>::printArray() {
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
bool data3d<DataType>::checkArray() {
    if (ny <= 0 || nx <= 0 || nz <= 0 || p == NULL) {
        return false;
    }
    return true;
}

template<class DataType>
void data3d<DataType>::saveArrayData(const unsigned num, unsigned leap) {
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
DataType data3d<DataType>::operator[](const Point index) const {
    return p[index.x][index.y][index.z];
}

template<class DataType>
void data3d<DataType>::initArray(DataType initVal) {
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
void data3d<DataType>::saveData(unsigned leap, unsigned step) {
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
    savePlain(k, leap, step, type);
}

/**
 *
 * @param k
 * @param leap
 */
template<class DataType>
void data3d<DataType>::savePlain(unsigned k, unsigned leap, unsigned step, int type = 3) {
    switch (type) {
        case 1:
            saveXPlain(k, leap, step);
            break;
        case 2:
            saveYPlain(k, leap, step);
            break;
        case 3:
            saveZPlain(k, leap, step);
            break;
        default:
            saveZPlain(k, leap, step);
    }
}

/**
 *
 * @param k
 * @param leap
 * @param step
 */
template<class DataType>
void data3d<DataType>::saveZPlain(unsigned k, unsigned leap, unsigned step) {
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
void data3d<DataType>::saveYPlain(unsigned k, unsigned leap, unsigned step) {
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
void data3d<DataType>::saveXPlain(unsigned k, unsigned leap, unsigned step) {
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
int data3d<DataType>::initMatlabEngine() {

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
int data3d<DataType>::closeMatlabEngine() {
#ifdef MATLAB_SIMULATION
    if (isMatlabEngineStarted) {
        engEvalString(ep, "close all;clear;");
        engClose(ep);
    }
#endif
    return 0;
}

template<class DataType>
int data3d<DataType>::create3DArray(const data3d< DataType > &stru) {
    return create3DArray(stru.nx, stru.ny, stru.nz);
}

template<class DataType>
int data3d<DataType>::create3DArray(const data3d< DataType > &stru, DataType initVal) {
    return create3DArray(stru.nx, stru.ny, stru.nz, initVal);
}

template<class DataType>
void data3d<DataType>::clearMatlabEngineArray() {
#ifdef MATLAB_SIMULATION
    if (!isMatlabEngineStarted)return;
    mxDestroyArray(MyArray);
    mxDestroyArray(num);
#endif
}

template<class DataType>
void data3d<DataType>::plotArrays() {
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
void data3d<DataType>::preparePlotting() {
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

/**
 * check if is nan at point (i,j,k)
 *
 */
template<class DataType>
bool data3d<DataType>::isNaN(unsigned i, unsigned j, unsigned k) {
    if (isnan(p[i][j][k])) {
        cout << "nan var found for " << getName() << " at:(" << i << "," << j << "," << k << ")" << endl;
        return true;
    }
    return false;
}

/**
 * check if is inf at point (i,j,k)
 *
 */
template<class DataType>
bool data3d<DataType>::isInf(unsigned i, unsigned j, unsigned k) {
    if (isinf(p[i][j][k])) {
        cout << "inf var found for " << getName() << " at:(" << i << "," << j << "," << k << ")" << endl;
        return true;
    }
    return false;
}

/**
 * check if is valid var at point (i,j,k)
 */
template<class DataType>
bool data3d<DataType>::isValid(unsigned i, unsigned j, unsigned k) {
    return isNaN(i, j, k) | isInf(i, j, k);
}



#endif // DATA3D_H

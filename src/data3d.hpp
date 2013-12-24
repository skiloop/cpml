
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>

#include "data3d.h"

using namespace std;

#ifdef MATLAB_SIMULATION
template<class DataType> Engine* data3d<DataType>::ep = NULL;
template<class DataType> bool data3d<DataType>::mIsMatlabEngineStarted = false;
#endif

template<class DataType> const string data3d<DataType>::OUTPUT_FILE_NAME_TAIL = ".dat";

#ifdef MATLAB_SIMULATION
template<class DataType> unsigned int data3d<DataType>::mMatlabFigureCount = 0;
#endif

template<class DataType> data3d<DataType>::data3d(const data3d< DataType >& obj) : p(NULL) {
#ifdef MATLAB_SIMULATION
    mMatlabFigureIndex = -1;
#endif
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
    freeArray();
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
    ss << mName << "_" << step << OUTPUT_FILE_NAME_TAIL;
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

template<class DataType>
void data3d<DataType>::save(int leap) {
    stringstream ss;
    ss << mName << OUTPUT_FILE_NAME_TAIL;
    string fname = ss.str();

    ofstream out(fname.c_str(), ios_base::binary);
    if (!out.is_open()) {
        cerr << "File " << fname << "cannot be opened!" << endl;
        return;
    }
    if (leap <= 0)leap = 1;
    unsigned i, j, k;
    for (i = 0; i < nx; i += leap) {
        for (j = 0; j < ny; j += leap) {
            for (k = 0; k < nz; k += leap) {
                out << p[i][j][k] << '\t';
            }
            out << endl;
        }
        // out << endl;
    }
    out.close();
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
    ss << mName << "_z_" << step << OUTPUT_FILE_NAME_TAIL;
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
    ss << mName << "_y_" << step << OUTPUT_FILE_NAME_TAIL;
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
    ss << mName << "_x_" << step << OUTPUT_FILE_NAME_TAIL;
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
    if (!mIsMatlabEngineStarted) {
        if (ep != NULL) {
            return -2;
        }
        if ((ep = engOpen(NULL)) == NULL) {
            cerr << "Can't start matlab engine!" << endl;
            return -1;
        }
        setMatlabEngineStarted(true);
    }
#endif
    return 0;

}

template<class DataType>
int data3d<DataType>::closeMatlabEngine() {
#ifdef MATLAB_SIMULATION
    if (mIsMatlabEngineStarted) {
        engEvalString(ep, "close all;clear;");
        engClose(ep);
        setMatlabEngineStarted(false);
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
    if (mIsMatlabEngineStarted && mMatlabFigureIndex > 0) {
        mxDestroyArray(mMatlabMXArray);
        mxDestroyArray(num);
    }
#endif
}

template<class DataType>
void data3d<DataType>::plotArrays() {
#ifdef MATLAB_SIMULATION
    if (mIsMatlabEngineStarted) {
        DataType *pData = (DataType*) malloc(nx * ny * sizeof (DataType));
        for (unsigned i = 0; i < nx; i++) {
            for (unsigned j = 0; j < ny; j++) {
                pData[i * ny + j] = p[i][j][nz / 2];
            }
        }
        engPutVariable(ep, "ind", num);
        engEvalString(ep, "ind=int32(ind);");
        memcpy(mxGetPr(mMatlabMXArray), pData, nx * ny * sizeof (DataType));
        engPutVariable(ep, "array", mMatlabMXArray);
        engEvalString(ep, "obj(ind).array=array;clear array;");
        engEvalString(ep, "set(obj(ind).img,'CData',obj(ind).array);drawnow;");
        free(pData);
    }
#endif
}

template<class DataType>
void data3d<DataType>::preparePlotting() {
#ifdef MATLAB_SIMULATION
    if (mIsMatlabEngineStarted) {
        string filename = mName + OUTPUT_FILE_NAME_TAIL;
        mMatlabFigureIndex = ++mMatlabFigureCount;

        mxArray *mxStr = mxCreateString(filename.c_str());
        DataType *pData = (DataType*) malloc(nx * ny * sizeof (DataType));
        if (pData == NULL)return;
        mMatlabMXArray = mxCreateDoubleMatrix(ny, nx, mxREAL);
        num = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        memcpy(mxGetPr(num), &mMatlabFigureIndex, sizeof (unsigned));
        engPutVariable(ep, "ind", num);
        engEvalString(ep, "ind=int32(ind);");
        engEvalString(ep, "obj(ind).fig=figure('NumberTitle','OFF');");

        engPutVariable(ep, "name", mxStr);
        engEvalString(ep, "obj(ind).name=name;");

        for (unsigned i = 0; i < nx; i++) {
            for (unsigned j = 0; j < ny; j++) {
                pData[i * ny + j] = p[i][j][nz / 2];
            }
        }
        memcpy(mxGetPr(mMatlabMXArray), pData, nx * ny * sizeof (DataType));
        engPutVariable(ep, "array", mMatlabMXArray);

        engEvalString(ep, "obj(ind).array=array;clear array;");
        engEvalString(ep, "obj(ind).img=imagesc(obj(ind).array);obj(ind).ax=gca;title(obj(ind).ax,obj(ind).name);drawnow;");
        engEvalString(ep, "set(gca,'YDir','Normal');colorbar;");
        free(pData);
        mxDestroyArray(mxStr);
    }
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

template<class DataType>
void data3d<DataType>::whenLargerThan(unsigned i, unsigned j, unsigned k, MyDataF limit, void (*fun)()) {
    if (p[i][j][k] > limit) {
        //cout << "larger var found for " << getName() << " at:(" << i << "," << j << "," << k << ")" << endl;
        if (NULL != fun) {
            (*fun)();
        }
    }
}

template<class DataType>
void data3d<DataType>::freeArray() {
    if (p != NULL) {
        for (unsigned i = 0; i < nx; i++) {
            if (p[i] != NULL) {
                for (unsigned j = 0; j < ny; j++) {
                    if (p[i][j] != NULL) delete [] p[i][j];
                }
                delete [] p[i];
            }
        }
        delete []p;
    }
}
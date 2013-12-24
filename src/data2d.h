/* 
 * File:   data2d.h
 * Author: skiloop
 *
 * Created on 2013年12月17日, 下午5:08
 */

#ifndef DATA2D_H
#define	DATA2D_H


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;

//#include "microdef.h"
#define MAX_ARRAY_SIZE 300000
#include "common.h"
#include "Point.h"

/**
 * nx:point count in x direction
 * ny:point count in x direction
 * p:pointer to store data;
 */
template<class DataType>
class data2d {
public:
    unsigned int nx;
    unsigned int ny;
    DataType** p;

public:
    static const std::string OUTPUT_FILE_NAME_TAIL;
private:
    std::string mName; // name for this to save file

private:

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
    data2d(unsigned int cx, unsigned int cy)
    : nx(cx), ny(cy),  p(NULL)
    {
        unsigned i, j;
        if (cx == 0 || cy == 0 ) {
            return;
        }
        try {

            p = new DataType*[cx];
            for (i = 0; i < cx; i++) {
                p[i] = new DataType[cy];
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
     *  @c p=NULL
     */
    data2d() : nx(0), ny(0),  p(NULL)

    {
    };

    /**
     * copy constructor
     * @param obj
     */
    data2d(const data2d< DataType > &obj);

    /**
     * deconstructor
     */
    ~data2d();

    /**
     * print p in struct data2d
     */
    void printArray();

    /**
     * free space created for data2d @c mst
     */
    void freeArray();

    /**
     * Set Data to val
     */
    int resetArray(DataType val = 0);

    /**
     * check p of data2d @c mst is valid
     * if p is not NULL and none of its subpointers,then
     * return true,otherwise false
     */
    bool checkArray();

    /**
     * Create Space for struct data2d and initialize its @c nx and @c ny
     */
    int create2DArray(unsigned nnx, unsigned nny);

    /**
     * Create Space for struct data2d and initialize its @c nx and @c ny
     */
    int create2DArray(unsigned nnx, unsigned nny, DataType initVal);

    /**
     * Copy all p in st to stpre
     * Dimensions of @c st and that of @c pstruct must macth,and both with valid
     * p
     */
    int backup2DArray(const data2d< DataType > &mstru);

    /**
     * @brief Save p of  data2d data skipping p rows and p columns
     *
     */
    void saveArrayData(const unsigned num, unsigned leap = 0);

    /**
     *
     * @param other
     */
    void operator=(data2d< DataType > const &other);

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
     * @param type
     */
    void savePlain(unsigned k, unsigned leap, unsigned step, int type);

    /**
     *  Save data at plain s=@c k where s=x,y which define by @c type
     * @param k
     * @param leap
     * @param step
     * @param type
     */
    void save(unsigned k, unsigned leap, unsigned step, int type);

    /**
     * save every @c leap cells data to file 
     * @param leap
     */
    void save(int leap = 1);

    /**
     * @brief Create a data2d with the same size;
     * @param stru the source data2d to be copied.
     * @return
     */
    int create2DArray(const data2d< DataType > &stru);

    /**
     * @brief Create a data2d with the same size as @c stru and initial all var to @c initVal;
     * @param stru the source data2d to be copied.
     */
    int create2DArray(const data2d< DataType > &stru, DataType initVal);

    /**
     * set @name to @sn
     * @param sn
     */
    void setName(const std::string &sn) {
        mName = sn;
    }

    /**
     * get name
     * @return @c name
     */
    string getName() {
        return mName;
    }

public:

    void clearMatlabEngineArray();
    void plotArrays();
    void preparePlotting();
public:
    static int initMatlabEngine();
    static int closeMatlabEngine();
    bool isNaN(unsigned i, unsigned j);
    bool isInf(unsigned i, unsigned j);
    bool isValid(unsigned i, unsigned j);
    /**
     * when value at (i,j,k) is larger than limit do something define by fun
     * @param i
     * @param j
     * @param k
     * @param limit
     * @param fun
     */
    void whenLargerThan(unsigned i, unsigned j, MyDataF limit, void(*fun)());

};

#include "data2d.hpp"

#endif	/* DATA2D_H */


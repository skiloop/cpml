
#include<iostream>
#ifdef _OPENMP
#include <cstdlib>
#include <omp.h>
int thread_count = 2; // thread count for openmp
#endif

#include "fdtd.h"

MyDataF epsR;
MyDataF dx, dy, dz;
MyDataF tw;
MyDataF T; // ns
MyDataF frequency=110E9;
MyDataF amp=100.0;

int main(int argc, char*argv[]) {

    int cellSize=100;
    int pmlWidth=10;
    epsR = 1.0;

    tw = T = 1 / frequency;
    dx = C * T / cellSize;
    dy = C * T / cellSize;
    dz = C * T / cellSize;

    unsigned xlen, ylen, zlen, tlen;
    

    //    MyDataF dt = 0.99 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1 / (dz * dz)));
    MyDataF dt = dx / 2 / C;
    xlen = 100;
    ylen = 100;
    zlen = 50;
    tlen=2000;

    cout << "xlen=" << xlen << endl;
    cout << "ylen=" << ylen << endl;
    cout << "zlen=" << zlen << endl;
    cout << "tlen=" << tlen << endl;
    cout << "dx=" << dx << endl;
    cout << "dt=" << dt << endl;

    fdtd cpmltest(tlen, xlen, ylen, zlen, tw, dx, dy, dz, amp, 10, 12, 4, 1, pmlWidth);

    //hpw.initialize();
    cpmltest.StartUp();
    return 0;
}

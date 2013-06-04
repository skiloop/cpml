#pragma once

#include "../src/cpml.h"

#ifndef WITH_DENSITY
//#define WITH_DENSITY
#endif

class fdtd {
public:

    fdtd(unsigned _totalTimeSteps = 500, unsigned _imax = 40, unsigned _jmax = 40, unsigned _kmax = 26,
            MyDataF _tw = 53.0e-12, MyDataF _dx = 1e-3, MyDataF _dy = 1e-3, MyDataF _dz = 1e-3,
            MyDataF _amp = 1000, unsigned _savemodulus = 10, unsigned _ksource = 12,
            unsigned _m = 4, unsigned _ma = 1, unsigned pmlw = 10);

    ~fdtd(void);

    //Function prototype definitions
    void initialize(); //Memory initialization
    void setUp(); //Coefficients, parameters etc will get computed
    void compute(); //E & H Field update equation
    void StartUp();

private:
    void yeeCube(unsigned, unsigned, unsigned, unsigned); //Sets material properties to a cell
    void writeField(unsigned); //Writes output
    void putvars();
private:
    //  Specify Number of Time Steps and Grid Size Parameters
    unsigned totalTimeSteps; // total number of time steps

    // grid size corresponding to the number of Ez field components
    unsigned Imax;
    unsigned Jmax;
    unsigned Kmax;
    //  Specify the Impulsive Source (Differentiated Gaussian) parameters
    MyDataF tw; //pulse width
    MyDataF dt, dx, dy, dz;

    MyDataF amp; // Amplitude
    // Specify the Time Step at which the data has to be saved for Visualization
    unsigned save_modulus;

    // Output recording 
    unsigned ksource;
    // Specify the CPML Order and Other Parameters:
    unsigned m;
    unsigned ma;
    unsigned pmlWidth;

    MyDataF t0; //delay
    MyDataF source; //Differentiated Gaussian source    

    // source position
    unsigned isp, jsp, ksp;

    // H & E Field components
    data3d<MyDataF> Hx;
    data3d<MyDataF> Hy;
    data3d<MyDataF> Hz;
    data3d<MyDataF> Ex;
    data3d<MyDataF> Ey;
    data3d<MyDataF> Ez;

    data3d<MyDataF> Cexe, Ceye, Ceze;
    data3d<MyDataF> Chxh, Chyh, Chzh;
    data3d<MyDataF> Cexhy, Ceyhz, Cezhx, Cexhz, Ceyhx, Cezhy;
    data3d<MyDataF> Chxey, Chyez, Chzex, Chxez, Chyex, Chzey;
    void initCoeficients();

    void updateHx();
    void updateHy();
    void updateHz();
    void updateMagneitcFields();
    void updateEx();
    void updateEy();
    void updateEz();
    void updateElectricFields();
    void updateSource(unsigned n);
    cpml<MyDataF> pml;

};



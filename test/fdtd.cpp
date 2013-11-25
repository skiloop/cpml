

#include <math.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>     
#include <assert.h>
#include <fstream>
#include <sstream>
#include <exception>

#ifdef _OPENMP
#include <omp.h>
extern int thread_count;
#endif

#include "common.h"
#include "../src/cpml.h"
#include "fdtd.h"

extern MyDataF epsR;
//extern MyDataF dt, dx, dy, dz;
extern MyDataF T;

using namespace std;

fdtd::fdtd(unsigned _totalTimeSteps, unsigned _imax, unsigned _jmax, unsigned _kmax,
        MyDataF _tw, MyDataF _dx, MyDataF _dy, MyDataF _dz,
        MyDataF _amp, unsigned _savemodulus, unsigned _ksource,
        unsigned _m, unsigned _ma, unsigned pmlw)
: totalTimeSteps(_totalTimeSteps), Imax(_imax), Jmax(_jmax), Kmax(_kmax)
, tw(_tw), dx(_dx), dy(_dy), dz(_dz)
, amp(_amp), save_modulus(_savemodulus), ksource(_ksource)
, m(_m), ma(_ma), pmlWidth(pmlw) {
}

fdtd::~fdtd(void) {
}

void fdtd::initialize() {

    Ez.create3DArray(Imax + 1, Jmax + 1, Kmax, 0);
    Ey.create3DArray(Imax + 1, Jmax, Kmax + 1, 0);
    Ex.create3DArray(Imax, Jmax + 1, Kmax + 1, 0);

    Hx.create3DArray(Imax + 1, Jmax, Kmax, 0);
    Hy.create3DArray(Imax, Jmax + 1, Kmax, 0);
    Hz.create3DArray(Imax, Jmax, Kmax + 1, 0);

    //coefficients
    Cexe.create3DArray(Ex, 0.0);
    Ceye.create3DArray(Ey, 0.0);
    Ceze.create3DArray(Ez, 0.0);
    Chxh.create3DArray(Hx, 0.0);
    Chyh.create3DArray(Hy, 0.0);
    Chzh.create3DArray(Hz, 0.0);
    Cexhy.create3DArray(Ex, 0.0);
    Cexhz.create3DArray(Ex, 0.0);
    Chxey.create3DArray(Hx, 0.0);
    Chxez.create3DArray(Hx, 0.0);
    Ceyhx.create3DArray(Ey, 0.0);
    Ceyhz.create3DArray(Ey, 0.0);
    Chyex.create3DArray(Hy, 0.0);
    Chyez.create3DArray(Hy, 0.0);
    Cezhy.create3DArray(Ez, 0.0);
    Cezhx.create3DArray(Ez, 0.0);
    Chzey.create3DArray(Hz, 0.0);
    Chzex.create3DArray(Hz, 0.0);

    Ez.setName("Ez");
    Ex.setName("Ex");
    Ey.setName("Ey");
    Hz.setName("Hz");
    Hx.setName("Hx");
    Hy.setName("Hy");
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void fdtd::setUp() {
    //Time step
    //    dt = 0.99 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) +
    //            1.0 / (dz * dz)));
    dt = dx / 2 / C;

    //delay
    t0 = 3.0 * tw;

    // source position
    isp = Imax / 2;
    jsp = Jmax / 2;
    ksp = Kmax / 2;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  COMPUTING FIELD UPDATE EQUATION COEFFICIENTS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    initCoeficients();
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  PML parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MyDataF sigmaMax = 0.75;
    MyDataF kappaMax = 10;
    MyDataF alphaMax = 0.025;
    int pmlOrder = 4;

    //pml.initParmeters(dx, dy, dz, m, ma);
    pml.setCPMLRegion(pmlWidth);
    pml.createCPMLArrays(Imax, Jmax, Kmax);
    pml.initCoefficientArrays(pmlOrder, sigmaMax, kappaMax, alphaMax, epsR, dt, dx, dy, dz,
            Ceyhz, Cezhy, Chyez, Chzey,
            Cexhz, Cezhx, Chxez, Chzex,
            Ceyhx, Cexhy, Chyex, Chxey);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initial Coefficients
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    putvars();

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SET CPML PARAMETERS IN EACH DIRECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void fdtd::compute() {

    unsigned n;
    unsigned ic, jc, kc; //capture field
    //    ic = isp + 1;
    //    jc = jsp + 1;
    //    kc = ksp + 2;
    ic = isp; // + (Imax - isp) / 2;
    jc = jsp + (Jmax - jsp) / 2;
    kc = ksp;
    assert(ic < Imax && jc < Jmax && kc < Kmax);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  BEGIN TIME STEP
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cout << endl;
    cout << "Begin time-stepping..." << endl;
    for (n = 1; n <= totalTimeSteps; ++n) {

        cout << "Ez at time step " << n << " at (" << ic << ", " << jc << ", " << kc;
        cout << ") :  " << Ez.p[ic][jc][kc] << endl;

        updateMagneitcFields();
        pml.updateCPML_M_Fields(Hx, Hy, Hz, Ex, Ey, Ez);
        updateElectricFields();
        pml.updateCPML_E_Fields(Ex, Ey, Ez, Hx, Hy, Hz);
        //====================================
        // update Source
        //====================================
        updateSource(n);


        if ((n % save_modulus) == 0) {
            writeField(n);
            Ez.save(isp, 1, n, 1);
            Ez.save(jsp, 1, n, 2);
            Ez.save(ksp, 1, n, 3);
        }
    }

    //  END TIME STEP
    cout << "Done time-stepping..." << endl;

}

void fdtd::updateSource(unsigned n) {
    MyDataF source;

    source = amp * -2.0 * ((n * dt - t0) / tw / tw)
            * exp(-pow(((n * dt - t0) / tw), 2)); //Differentiated Gaussian pulse

    Ez.p[isp][jsp][ksp] = Ez.p[isp][jsp][ksp] + dt / eps_0 * source / dx / dy / dz;
    //cout<<"source="<<source<<"\t"<<amp<<"\t"<<n<<"\t"<<dt<<"\t"<<
    //        amp * -2.0 * ((n * dt - t0) / tw / tw) * exp(-pow(((n * dt - t0) / tw), 2))<<endl;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Saving Output Data to files
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void fdtd::writeField(unsigned iteration) {
    unsigned i, j;
    ofstream out;
    stringstream ss;
    // form file name
    ss << "E_Field_" << iteration << ".txt";
    string fileBaseName(ss.str());
    // open file
    out.open(fileBaseName.c_str());
    if (out.is_open()) {

        for (i = 0; i < Imax - 1; i++) {
            for (j = 0; j < Jmax - 1; j++) { // |E|
                out << sqrt(pow(Ex.p[i][j][ksource], 2) +
                        pow(Ey.p[i][j][ksource], 2) + pow(Ez.p[i][j][ksource], 2)) << '\t';
            }
            out << endl;
        }
        out.close();
    }
}

//start up

void fdtd::StartUp() {
    cout << "initializing(in Startup)..." << endl;
    initialize();
    //    cout << "initial pml (in Statup)" << endl;
    //    pml.Initial(Imax, Jmax, Kmax, 11);
    cout << "setUp (in Startup)" << endl;
    setUp();
    cout << "computing (in Startup)" << endl;
    compute();
    cout << "exit Startup" << endl;
}

void fdtd::putvars() {
    cout << "dx = " << dx << endl;
    cout << "dy = " << dy << endl;
    cout << "dz = " << dz << endl;
    cout << "(Imax,Jmax,Kmax) = (" << Imax << "," << Jmax << "," << Kmax << ")" << endl;
    // time step increment
    cout << "dt = " << dt << endl;

    //  Specify the Impulsive Source (Differentiated Gaussian) parameters
    cout << "tw = " << tw << endl; //pulse width
    cout << "t0 = " << t0 << endl; //delay
    cout << "source = " << source << endl; //Differentiated Gaussian source
    cout << "amp = " << amp << endl; // Amplitude

    //Specify the Time Step at which the data has to be saved for Visualization
    cout << "save_modulus = " << save_modulus << endl;

    //Output recording point
    cout << "ksource = " << ksource << endl;

    //  Specify the CPML Order and Other Parameters:
    cout << " m = " << m << endl;
    cout << " ma = " << ma << endl;

    cout << endl << "TIme step = " << dt << endl;
    cout << endl << "Number of steps = " << totalTimeSteps << endl;
    cout << endl << "Total Simulation time = " << totalTimeSteps * dt << " Seconds" << endl;
}

void fdtd::updateHx() {
    unsigned i, j, k;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k)//shared(Hx,Ez,Ey,pml,DA,DB,dy)
#endif
    for (k = 0; k < Kmax; ++k) {
        for (i = 0; i < Imax + 1; ++i) {
            for (j = 0; j < Jmax; ++j) {
                Hx.p[i][j][k] = Chxh.p[i][j][k] * Hx.p[i][j][k] + Chxez.p[i][j][k]*(Ez.p[i][j + 1][k] - Ez.p[i][j][k]) +
                        Chxey.p[i][j][k]*(Ey.p[i][j][k + 1] - Ey.p[i][j][k]);
            }
        }
    }
}

void fdtd::updateHy() {
    unsigned i, j, k;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k)//shared(Hy,Ez,Ex,pml,DA,DB,dx,dz)
#endif
    for (k = 0; k < Kmax; ++k) {
        for (i = 0; i < Imax; ++i) {
            for (j = 0; j < Jmax + 1; ++j) {
                Hy.p[i][j][k] = Chyh.p[i][j][k] * Hy.p[i][j][k] + Chyez.p[i][j][k]*(Ez.p[i + 1][j][k] - Ez.p[i][j][k]) +
                        Chyex.p[i][j][k]*(Ex.p[i][j][k + 1] - Ex.p[i][j][k]);
            }
        }
    }
}

void fdtd::updateHz() {
    unsigned i, j, k;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  UPDATE Hz
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k)//shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (k = 0; k < Kmax + 1; ++k) {
        for (i = 0; i < Imax; ++i) {
            for (j = 0; j < Jmax; ++j) {
                Hz.p[i][j][k] = Chzh.p[i][j][k] * Hz.p[i][j][k] + Chzey.p[i][j][k]
                        * (Ey.p[i + 1][j][k] - Ey.p[i][j][k]) +
                        (Ex.p[i][j + 1][k] - Ex.p[i][j][k]) * Chzex.p[i][j][k];
            }
        }
    }
}

void fdtd::updateEx() {
    unsigned i, j, k;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  UPDATE Ex
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k)
#endif
    for (k = 1; k < Kmax; ++k) {
        for (i = 0; i < Imax; ++i) {
            for (j = 1; j < Jmax; ++j) {
                Ex.p[i][j][k] = Cexe.p[i][j][k] * Ex.p[i][j][k] +
                        (Hz.p[i][j][k] - Hz.p[i][j - 1][k]) * Cexhz.p[i][j][k] +
                        (Hy.p[i][j][k] - Hy.p[i][j][k - 1]) * Cexhy.p[i][j][k];
            }
        }
    }
}

void fdtd::updateEy() {
    unsigned i, j, k;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  UPDATE Ey
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k)
#endif
    for (k = 1; k < Kmax; ++k) {
        for (i = 1; i < Imax; ++i) {
            for (j = 0; j < Jmax; ++j) {
                Ey.p[i][j][k] = Ceye.p[i][j][k] * Ey.p[i][j][k] +
                        (Hz.p[i][j][k] - Hz.p[i - 1][j][k]) * Ceyhz.p[i][j][k] +
                        (Hx.p[i][j][k] - Hx.p[i][j][k - 1]) * Ceyhx.p[i][j][k];
            }
        }
    }
}

void fdtd::updateEz() {
    unsigned i, j, k;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  UPDATE Ez
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) private(i,j,k)
#endif
    for (k = 0; k < Kmax; ++k) {
        for (i = 1; i < Imax; ++i) {
            for (j = 1; j < Jmax; ++j) {
                Ez.p[i][j][k] = Ceze.p[i][j][k] * Ez.p[i][j][k] + (Hy.p[i][j][k] - Hy.p[i - 1][j][k]) * Cezhy.p[i][j][k] +
                        (Hx.p[i][j][k] - Hx.p[i][j - 1][k]) * Cezhx.p[i][j][k];
            }
        }
    }
}

void fdtd::initCoeficients() {
    //////////////////////////
    // E Field coefficients
    //////////////////////////
    for (unsigned i = 0; i < Ex.nx; i++) {
        for (unsigned j = 0; j < Ex.ny; j++) {
            for (unsigned k = 0; k < Ex.nz; k++) {
                Cexhy.p[i][j][k] = -dt / eps_0 / dz;
                Cexhz.p[i][j][k] = dt / eps_0 / dy;
                Cexe.p[i][j][k] = 1;
            }
        }
    }
    for (unsigned i = 0; i < Ey.nx; i++) {
        for (unsigned j = 0; j < Ey.ny; j++) {
            for (unsigned k = 0; k < Ey.nz; k++) {
                Ceyhx.p[i][j][k] = dt / eps_0 / dz;
                Ceyhz.p[i][j][k] = -dt / eps_0 / dx;
                Ceye.p[i][j][k] = 1;
            }
        }
    }
    for (unsigned i = 0; i < Ez.nx; i++) {
        for (unsigned j = 0; j < Ez.ny; j++) {
            for (unsigned k = 0; k < Ez.nz; k++) {
                Cezhx.p[i][j][k] = -dt / eps_0 / dy;
                Cezhy.p[i][j][k] = dt / eps_0 / dx;
                Ceze.p[i][j][k] = 1;
            }
        }
    }
    //////////////////////////
    // M Field coefficients
    //////////////////////////
    for (unsigned i = 0; i < Hx.nx; i++) {
        for (unsigned j = 0; j < Hx.ny; j++) {
            for (unsigned k = 0; k < Hx.nz; k++) {
                Chxey.p[i][j][k] = dt / mu_0 / dz;
                Chxez.p[i][j][k] = -dt / mu_0 / dy;
                Chxh.p[i][j][k] = 1;
            }
        }
    }
    for (unsigned i = 0; i < Hy.nx; i++) {
        for (unsigned j = 0; j < Hy.ny; j++) {
            for (unsigned k = 0; k < Hy.nz; k++) {
                Chyex.p[i][j][k] = -dt / mu_0 / dz;
                Chyez.p[i][j][k] = dt / mu_0 / dx;
                Chyh.p[i][j][k] = 1;
            }
        }
    }
    for (unsigned i = 0; i < Hz.nx; i++) {
        for (unsigned j = 0; j < Hz.ny; j++) {
            for (unsigned k = 0; k < Hz.nz; k++) {
                Chzex.p[i][j][k] = dt / mu_0 / dy;
                Chzey.p[i][j][k] = -dt / mu_0 / dx;
                Chzh.p[i][j][k] = 1;
            }
        }
    }
}

void fdtd::updateMagneitcFields() {
    updateHx();
    updateHy();
    updateHz();
}

void fdtd::updateElectricFields() {
    updateEx();
    updateEy();
    updateEz();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// END OF PROGRAM CPMLFDTD3D
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

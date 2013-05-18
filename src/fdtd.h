#pragma once

#include "cpml.h"

#ifndef WITH_DENSITY
//#define WITH_DENSITY
#endif

class fdtd {
public:
#ifdef WITH_DENSITY
    fdtd(unsigned _totalTimeSteps = 500, unsigned _imax = 40, unsigned _jmax = 40, unsigned _kmax = 26,
            MyDataF _tw = 53.0e-12, MyDataF _dx = 1e-3, MyDataF _dy = 1e-3, MyDataF _dz = 1e-3,
            MyDataF _amp = 1000, unsigned _savemodulus = 10, unsigned _ksource = 12,
            unsigned _m = 3, unsigned _ma = 1, unsigned pmlw = 6, unsigned _nmatrial = 50, unsigned _neGrid = 16);

    void SetPlasmaVar(MyDataF _rei, MyDataF _vm, MyDataF _p, int _ftype);
#else
    fdtd(unsigned _totalTimeSteps = 500, unsigned _imax = 40, unsigned _jmax = 40, unsigned _kmax = 26,
            MyDataF _tw = 53.0e-12, MyDataF _dx = 1e-3, MyDataF _dy = 1e-3, MyDataF _dz = 1e-3,
            MyDataF _amp = 1000, unsigned _savemodulus = 10, unsigned _ksource = 12,
            unsigned _m = 3, unsigned _ma = 1, unsigned pmlw = 6, unsigned _nmatrial = 50);
#endif
    ~fdtd(void);

    static const int SOURCE_GAUSSIAN = GAUSSIAN_WAVE_TYPE;
    static const int SOURCE_SINE = SINE_WAVE_TYPE;
    static const int SOURCE_DERIVE_GAUSSIAN = DERIVE_GAUSSIAN_TYPE;
    static const int SOURCE_ZERO = ZERO_TYPE;

    //Function prototype definitions
    void initialize(); //Memory initialization
    void intSourceSinePulse(MyDataF t_0, MyDataF omega_, MyDataF tUp, MyDataF tDown, MyDataF amptidute);
    void setUp(); //Coefficients, parameters etc will get computed
    void compute(); //E & H Field update equation
    void StartUp();
    void SetSineSource(MyDataF omega_);
    void setSourceType(int sourceType);
#ifdef MATLAB_SIMULATION
    int initMatlabSimulation();
    void doMatlabSimulation();
    void finishMatlabSimulation();
#endif
private:
    void buildObject(); //Creates the object geometry
    void yeeCube(unsigned, unsigned, unsigned, unsigned); //Sets material properties to a cell
    void writeField(unsigned); //Writes output
    void buildSphere(); //Builds a spherical object
    void buildDipole(); //Builds a dipole
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

    /////////////////////////////////////////////////////////
    // SOURCE PARAMETRS
    ////////////////////////////////////////////////////////
    MyDataF t_up; // up bound of Pulse for Sine Pulse
    MyDataF t_down; // down bound of Pulse for Sine Pulse

    //Max number of materials allowed
    unsigned numMaterials;
#ifdef WITH_DENSITY
    //how many fine grids per coarse grid 
    unsigned neGrid;
    //initial plasma value
    MyDataF Ne0;
#endif
    // source type
    int srcType;

    //permittivity, permeability and conductivity of different materials
    MyDataF *epsilon;
    MyDataF *sigma;
    MyDataF *mu;


    MyDataF t0; //delay
    MyDataF source; //Differentiated Gaussian source    
    MyDataF omega; //wave angle frequency

    //  Specify the dipole Boundaries(A cuboidal rode- NOT as a cylinder)
    unsigned istart, iend, jstart;
    unsigned jend, kstart, kend;

    // source position
    unsigned isp, jsp, ksp;

    // H & E Field components
    data3d<MyDataF> Hx;
    data3d<MyDataF> Hy;
    data3d<MyDataF> Hz;
    data3d<MyDataF> Ex;
    data3d<MyDataF> Ey;
    data3d<MyDataF> Ez;

    data3d<unsigned> ID1; //medium definition array for Ex
    data3d<unsigned> ID2; //medium definition array for Ey
    data3d<unsigned> ID3; //medium definition array for Ez

    //E field update coefficients
    MyDataF *CA;
    MyDataF *CB;

    //H field update coefficients
    MyDataF DA;
    MyDataF DB;

    data3d<MyDataF> Cexe, Ceye, Ceze;
    data3d<MyDataF> Chxh, Chyh, Chzh;
    data3d<MyDataF> Cexhy, Ceyhz, Cezhx, Cexhz, Ceyhx, Cezhy;
    data3d<MyDataF> Chxey, Chyez, Chzex, Chxez, Chyex, Chzey;
    void initCoeficients();

#ifdef WITH_DENSITY

    int niutype;
    //Fine Grid size
    MyDataF dsf;
    //time step of plasma
    MyDataF dtf;
    int neSkipStep;

    // tempary that offen used
    MyDataF half_dt;
    MyDataF half_e;
    MyDataF Coeff_velocity;
    unsigned halfNeGrid;

    //plasma variables
    MyDataF vm; //collision frequency
    MyDataF p; // air pressure
    MyDataF De;
    MyDataF Da;
    MyDataF rei;
    MyDataF mu_i;
    MyDataF mu_e;
    MyDataF gamma;
    MyDataF a;
    MyDataF alpha;

    // update coefficients
    //    MyDataF Chxey,Chxez,Chyez,Chyex,Chzex,Chzey;
    data3d<MyDataF> Cezvz, Ceyvy, Cexvx;    
    data3d<MyDataF> Cvxex_guassian, Cvyey_guassian, Cvzez_guassian;
    MyDataF Cvxex, Cvyey, Cvzez;

    //Plasma
    data3d<MyDataF> Ne, Ne_pre;
    //
    data3d<MyDataF> Erms;
    // Beta
    data3d<MyDataF> beta;

    // collision frequency
    data3d<MyDataF> Nu_c;
    //
    data3d<MyDataF> Vx;
    data3d<MyDataF> Vy;
    data3d<MyDataF> Vz;

    //initials
    void initCoeff();
    void initDensity();
    void createCoeff();
    void updateCoeff();
    void updateBeta();

    // Erms or Eeff operation
    void IntegerEeff();
    void UpdateErms(void);
    void updateCollisionFrequency();
    void InterpErms();
    void UpdateDensity(void);
    void UpdateVeloity(void);
    void WallCircleBound(data3d<MyDataF> &stru);
#endif
    void updateHx();
    void updateHy();
    void updateHz();
    void updateMagneitcFields();
    void updateEx();
    void updateEy();
    void updateEz();
    void updateElectricAndVeloityFields();
    void updateSource(unsigned n);
    cpml<MyDataF> pml;

};



/* 
 * File:   cpml.h
 * Author: skiloop
 *
 * Created on April 3, 2013, 8:41 AM
 */

#ifndef CPML_H
#define	CPML_H

#include "datastruct.h"
#include <math.h>

#ifndef COMMON_H
// some common constants
const double C = 2.99792458E8; // speed of light in free space
const double me = 9.110e-31; // electricity mass
const double e = 1.602e-19; // electricity charge
const double mu_0 = 4.0 * M_PI * 1.0E-7;
const double eps_0 = 1.0 / (C * C * mu_0);
const double M_PI_TWO = M_PI * 2;
#endif
const double Mu0DivEps0 = mu_0 / eps_0;

#ifdef _OPENMP
#include <omp.h>
extern int thread_count;
#endif

/***************************************************************
 * FDTD region 
 * Field    |    x     |    y    |   z
 * ----------------------------------------------
 * Ez       |   I      |    J    |   K-1
 * Ex       |   I-1    |    J    |   K
 * Ey       |   I      |   J-1   |   K
 * Hx       |   I      |   J-1   |   K-1
 * Hy       |   I-1    |    J    |   K-1
 * Hz       |   I-1    |   J-1   |   K
 * ----------------------------------------------
 */
template<class T>
class cpml {
public:
    /**
     * default constructor
     */
    cpml();

    /**
     * constructor
     * @param width_xn width of cpml layer on x negetive side
     * @param width_xp width of cpml layer on x positive side
     * @param width_yn width of cpml layer on y negetive side
     * @param width_yp width of cpml layer on y positive side
     * @param width_zn width of cpml layer on z negetive side
     * @param width_zp width of cpml layer on z positive side
     * @param imax max number of cells in x direction
     * @param jmax max number of cells in y direction
     * @param kmax max number of cells in z direction
     * @param pmlOrder pml order
     * @param sigmaMax 
     * @param kappaMax
     * @param alphaMax
     * @param epsR
     */
    cpml(unsigned short width_xn, unsigned short width_xp, unsigned short width_yn,
            unsigned short width_yp, unsigned short width_zn, unsigned short width_zp,
            unsigned imax, unsigned jmax, unsigned kmax);

    /**
     * constructor
     * @param cpmlWidth width of cpml layer on all sides
     * @param imax max number of cells in x direction
     * @param jmax max number of cells in y direction
     * @param kmax max number of cells in z direction
     * @param pmlOrder pml order
     * @param sigmaMax 
     * @param kappaMax
     * @param alphaMax
     * @param epsR
     */
    cpml(unsigned short cpmlWidth, unsigned imax, unsigned jmax, unsigned kmax);

    /**
     * 
     * @param orig
     */
    cpml(const cpml& orig);

    /**
     * deconstructor
     */
    virtual ~cpml();
private:
    // cpml layer width on each size
    unsigned short n_cpml_xn;
    unsigned short n_cpml_xp;
    unsigned short n_cpml_yn;
    unsigned short n_cpml_yp;
    unsigned short n_cpml_zn;
    unsigned short n_cpml_zp;

    // cpml flags
    bool is_cpml_xn;
    bool is_cpml_xp;
    bool is_cpml_yn;
    bool is_cpml_yp;
    bool is_cpml_zn;
    bool is_cpml_zp;

    //============================================
    // cpml arrays
    //============================================
    // xn arrays
    data3d<T> Psi_eyx_xn;
    data3d<T> Psi_ezx_xn;
    data3d<T> Psi_hyx_xn;
    data3d<T> Psi_hzx_xn;
    // xp arrays
    data3d<T> Psi_eyx_xp;
    data3d<T> Psi_ezx_xp;
    data3d<T> Psi_hyx_xp;
    data3d<T> Psi_hzx_xp;
    // yn arrays
    data3d<T> Psi_exy_yn;
    data3d<T> Psi_ezy_yn;
    data3d<T> Psi_hxy_yn;
    data3d<T> Psi_hzy_yn;
    // yp arrays
    data3d<T> Psi_exy_yp;
    data3d<T> Psi_ezy_yp;
    data3d<T> Psi_hxy_yp;
    data3d<T> Psi_hzy_yp;
    // zn arrays
    data3d<T> Psi_exz_zn;
    data3d<T> Psi_eyz_zn;
    data3d<T> Psi_hxz_zn;
    data3d<T> Psi_hyz_zn;
    // zp arrays
    //    data3d<T> Psi_exz_zp;
    data3d<T> Psi_eyz_zp;
    data3d<T> Psi_hxz_zp;
    data3d<T> Psi_hyz_zp;
    //============================================
    // cpml coefficient arrays
    //============================================
    // xn arrays
    data3d<T> CPsi_eyx_xn;
    data3d<T> CPsi_ezx_xn;
    data3d<T> CPsi_hyx_xn;
    data3d<T> CPsi_hzx_xn;
    // xp arrays
    data3d<T> CPsi_eyx_xp;
    data3d<T> CPsi_ezx_xp;
    data3d<T> CPsi_hyx_xp;
    data3d<T> CPsi_hzx_xp;
    // yn arrays
    data3d<T> CPsi_exy_yn;
    data3d<T> CPsi_ezy_yn;
    data3d<T> CPsi_hxy_yn;
    data3d<T> CPsi_hzy_yn;
    // yp arrays
    data3d<T> CPsi_exy_yp;
    data3d<T> CPsi_ezy_yp;
    data3d<T> CPsi_hxy_yp;
    data3d<T> CPsi_hzy_yp;
    // zn arrays
    data3d<T> CPsi_exz_zn;
    data3d<T> CPsi_eyz_zn;
    data3d<T> CPsi_hxz_zn;
    data3d<T> CPsi_hyz_zn;
    // zp arrays
    data3d<T> CPsi_exz_zp;
    data3d<T> CPsi_eyz_zp;
    data3d<T> CPsi_hxz_zp;
    data3d<T> CPsi_hyz_zp;

    // a and b for cpml to update Psi
    // x direction
    data1d<T> cpml_a_ex_xn;
    data1d<T> cpml_b_ex_xn;
    data1d<T> cpml_a_mx_xn;
    data1d<T> cpml_b_mx_xn;
    data1d<T> cpml_a_ex_xp;
    data1d<T> cpml_b_ex_xp;
    data1d<T> cpml_a_mx_xp;
    data1d<T> cpml_b_mx_xp;
    // y direction
    data1d<T> cpml_a_ey_yn;
    data1d<T> cpml_b_ey_yn;
    data1d<T> cpml_a_my_yn;
    data1d<T> cpml_b_my_yn;
    data1d<T> cpml_a_ey_yp;
    data1d<T> cpml_b_ey_yp;
    data1d<T> cpml_a_my_yp;
    data1d<T> cpml_b_my_yp;
    // z direction
    data1d<T> cpml_a_ez_zn;
    data1d<T> cpml_b_ez_zn;
    data1d<T> cpml_a_mz_zn;
    data1d<T> cpml_b_mz_zn;
    data1d<T> cpml_a_ez_zp;
    data1d<T> cpml_b_ez_zp;
    data1d<T> cpml_a_mz_zp;
    data1d<T> cpml_b_mz_zp;

    //=================================================
    // PUBLIC interface 
    //=================================================
public:
    data3d<T> Psi_exz_zp;
    /**
     * update electric fields in CPML region
     * @param Ex
     * @param Ey
     * @param Ez
     * @param Hx
     * @param Hy
     * @param Hz
     */
    void updateCPML_E_Fields(data3d<T> &Ex, data3d<T>& Ey, data3d<T> &Ez,
            const data3d<T> & Hx, const data3d<T> & Hy, const data3d<T> & Hz);
    /**
     * update magnetic fields in CPML region
     * @param Hx
     * @param Hy
     * @param Hz
     * @param Ex
     * @param Ey
     * @param Ez
     */
    void updateCPML_M_Fields(data3d<T> &Hx, data3d<T>& Hy, data3d<T> &Hz,
            const data3d<T> & Ex, const data3d<T> & Ey, const data3d<T> & Ez);

    /**
     * 
     * @param pmlWidth
     */
    void setCPMLRegion(short pmlWidth);
    /**
     * 
     * @param width_xn
     * @param width_xp
     * @param width_yn
     * @param width_yp
     * @param width_zn
     * @param width_zp
     */
    void setCPMLRegion(short width_xn, short width_xp, short width_yn, short width_yp, short width_zn, short width_zp);
    /**
     * 
     * @param nx
     * @param ny
     * @param nz
     */
    void createCPMLArrays(unsigned nx, unsigned ny, unsigned nz);
    /**
     * 
     * @param pmlOrder
     * @param sigmaMax
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dx
     * @param dy
     * @param dz
     * @param Ceyhz
     * @param Cezhy
     * @param Chyez
     * @param Chzey
     * @param Cexhz
     * @param Cezhx
     * @param Chxez
     * @param Chzex
     * @param Ceyhx
     * @param Cexhy
     * @param Chyex
     * @param Chxey
     */
    void initCoefficientArrays(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dx, T dy, T dz,
            data3d<T>&Ceyhz, data3d<T>&Cezhy, data3d<T>&Chyez, data3d<T>&Chzey,
            data3d<T>&Cexhz, data3d<T>&Cezhx, data3d<T>&Chxez, data3d<T>&Chzex,
            data3d<T>&Ceyhx, data3d<T>&Cexhy, data3d<T>&Chyex, data3d<T>&Chxey);
    //=======================================================
    // private functions
    //=======================================================
private:
    /**
     * 
     * @param pmlOrder
     * @param sigmaMax
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dx
     * @param Ceyhz
     * @param Cezhy
     * @param Chyez
     * @param Chzey
     */
    void initCoefficientArraysXN(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dx, data3d<T>&Ceyhz, data3d<T>&Cezhy, data3d<T>&Chyez, data3d<T>&Chzey);

    /**
     * 
     * @param pmlOrder
     * @param sigmaMax
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dx
     * @param Ceyhz
     * @param Cezhy
     * @param Chyez
     * @param Chzey
     */
    void initCoefficientArraysXP(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dx, data3d<T>&Ceyhz, data3d<T>&Cezhy, data3d<T>&Chyez, data3d<T>&Chzey);

    /**
     * 
     * @param pmlOrder
     * @param sigmaMax
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dy
     * @param Cexhz
     * @param Cezhx
     * @param Chxez
     * @param Chzex
     */
    void initCoefficientArraysYN(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dy, data3d<T>&Cexhz, data3d<T>&Cezhx, data3d<T>&Chxez, data3d<T>&Chzex);
    /**
     * 
     * @param pmlOrder
     * @param sigmaMax
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dy
     * @param Cexhz
     * @param Cezhx
     * @param Chxez
     * @param Chzex
     */
    void initCoefficientArraysYP(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dy, data3d<T>&Cexhz, data3d<T>&Cezhx, data3d<T>&Chxez, data3d<T>&Chzex);
    /**
     * 
     * @param pmlOrder
     * @param sigmaMax
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dz
     * @param Ceyhx
     * @param Cexhy
     * @param Chyex
     * @param Chxey
     */
    void initCoefficientArraysZN(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dz, data3d<T>&Ceyhx, data3d<T>&Cexhy, data3d<T>&Chyex, data3d<T>&Chxey);
    /**
     * 
     * @param pmlOrder
     * @param sigmaMax
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dz
     * @param Ceyhx
     * @param Cexhy
     * @param Chyex
     * @param Chxey
     */
    void initCoefficientArraysZP(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dz, data3d<T>&Ceyhx, data3d<T>&Cexhy, data3d<T>&Chyex, data3d<T>&Chxey);

    /**
     * 
     * @param Hx
     * @param Hy
     * @param Hz
     */
    void updatePsiForEFields(const data3d<T>& Hx, const data3d<T>& Hy, const data3d<T>& Hz);
    /**
     * 
     * @param Ex
     * @param Ey
     * @param Ez
     */
    void updatePsiForMFields(const data3d<T>& Ex, const data3d<T>& Ey, const data3d<T>& Ez);

    void updatePsi_eyz_zn(const data3d<T>& Hx);
    void updatePsi_ezy_yn(const data3d<T>& Hx);
    void updatePsi_exz_zn(const data3d<T>& Hy);
    void updatePsi_ezx_xn(const data3d<T>& Hy);
    void updatePsi_eyx_xn(const data3d<T>& Hz);
    void updatePsi_exy_yn(const data3d<T>& Hz);

    void updatePsi_hyz_zn(const data3d<T>& Ex);
    void updatePsi_hzy_yn(const data3d<T>& Ex);
    void updatePsi_hxz_zn(const data3d<T>& Ey);
    void updatePsi_hzx_xn(const data3d<T>& Ey);
    void updatePsi_hyx_xn(const data3d<T>& Ez);
    void updatePsi_hxy_yn(const data3d<T>& Ez);

    void updatePsi_eyz_zp(const data3d<T>& Hx);
    void updatePsi_ezy_yp(const data3d<T>& Hx);
    void updatePsi_exz_zp(const data3d<T>& Hy);
    void updatePsi_ezx_xp(const data3d<T>& Hy);
    void updatePsi_eyx_xp(const data3d<T>& Hz);
    void updatePsi_exy_yp(const data3d<T>& Hz);

    void updatePsi_hyz_zp(const data3d<T>& Ex);
    void updatePsi_hzy_yp(const data3d<T>& Ex);
    void updatePsi_hxz_zp(const data3d<T>& Ey);
    void updatePsi_hzx_xp(const data3d<T>& Ey);
    void updatePsi_hyx_xp(const data3d<T>& Ez);
    void updatePsi_hxy_yp(const data3d<T>& Ez);

    void updateEFieldCPML_x(data3d<T>&Ey, data3d<T>&Ez);
    void updateEFieldCPML_y(data3d<T>&Ex, data3d<T>&Ez);
    void updateEFieldCPML_z(data3d<T>&Ex, data3d<T>&Ey);

    void updateMFieldCPML_x(data3d<T>&Hy, data3d<T>&Hz);
    void updateMFieldCPML_y(data3d<T>&Hx, data3d<T>&Hz);
    void updateMFieldCPML_z(data3d<T>&Hx, data3d<T>&Hy);
};

template<class T>
cpml<T>::cpml()
: n_cpml_xn(0)
, n_cpml_xp(0)
, n_cpml_yn(0)
, n_cpml_yp(0)
, n_cpml_zn(0)
, n_cpml_zp(0)
, is_cpml_xn(false)
, is_cpml_xp(false)
, is_cpml_yn(false)
, is_cpml_yp(false)
, is_cpml_zn(false)
, is_cpml_zp(false) {
}

template<class T> cpml<T>::cpml(unsigned short width_xn, unsigned short width_xp,
        unsigned short width_yn, unsigned short width_yp,
        unsigned short width_zn, unsigned short width_zp,
        unsigned imax, unsigned jmax, unsigned kmax) {
    setCPMLRegion(width_xn, width_xp, width_yn, width_yp, width_zn, width_zp);
    createCPMLArrays(imax, jmax, kmax);
}

template<class T> cpml<T>::cpml(unsigned short cpmlWidth, unsigned imax, unsigned jmax, unsigned kmax) {
    setCPMLRegion(cpmlWidth);
    createCPMLArrays(imax, jmax, kmax);
}

template<class T>
cpml<T>::cpml(const cpml& orig)
: n_cpml_xn(orig.n_cpml_xn)
, n_cpml_xp(orig.n_cpml_xp)
, n_cpml_yn(orig.n_cpml_yn)
, n_cpml_yp(orig.n_cpml_yp)
, n_cpml_zn(orig.n_cpml_zn)
, n_cpml_zp(orig.n_cpml_zp)
, is_cpml_xn(orig.is_cpml_xn)
, is_cpml_xp(orig.is_cpml_xp)
, is_cpml_yn(orig.is_cpml_yn)
, is_cpml_yp(orig.is_cpml_yp)
, is_cpml_zn(orig.is_cpml_zn)
, is_cpml_zp(orig.is_cpml_zp) {
    //============================================
    // cpml arrays
    //============================================
    // xn arrays
    Psi_eyx_xn.CreateStruct(orig. Psi_eyx_xn);
    Psi_eyx_xn.Backupdata3d(orig. Psi_eyx_xn);
    Psi_ezx_xn.CreateStruct(orig. Psi_ezx_xn);
    Psi_ezx_xn.Backupdata3d(orig. Psi_ezx_xn);
    Psi_hyx_xn.CreateStruct(orig. Psi_hyx_xn);
    Psi_hyx_xn.Backupdata3d(orig. Psi_hyx_xn);
    Psi_hzx_xn.CreateStruct(orig. Psi_hzx_xn);
    Psi_hzx_xn.Backupdata3d(orig. Psi_hzx_xn);
    // xp arrays
    Psi_eyx_xp.CreateStruct(orig. Psi_eyx_xp);
    Psi_eyx_xp.Backupdata3d(orig. Psi_eyx_xp);
    Psi_ezx_xp.CreateStruct(orig. Psi_ezx_xp);
    Psi_ezx_xp.Backupdata3d(orig. Psi_ezx_xp);
    Psi_hyx_xp.CreateStruct(orig. Psi_hyx_xp);
    Psi_hyx_xp.Backupdata3d(orig. Psi_hyx_xp);
    Psi_hzx_xp.CreateStruct(orig. Psi_hzx_xp);
    Psi_hzx_xp.Backupdata3d(orig. Psi_hzx_xp);
    // yn arrays
    Psi_exy_yn.CreateStruct(orig. Psi_exy_yn);
    Psi_exy_yn.Backupdata3d(orig. Psi_exy_yn);
    Psi_ezy_yn.CreateStruct(orig. Psi_ezy_yn);
    Psi_ezy_yn.Backupdata3d(orig. Psi_ezy_yn);
    Psi_hxy_yn.CreateStruct(orig. Psi_hxy_yn);
    Psi_hxy_yn.Backupdata3d(orig. Psi_hxy_yn);
    Psi_hzy_yn.CreateStruct(orig. Psi_hzy_yn);
    Psi_hzy_yn.Backupdata3d(orig. Psi_hzy_yn);
    // yp arrays
    Psi_exy_yp.CreateStruct(orig. Psi_exy_yp);
    Psi_exy_yp.Backupdata3d(orig. Psi_exy_yp);
    Psi_ezy_yp.CreateStruct(orig. Psi_ezy_yp);
    Psi_ezy_yp.Backupdata3d(orig. Psi_ezy_yp);
    Psi_hxy_yp.CreateStruct(orig. Psi_hxy_yp);
    Psi_hxy_yp.Backupdata3d(orig. Psi_hxy_yp);
    Psi_hzy_yp.CreateStruct(orig. Psi_hzy_yp);
    Psi_hzy_yp.Backupdata3d(orig. Psi_hzy_yp);
    // zn arrays
    Psi_exz_zn.CreateStruct(orig. Psi_exz_zn);
    Psi_exz_zn.Backupdata3d(orig. Psi_exz_zn);
    Psi_eyz_zn.CreateStruct(orig. Psi_eyz_zn);
    Psi_eyz_zn.Backupdata3d(orig. Psi_eyz_zn);
    Psi_hxz_zn.CreateStruct(orig. Psi_hxz_zn);
    Psi_hxz_zn.Backupdata3d(orig. Psi_hxz_zn);
    Psi_hyz_zn.CreateStruct(orig. Psi_hyz_zn);
    Psi_hyz_zn.Backupdata3d(orig. Psi_hyz_zn);
    // zp arrays
    Psi_exz_zp.CreateStruct(orig. Psi_exz_zp);
    Psi_exz_zp.Backupdata3d(orig. Psi_exz_zp);
    Psi_eyz_zp.CreateStruct(orig. Psi_eyz_zp);
    Psi_eyz_zp.Backupdata3d(orig. Psi_eyz_zp);
    Psi_hxz_zp.CreateStruct(orig. Psi_hxz_zp);
    Psi_hxz_zp.Backupdata3d(orig. Psi_hxz_zp);
    Psi_hyz_zp.CreateStruct(orig. Psi_hyz_zp);
    Psi_hyz_zp.Backupdata3d(orig. Psi_hyz_zp);
    //============================================
    // cpml coefficient arrays
    //============================================
    // xn arrays
    CPsi_eyx_xn.CreateStruct(orig. CPsi_eyx_xn);
    CPsi_eyx_xn.Backupdata3d(orig. CPsi_eyx_xn);
    CPsi_ezx_xn.CreateStruct(orig. CPsi_ezx_xn);
    CPsi_ezx_xn.Backupdata3d(orig. CPsi_ezx_xn);
    CPsi_hyx_xn.CreateStruct(orig. CPsi_hyx_xn);
    CPsi_hyx_xn.Backupdata3d(orig. CPsi_hyx_xn);
    CPsi_hzx_xn.CreateStruct(orig. CPsi_hzx_xn);
    CPsi_hzx_xn.Backupdata3d(orig. CPsi_hzx_xn);
    // xp arrays
    CPsi_eyx_xp.CreateStruct(orig. CPsi_eyx_xp);
    CPsi_eyx_xp.Backupdata3d(orig. CPsi_eyx_xp);
    CPsi_ezx_xp.CreateStruct(orig. CPsi_ezx_xp);
    CPsi_ezx_xp.Backupdata3d(orig. CPsi_ezx_xp);
    CPsi_hyx_xp.CreateStruct(orig. CPsi_hyx_xp);
    CPsi_hyx_xp.Backupdata3d(orig. CPsi_hyx_xp);
    CPsi_hzx_xp.CreateStruct(orig. CPsi_hzx_xp);
    CPsi_hzx_xp.Backupdata3d(orig. CPsi_hzx_xp);
    // yn arrays
    CPsi_exy_yn.CreateStruct(orig. CPsi_exy_yn);
    CPsi_exy_yn.Backupdata3d(orig. CPsi_exy_yn);
    CPsi_ezy_yn.CreateStruct(orig. CPsi_ezy_yn);
    CPsi_ezy_yn.Backupdata3d(orig. CPsi_ezy_yn);
    CPsi_hxy_yn.CreateStruct(orig. CPsi_hxy_yn);
    CPsi_hxy_yn.Backupdata3d(orig. CPsi_hxy_yn);
    CPsi_hzy_yn.CreateStruct(orig. CPsi_hzy_yn);
    CPsi_hzy_yn.Backupdata3d(orig. CPsi_hzy_yn);
    // yp arrays
    CPsi_exy_yp.CreateStruct(orig. CPsi_exy_yp);
    CPsi_exy_yp.Backupdata3d(orig. CPsi_exy_yp);
    CPsi_ezy_yp.CreateStruct(orig. CPsi_ezy_yp);
    CPsi_ezy_yp.Backupdata3d(orig. CPsi_ezy_yp);
    CPsi_hxy_yp.CreateStruct(orig. CPsi_hxy_yp);
    CPsi_hxy_yp.Backupdata3d(orig. CPsi_hxy_yp);
    CPsi_hzy_yp.CreateStruct(orig. CPsi_hzy_yp);
    CPsi_hzy_yp.Backupdata3d(orig. CPsi_hzy_yp);
    // zn arrays
    CPsi_exz_zn.CreateStruct(orig. CPsi_exz_zn);
    CPsi_exz_zn.Backupdata3d(orig. CPsi_exz_zn);
    CPsi_eyz_zn.CreateStruct(orig. CPsi_eyz_zn);
    CPsi_eyz_zn.Backupdata3d(orig. CPsi_eyz_zn);
    CPsi_hxz_zn.CreateStruct(orig. CPsi_hxz_zn);
    CPsi_hxz_zn.Backupdata3d(orig. CPsi_hxz_zn);
    CPsi_hyz_zn.CreateStruct(orig. CPsi_hyz_zn);
    CPsi_hyz_zn.Backupdata3d(orig. CPsi_hyz_zn);
    // zp arrays
    CPsi_exz_zp.CreateStruct(orig. CPsi_exz_zp);
    CPsi_exz_zp.Backupdata3d(orig. CPsi_exz_zp);
    CPsi_eyz_zp.CreateStruct(orig. CPsi_eyz_zp);
    CPsi_eyz_zp.Backupdata3d(orig. CPsi_eyz_zp);
    CPsi_hxz_zp.CreateStruct(orig. CPsi_hxz_zp);
    CPsi_hxz_zp.Backupdata3d(orig. CPsi_hxz_zp);
    CPsi_hyz_zp.CreateStruct(orig. CPsi_hyz_zp);
    CPsi_hyz_zp.Backupdata3d(orig. CPsi_hyz_zp);
}

template<class T>
cpml<T>::~cpml() {
}

template<class T>
void cpml<T>::createCPMLArrays(unsigned nx, unsigned ny, unsigned nz) {
    //============================================
    // cpml arrays
    //============================================
    // xn arrays
    unsigned nxp1 = nx + 1;
    unsigned nyp1 = ny + 1;
    unsigned nzp1 = nz + 1;
    if (is_cpml_xn) {
        // x direction
        cpml_a_ex_xn.CreateStruct(n_cpml_xn, 0.0);
        cpml_b_ex_xn.CreateStruct(n_cpml_xn, 0.0);
        cpml_a_mx_xn.CreateStruct(n_cpml_xn, 0.0);
        cpml_b_mx_xn.CreateStruct(n_cpml_xn, 0.0);
        Psi_eyx_xn.CreateStruct(n_cpml_xn, ny, nzp1, 0.0);
        Psi_ezx_xn.CreateStruct(n_cpml_xn, nyp1, nz, 0.0);
        Psi_hyx_xn.CreateStruct(n_cpml_xn, nyp1, nz, 0.0);
        Psi_hzx_xn.CreateStruct(n_cpml_xn, ny, nzp1, 0.0);
        CPsi_eyx_xn.CreateStruct(n_cpml_xn, ny, nzp1, 0.0);
        CPsi_ezx_xn.CreateStruct(n_cpml_xn, nyp1, nz, 0.0);
        CPsi_hyx_xn.CreateStruct(n_cpml_xn, nyp1, nz, 0.0);
        CPsi_hzx_xn.CreateStruct(n_cpml_xn, ny, nzp1, 0.0);
    }
    // xp arrays
    if (is_cpml_xp) {
        cpml_a_ex_xp.CreateStruct(n_cpml_xp, 0.0);
        cpml_b_ex_xp.CreateStruct(n_cpml_xp, 0.0);
        cpml_a_mx_xp.CreateStruct(n_cpml_xp, 0.0);
        cpml_b_mx_xp.CreateStruct(n_cpml_xp, 0.0);
        Psi_eyx_xp.CreateStruct(n_cpml_xp, ny, nzp1, 0.0);
        Psi_ezx_xp.CreateStruct(n_cpml_xp, nyp1, nz, 0.0);
        Psi_hyx_xp.CreateStruct(n_cpml_xp, nyp1, nz, 0.0);
        Psi_hzx_xp.CreateStruct(n_cpml_xp, ny, nzp1, 0.0);
        CPsi_eyx_xp.CreateStruct(n_cpml_xp, ny, nzp1, 0.0);
        CPsi_ezx_xp.CreateStruct(n_cpml_xp, nyp1, nz, 0.0);
        CPsi_hyx_xp.CreateStruct(n_cpml_xp, nyp1, nz, 0.0);
        CPsi_hzx_xp.CreateStruct(n_cpml_xp, ny, nzp1, 0.0);
    }
    // yn arrays
    if (is_cpml_yn) {
        // y direction
        cpml_a_ey_yn.CreateStruct(n_cpml_yn, 0.0);
        cpml_b_ey_yn.CreateStruct(n_cpml_yn, 0.0);
        cpml_a_my_yn.CreateStruct(n_cpml_yn, 0.0);
        cpml_b_my_yn.CreateStruct(n_cpml_yn, 0.0);
        Psi_exy_yn.CreateStruct(nx, n_cpml_yn, nzp1, 0.0);
        Psi_ezy_yn.CreateStruct(nxp1, n_cpml_yn, nz, 0.0);
        Psi_hxy_yn.CreateStruct(nxp1, n_cpml_yn, nz, 0.0);
        Psi_hzy_yn.CreateStruct(nx, n_cpml_yn, nzp1, 0.0);
        CPsi_exy_yn.CreateStruct(nx, n_cpml_yn, nzp1, 0.0);
        CPsi_ezy_yn.CreateStruct(nxp1, n_cpml_yn, nz, 0.0);
        CPsi_hxy_yn.CreateStruct(nxp1, n_cpml_yn, nz, 0.0);
        CPsi_hzy_yn.CreateStruct(nx, n_cpml_yn, nzp1, 0.0);
    }
    // yp arrays
    if (is_cpml_yp) {
        cpml_a_ey_yp.CreateStruct(n_cpml_yp, 0.0);
        cpml_b_ey_yp.CreateStruct(n_cpml_yp, 0.0);
        cpml_a_my_yp.CreateStruct(n_cpml_yp, 0.0);
        cpml_b_my_yp.CreateStruct(n_cpml_yp, 0.0);
        Psi_exy_yp.CreateStruct(nx, n_cpml_yp, nzp1, 0.0);
        Psi_ezy_yp.CreateStruct(nxp1, n_cpml_yp, nz, 0.0);
        Psi_hxy_yp.CreateStruct(nxp1, n_cpml_yp, nz, 0.0);
        Psi_hzy_yp.CreateStruct(nx, n_cpml_yp, nzp1, 0.0);
        CPsi_exy_yp.CreateStruct(nx, n_cpml_yp, nzp1, 0.0);
        CPsi_ezy_yp.CreateStruct(nxp1, n_cpml_yp, nz, 0.0);
        CPsi_hxy_yp.CreateStruct(nxp1, n_cpml_yp, nz, 0.0);
        CPsi_hzy_yp.CreateStruct(nx, n_cpml_yp, nzp1, 0.0);
    }
    // zn arrays
    if (is_cpml_zn) {
        // z direction
        cpml_a_ez_zn.CreateStruct(n_cpml_zn, 0.0);
        cpml_b_ez_zn.CreateStruct(n_cpml_zn, 0.0);
        cpml_a_mz_zn.CreateStruct(n_cpml_zn, 0.0);
        cpml_b_mz_zn.CreateStruct(n_cpml_zn, 0.0);
        Psi_exz_zn.CreateStruct(nx, nyp1, n_cpml_zn, 0.0);
        Psi_eyz_zn.CreateStruct(nxp1, ny, n_cpml_zn, 0.0);
        Psi_hxz_zn.CreateStruct(nxp1, ny, n_cpml_zn, 0.0);
        Psi_hyz_zn.CreateStruct(nx, nyp1, n_cpml_zn, 0.0);
        CPsi_exz_zn.CreateStruct(nx, nyp1, n_cpml_zn, 0.0);
        CPsi_eyz_zn.CreateStruct(nxp1, ny, n_cpml_zn, 0.0);
        CPsi_hxz_zn.CreateStruct(nxp1, ny, n_cpml_zn, 0.0);
        CPsi_hyz_zn.CreateStruct(nx, nyp1, n_cpml_zn, 0.0);
    }
    // zp arrays
    if (is_cpml_zp) {
        cpml_a_ez_zp.CreateStruct(n_cpml_zp, 0.0);
        cpml_b_ez_zp.CreateStruct(n_cpml_zp, 0.0);
        cpml_a_mz_zp.CreateStruct(n_cpml_zp, 0.0);
        cpml_b_mz_zp.CreateStruct(n_cpml_zp, 0.0);
        Psi_exz_zp.CreateStruct(nx, nyp1, n_cpml_zp, 0.0);
        Psi_eyz_zp.CreateStruct(nxp1, ny, n_cpml_zp, 0.0);
        Psi_hxz_zp.CreateStruct(nxp1, ny, n_cpml_zp, 0.0);
        Psi_hyz_zp.CreateStruct(nx, nyp1, n_cpml_zp, 0.0);
        CPsi_exz_zp.CreateStruct(nx, nyp1, n_cpml_zp, 0.0);
        CPsi_eyz_zp.CreateStruct(nxp1, ny, n_cpml_zp, 0.0);
        CPsi_hxz_zp.CreateStruct(nxp1, ny, n_cpml_zp, 0.0);
        CPsi_hyz_zp.CreateStruct(nx, nyp1, n_cpml_zp, 0.0);
    }
}

template<class T>
void cpml<T>::setCPMLRegion(short pmlWidth) {
    setCPMLRegion(pmlWidth, pmlWidth, pmlWidth, pmlWidth, pmlWidth, pmlWidth);
}

template<class T>
void cpml<T>::setCPMLRegion(short width_xn, short width_xp, short width_yn, short width_yp, short width_zn, short width_zp) {
    n_cpml_xn = width_xn;
    n_cpml_xp = width_xp;
    n_cpml_yn = width_yn;
    n_cpml_yp = width_yp;
    n_cpml_zn = width_zn;
    n_cpml_zp = width_zp;
    if (n_cpml_xn > 0) {
        is_cpml_xn = true;
    }
    if (n_cpml_xp > 0) {
        is_cpml_xp = true;
    }
    if (n_cpml_yn > 0) {
        is_cpml_yn = true;
    }
    if (n_cpml_yp > 0) {
        is_cpml_yp = true;
    }
    if (n_cpml_zn > 0) {
        is_cpml_zn = true;
    }
    if (n_cpml_zp > 0) {
        is_cpml_zp = true;
    }
}

template<class T>
void cpml<T>::initCoefficientArrays(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dx, T dy, T dz,
        data3d<T>&Ceyhz, data3d<T>&Cezhy, data3d<T>&Chyez, data3d<T>&Chzey,
        data3d<T>&Cexhz, data3d<T>&Cezhx, data3d<T>&Chxez, data3d<T>&Chzex,
        data3d<T>&Ceyhx, data3d<T>&Cexhy, data3d<T>&Chyex, data3d<T>&Chxey) {
    initCoefficientArraysXN(pmlOrder, sigmaMax, kappaMax, alphaMax, epsR, dt, dx, Ceyhz, Cezhy, Chyez, Chzey);
    initCoefficientArraysXP(pmlOrder, sigmaMax, kappaMax, alphaMax, epsR, dt, dx, Ceyhz, Cezhy, Chyez, Chzey);
    initCoefficientArraysYN(pmlOrder, sigmaMax, kappaMax, alphaMax, epsR, dt, dy, Cexhz, Cezhx, Chxez, Chzex);
    initCoefficientArraysYP(pmlOrder, sigmaMax, kappaMax, alphaMax, epsR, dt, dy, Cexhz, Cezhx, Chxez, Chzex);
    initCoefficientArraysZN(pmlOrder, sigmaMax, kappaMax, alphaMax, epsR, dt, dz, Ceyhx, Cexhy, Chyex, Chxey);
    initCoefficientArraysZP(pmlOrder, sigmaMax, kappaMax, alphaMax, epsR, dt, dz, Ceyhx, Cexhy, Chyex, Chxey);
}

template<class T>
void cpml<T>::initCoefficientArraysXN(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dx,
        data3d<T>&Ceyhz, data3d<T>&Cezhy, data3d<T>&Chyez, data3d<T>&Chzey) {
    if (is_cpml_xn) {
        T sigmaOpt = sigmaMax * (pmlOrder + 1) / (150 * M_PI * dx * sqrt(epsR));

        for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
            T rho_e = (n_cpml_xn - i - 0.75) / (T) n_cpml_xn;
            //T rho_e = (n_cpml_xn - i) /(T)n_cpml_xn;	
            T rho_m = (n_cpml_xn - i - 0.25) / (T) n_cpml_xn;
            //T rho_m = (n_cpml_xn - i) /(T)n_cpml_xn;
            T rho_e_pmlOrder = pow(fabs(rho_e), pmlOrder);
            T rho_m_pmlOrder = pow(fabs(rho_m), pmlOrder);
            T sigma_pex = sigmaOpt*rho_e_pmlOrder;
            T sigma_pmx = sigmaOpt * rho_m_pmlOrder;
            T kappa_ex = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            T kappa_mx = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            T alpha_ex = alphaMax * fabs(rho_e);
            T alpha_mx = alphaMax * fabs(rho_m);
            cpml_b_ex_xn.p[i] = exp((-dt / eps_0) * (sigma_pex / kappa_ex + alpha_ex));
            cpml_b_mx_xn.p[i] = exp((-dt / eps_0) * (sigma_pmx / kappa_mx + alpha_mx));
            if (rho_e != 0) {
                cpml_a_ex_xn.p[i] = 1 / dx * (cpml_b_ex_xn.p[i] - 1.0) * sigma_pex / (kappa_ex * (sigma_pex + kappa_ex * alpha_ex));
            } else {
                cpml_a_ex_xn.p[i] = 0;
            }
            if (rho_m != 0) {
                cpml_a_mx_xn.p[i] = 1 / dx * (cpml_b_mx_xn.p[i] - 1.0) * sigma_pmx / (kappa_mx * (sigma_pmx + kappa_mx * alpha_mx));
            } else {
                cpml_a_mx_xn.p[i] = 0;
            }
            for (unsigned j = 0; j < Psi_eyx_xn.ny; j++) {
                for (unsigned k = 0; k < Psi_eyx_xn.nz; k++) {
                    CPsi_eyx_xn.p[i][j][k] = Ceyhz.p[iplus][j][k] * dx;
                    Ceyhz.p[iplus][j][k] = Ceyhz.p[iplus][j][k] / kappa_ex;
                }
            }
            for (unsigned j = 0; j < Psi_ezx_xn.nz; j++) {
                for (unsigned k = 0; k < Psi_ezx_xn.ny; k++) {
                    CPsi_ezx_xn.p[i][k][j] = Cezhy.p[iplus][k][j] * dx;
                    Cezhy.p[iplus][k][j] = Cezhy.p[iplus][k][j] / kappa_ex;
                }
            }
            for (unsigned j = 0; j < Psi_hyx_xn.ny; j++) {
                for (unsigned k = 0; k < Psi_hyx_xn.nz; k++) {
                    CPsi_hyx_xn.p[i][j][k] = Chyez.p[i][j][k] * dx;
                    Chyez.p[i][j][k] = Chyez.p[i][j][k] / kappa_mx;
                }
            }
            for (unsigned j = 0; j < Psi_hzx_xn.nz; j++) {
                for (unsigned k = 0; k < Psi_hzx_xn.ny; k++) {
                    CPsi_hzx_xn.p[i][k][j] = Chzey.p[i][k][j] * dx;
                    Chzey.p[i][k][j] = Chzey.p[i][k][j] / kappa_mx;
                }
            }
            //T cev = 1 / (eps_0 / dt + sigma_pex / 2) / dx;
            //T chv = 1 / (mu_0 / dt + sigma_pmx / 2) / dx;
            //for (unsigned j = 0; j < Psi_eyx_xn.ny; j++) {
            //    for (unsigned k = 0; k < Psi_eyx_xn.nz; k++) {
            //        Ceyhz.p[iplus][j][k] = -cev;
            //        CPsi_eyx_xn.p[i][j][k] = Ceyhz.p[iplus][j][k] * dx;
            //        Ceyhz.p[iplus][j][k] = Ceyhz.p[iplus][j][k] / kappa_ex;
            //    }
            //}
            //for (unsigned j = 0; j < Psi_ezx_xn.nz; j++) {
            //    for (unsigned k = 0; k < Psi_ezx_xn.ny; k++) {
            //        Cezhy.p[iplus][k][j] = cev;
            //        CPsi_ezx_xn.p[i][k][j] = Cezhy.p[iplus][k][j] * dx;
            //        Cezhy.p[iplus][k][j] = Cezhy.p[iplus][k][j] / kappa_ex;
            //    }
            //}
            //for (unsigned j = 0; j < Psi_hyx_xn.ny; j++) {
            //    for (unsigned k = 0; k < Psi_hyx_xn.nz; k++) {
            //        Chyez.p[i][j][k] = chv;
            //        CPsi_hyx_xn.p[i][j][k] = Chyez.p[i][j][k] * dx;
            //        Chyez.p[i][j][k] = Chyez.p[i][j][k] / kappa_mx;
            //    }
            //}
            //for (unsigned j = 0; j < Psi_hzx_xn.nz; j++) {
            //    for (unsigned k = 0; k < Psi_hzx_xn.ny; k++) {
            //        Chzey.p[i][k][j] = -chv;
            //        CPsi_hzx_xn.p[i][k][j] = Chzey.p[i][k][j] * dx;
            //        Chzey.p[i][k][j] = Chzey.p[i][k][j] / kappa_mx;
            //    }
            //}
        }
    }
}

template<class T>
void cpml<T>::initCoefficientArraysXP(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dx,
        data3d<T>&Ceyhz, data3d<T>&Cezhy, data3d<T>&Chyez, data3d<T>&Chzey) {
    if (is_cpml_xp) {
        T sigmaOpt = sigmaMax * (pmlOrder + 1) / (150 * M_PI * dx * sqrt(epsR));
        unsigned iex = Ceyhz.nx - n_cpml_xp - 1;
        unsigned ihx = Chyez.nx - n_cpml_xp;
        for (unsigned i = 0; i < n_cpml_xp; i++) {
            T rho_e = (i + 0.25) / (T) n_cpml_xp;
            //T rho_e = (i) /(T)n_cpml_xp;
            T rho_m = (i + 0.75) / (T) n_cpml_xp;
            //T rho_m = (i) /(T)n_cpml_xp;
            T rho_e_pmlOrder = pow(fabs(rho_e), pmlOrder);
            T rho_m_pmlOrder = pow(fabs(rho_m), pmlOrder);
            T sigma_pex = sigmaOpt*rho_e_pmlOrder;
            T sigma_pmx = sigmaOpt * rho_m_pmlOrder;
            T kappa_ex = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            T kappa_mx = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            T alpha_ex = alphaMax * fabs(rho_e);
            T alpha_mx = alphaMax * fabs(rho_m);
            cpml_b_ex_xp.p[i] = exp((-dt / eps_0) * (sigma_pex / kappa_ex + alpha_ex));
            cpml_b_mx_xp.p[i] = exp((-dt / eps_0) * (sigma_pmx / kappa_mx + alpha_mx));
            if (rho_e != 0) {
                cpml_a_ex_xp.p[i] = 1 / dx * (cpml_b_ex_xp.p[i] - 1.0) * sigma_pex / (kappa_ex * (sigma_pex + kappa_ex * alpha_ex));
            } else {
                cpml_a_ex_xp.p[i] = 0;
            }
            if (rho_m != 0) {
                cpml_a_mx_xp.p[i] = 1 / dx * (cpml_b_mx_xp.p[i] - 1.0) * sigma_pmx / (kappa_mx * (sigma_pmx + kappa_mx * alpha_mx));
            } else {
                cpml_a_mx_xp.p[i] = 0;
            }

            for (unsigned j = 0; j < Psi_eyx_xp.ny; j++) {
                for (unsigned k = 0; k < Psi_eyx_xp.nz; k++) {
                    CPsi_eyx_xp.p[i][j][k] = Ceyhz.p[iex][j][k] * dx;
                    Ceyhz.p[iex][j][k] = Ceyhz.p[iex][j][k] / kappa_ex;
                }
            }
            for (unsigned j = 0; j < Psi_ezx_xp.nz; j++) {
                for (unsigned k = 0; k < Psi_ezx_xp.ny; k++) {
                    CPsi_ezx_xp.p[i][k][j] = Cezhy.p[iex][k][j] * dx;
                    Cezhy.p[iex][k][j] = Cezhy.p[iex][k][j] / kappa_ex;
                }
            }
            for (unsigned j = 0; j < Psi_hyx_xp.ny; j++) {
                for (unsigned k = 0; k < Psi_hyx_xp.nz; k++) {
                    CPsi_hyx_xp.p[i][j][k] = Chyez.p[ihx][j][k] * dx;
                    Chyez.p[ihx][j][k] = Chyez.p[ihx][j][k] / kappa_mx;
                }
            }
            for (unsigned j = 0; j < Psi_hzx_xp.nz; j++) {
                for (unsigned k = 0; k < Psi_hzx_xp.ny; k++) {
                    CPsi_hzx_xp.p[i][k][j] = Chzey.p[ihx][k][j] * dx;
                    Chzey.p[ihx][k][j] = Chzey.p[ihx][k][j] / kappa_mx;
                }
            }
            //T cev = 1 / (eps_0 / dt + sigma_pex / 2) / dx;
            //T chv = 1 / (mu_0 / dt + sigma_pmx / 2) / dx;
            //for (unsigned j = 0; j < Psi_eyx_xp.ny; j++) {
            //    for (unsigned k = 0; k < Psi_eyx_xp.nz; k++) {
            //        Ceyhz.p[iex][j][k] = -cev;
            //        CPsi_eyx_xp.p[i][j][k] = Ceyhz.p[iex][j][k] * dx;
            //        Ceyhz.p[iex][j][k] = Ceyhz.p[iex][j][k] / kappa_ex;
            //    }
            //}
            //for (unsigned j = 0; j < Psi_ezx_xp.nz; j++) {
            //    for (unsigned k = 0; k < Psi_ezx_xp.ny; k++) {
            //        Cezhy.p[iex][k][j] = cev;
            //        CPsi_ezx_xp.p[i][k][j] = Cezhy.p[iex][k][j] * dx;
            //        Cezhy.p[iex][k][j] = Cezhy.p[iex][k][j] / kappa_ex;
            //    }
            //}
            //for (unsigned j = 0; j < Psi_hyx_xp.ny; j++) {
            //    for (unsigned k = 0; k < Psi_hyx_xp.nz; k++) {
            //        Chyez.p[ihx][j][k] = chv;
            //        CPsi_hyx_xp.p[i][j][k] = Chyez.p[ihx][j][k] * dx;
            //        Chyez.p[ihx][j][k] = Chyez.p[ihx][j][k] / kappa_mx;
            //    }
            //}
            //for (unsigned j = 0; j < Psi_hzx_xp.nz; j++) {
            //    for (unsigned k = 0; k < Psi_hzx_xp.ny; k++) {
            //        Chzey.p[ihx][k][j] = -chv;
            //        CPsi_hzx_xp.p[i][k][j] = Chzey.p[ihx][k][j] * dx;
            //        Chzey.p[ihx][k][j] = Chzey.p[ihx][k][j] / kappa_mx;
            //    }
            //}
            ihx++;
            iex++;
        }
    }
}

template<class T>
void cpml<T>::initCoefficientArraysYN(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dy,
        data3d<T>&Cexhz, data3d<T>&Cezhx, data3d<T>&Chxez, data3d<T>&Chzex) {
    if (is_cpml_yn) {
        T sigmaOpt = sigmaMax * (pmlOrder + 1) / (150 * M_PI * dy * sqrt(epsR));
        for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
            T rho_e = (n_cpml_yn - j - 0.75) / (T) n_cpml_yn;
            //T rho_e = (n_cpml_yn - j) /(T)n_cpml_yn;
            T rho_m = (n_cpml_yn - j - 0.25) / (T) n_cpml_yn;
            //T rho_m = (n_cpml_yn - j) /(T)n_cpml_yn;
            T rho_e_pmlOrder = pow(fabs(rho_e), pmlOrder);
            T rho_m_pmlOrder = pow(fabs(rho_m), pmlOrder);
            T sigma_pey = sigmaOpt*rho_e_pmlOrder;
            T sigma_pmy = sigmaOpt * rho_m_pmlOrder;
            T kappa_ey = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            T kappa_my = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            T alpha_ey = alphaMax * fabs(rho_e);
            T alpha_my = alphaMax * fabs(rho_m);
            cpml_b_ey_yn.p[j] = exp((-dt / eps_0) * (sigma_pey / kappa_ey + alpha_ey));
            cpml_b_my_yn.p[j] = exp((-dt / eps_0) * (sigma_pmy / kappa_my + alpha_my));
            if (rho_e != 0) {
                cpml_a_ey_yn.p[j] = 1 / dy * (cpml_b_ey_yn.p[j] - 1.0) * sigma_pey / (kappa_ey * (sigma_pey + kappa_ey * alpha_ey));
            } else {
                cpml_a_ey_yn.p[j] = 0;
            }
            if (rho_m != 0) {
                cpml_a_my_yn.p[j] = 1 / dy * (cpml_b_my_yn.p[j] - 1.0) * sigma_pmy / (kappa_my * (sigma_pmy + kappa_my * alpha_my));
            } else {
                cpml_a_my_yn.p[j] = 0;
            }

            for (unsigned i = 0; i < Psi_exy_yn.nx; i++) {
                for (unsigned k = 0; k < Psi_exy_yn.nz; k++) {
                    CPsi_exy_yn.p[i][j][k] = Cexhz.p[i][jplus][k] * dy;
                    Cexhz.p[i][jplus][k] = Cexhz.p[i][jplus][k] / kappa_ey;
                }
            }
            for (unsigned i = 0; i < Psi_ezy_yn.nz; i++) {
                for (unsigned k = 0; k < Psi_ezy_yn.nx; k++) {
                    CPsi_ezy_yn.p[k][j][i] = Cezhx.p[k][jplus][i] * dy;
                    Cezhx.p[k][jplus][i] = Cezhx.p[k][jplus][i] / kappa_ey;
                }
            }
            for (unsigned i = 0; i < Psi_hxy_yn.nx; i++) {
                for (unsigned k = 0; k < Psi_hxy_yn.nz; k++) {
                    CPsi_hxy_yn.p[i][j][k] = Chxez.p[i][j][k] * dy;
                    Chxez.p[i][j][k] = Chxez.p[i][j][k] / kappa_my;
                }
            }
            for (unsigned i = 0; i < Psi_hzy_yn.nz; i++) {
                for (unsigned k = 0; k < Psi_hzy_yn.nx; k++) {
                    CPsi_hzy_yn.p[k][j][i] = Chzex.p[k][j][i] * dy;
                    Chzex.p[k][j][i] = Chzex.p[k][j][i] / kappa_my;
                }
            }
            //T cev = 1 / (eps_0 / dt + sigma_pey / 2) / dy;
            //T chv = 1 / (mu_0 / dt + sigma_pmy / 2) / dy;
            //for (unsigned i = 0; i < Psi_exy_yn.nx; i++) {
            //    for (unsigned k = 0; k < Psi_exy_yn.nz; k++) {
            //        Cexhz.p[i][jplus][k] = cev;
            //        CPsi_exy_yn.p[i][j][k] = Cexhz.p[i][jplus][k] * dy;
            //        Cexhz.p[i][jplus][k] = Cexhz.p[i][jplus][k] / kappa_ey;
            //    }
            //}
            //for (unsigned i = 0; i < Psi_ezy_yn.nz; i++) {
            //    for (unsigned k = 0; k < Psi_ezy_yn.nx; k++) {
            //        Cezhx.p[k][jplus][i] = -cev;
            //        CPsi_ezy_yn.p[k][j][i] = Cezhx.p[k][jplus][i] * dy;
            //        Cezhx.p[k][jplus][i] = Cezhx.p[k][jplus][i] / kappa_ey;
            //    }
            //}
            //for (unsigned i = 0; i < Psi_hxy_yn.nx; i++) {
            //    for (unsigned k = 0; k < Psi_hxy_yn.nz; k++) {
            //        Chxez.p[i][j][k] = -chv;
            //        CPsi_hxy_yn.p[i][j][k] = Chxez.p[i][j][k] * dy;
            //        Chxez.p[i][j][k] = Chxez.p[i][j][k] / kappa_my;
            //    }
            //}
            //for (unsigned i = 0; i < Psi_hzy_yn.nz; i++) {
            //    for (unsigned k = 0; k < Psi_hzy_yn.nx; k++) {
            //        Chzex.p[k][j][i] = chv;
            //        CPsi_hzy_yn.p[k][j][i] = Chzex.p[k][j][i] * dy;
            //        Chzex.p[k][j][i] = Chzex.p[k][j][i] / kappa_my;
            //    }
            //}
        }
    }
}

template<class T>
void cpml<T>::initCoefficientArraysYP(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dy,
        data3d<T>&Cexhz, data3d<T>&Cezhx, data3d<T>&Chxez, data3d<T>&Chzex) {
    if (is_cpml_yp) {
        T sigmaOpt = sigmaMax * (pmlOrder + 1) / (150 * M_PI * dy * sqrt(epsR));
        unsigned iex = Cexhz.ny - n_cpml_yp - 1;
        unsigned ihx = Chxez.ny - n_cpml_yp;
        for (unsigned j = 0; j < n_cpml_yp; j++) {
            T rho_e = (j + 0.25) / (T) n_cpml_yp;
            //T rho_e = (j) /(T)n_cpml_yp;
            T rho_m = (j + 0.75) / (T) n_cpml_yp;
            //T rho_m = (j) /(T)n_cpml_yp;
            T rho_e_pmlOrder = pow(fabs(rho_e), pmlOrder);
            T rho_m_pmlOrder = pow(fabs(rho_m), pmlOrder);
            T sigma_pey = sigmaOpt*rho_e_pmlOrder;
            T sigma_pmy = sigmaOpt * rho_m_pmlOrder;
            T kappa_ey = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            T kappa_my = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            T alpha_ey = alphaMax * fabs(rho_e);
            T alpha_my = alphaMax * fabs(rho_m);
            cpml_b_ey_yp.p[j] = exp((-dt / eps_0) * (sigma_pey / kappa_ey + alpha_ey));
            cpml_b_my_yp.p[j] = exp((-dt / eps_0) * (sigma_pmy / kappa_my + alpha_my));
            if (rho_e != 0) {
                cpml_a_ey_yp.p[j] = 1 / dy * (cpml_b_ey_yp.p[j] - 1.0) * sigma_pey / (kappa_ey * (sigma_pey + kappa_ey * alpha_ey));
            } else {
                cpml_a_my_yp.p[j] = 0;
            }
            if (rho_m != 0) {
                cpml_a_my_yp.p[j] = 1 / dy * (cpml_b_my_yp.p[j] - 1.0) * sigma_pmy / (kappa_my * (sigma_pmy + kappa_my * alpha_my));
            } else {
                cpml_a_my_yp.p[j] = 0;
            }

            for (unsigned i = 0; i < Psi_exy_yp.nx; i++) {
                for (unsigned k = 0; k < Psi_exy_yp.nz; k++) {
                    CPsi_exy_yp.p[i][j][k] = Cexhz.p[i][iex][k] * dy;
                    Cexhz.p[i][iex][k] = Cexhz.p[i][iex][k] / kappa_ey;
                }
            }
            for (unsigned i = 0; i < Psi_ezy_yp.nz; i++) {
                for (unsigned k = 0; k < Psi_ezy_yp.nx; k++) {
                    CPsi_ezy_yp.p[k][j][i] = Cezhx.p[k][iex][i] * dy;
                    Cezhx.p[k][iex][i] = Cezhx.p[k][iex][i] / kappa_ey;
                }
            }
            for (unsigned i = 0; i < Psi_hxy_yp.nx; i++) {
                for (unsigned k = 0; k < Psi_hxy_yp.nz; k++) {
                    CPsi_hxy_yp.p[i][j][k] = Chxez.p[i][ihx][k] * dy;
                    Chxez.p[i][ihx][k] = Chxez.p[i][ihx][k] / kappa_my;
                }
            }
            for (unsigned i = 0; i < Psi_hzy_yp.nz; i++) {
                for (unsigned k = 0; k < Psi_hzy_yp.nx; k++) {
                    CPsi_hzy_yp.p[k][j][i] = Chzex.p[k][ihx][i] * dy;
                    Chzex.p[k][ihx][i] = Chzex.p[k][ihx][i] / kappa_my;
                }
            }
            //T cev = 1 / (eps_0 / dt + sigma_pey / 2) / dy;
            //T chv = 1 / (mu_0 / dt + sigma_pmy / 2) / dy;
            //for (unsigned i = 0; i < Psi_exy_yp.nx; i++) {
            //    for (unsigned k = 0; k < Psi_exy_yp.nz; k++) {
            //        Cexhz.p[i][iex][k] = cev;
            //        CPsi_exy_yp.p[i][j][k] = Cexhz.p[i][iex][k] * dy;
            //        Cexhz.p[i][iex][k] = Cexhz.p[i][iex][k] / kappa_ey;
            //    }
            //}
            //for (unsigned i = 0; i < Psi_ezy_yp.nz; i++) {
            //    for (unsigned k = 0; k < Psi_ezy_yp.nx; k++) {
            //        Cezhx.p[k][iex][i] = -cev;
            //        CPsi_ezy_yp.p[k][j][i] = Cezhx.p[k][iex][i] * dy;
            //        Cezhx.p[k][iex][i] = Cezhx.p[k][iex][i] / kappa_ey;
            //    }
            //}
            //for (unsigned i = 0; i < Psi_hxy_yp.nx; i++) {
            //    for (unsigned k = 0; k < Psi_hxy_yp.nz; k++) {
            //        Chxez.p[i][ihx][k] = -chv;
            //        CPsi_hxy_yp.p[i][j][k] = Chxez.p[i][ihx][k] * dy;
            //        Chxez.p[i][ihx][k] = Chxez.p[i][ihx][k] / kappa_my;
            //    }
            //}
            //for (unsigned i = 0; i < Psi_hzy_yp.nz; i++) {
            //    for (unsigned k = 0; k < Psi_hzy_yp.nx; k++) {
            //        Chzex.p[k][ihx][i] = chv;
            //        CPsi_hzy_yp.p[k][j][i] = Chzex.p[k][ihx][i] * dy;
            //        Chzex.p[k][ihx][i] = Chzex.p[k][ihx][i] / kappa_my;
            //    }
            //}
            ihx++;
            iex++;
        }
    }
}

template<class T>
void cpml<T>::initCoefficientArraysZN(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dz,
        data3d<T>&Ceyhx, data3d<T>&Cexhy, data3d<T>&Chyex, data3d<T>&Chxey) {
    if (is_cpml_zn) {
        T sigmaOpt = sigmaMax * (pmlOrder + 1) / (150 * M_PI * dz * sqrt(epsR));
        for (unsigned k = 0, iplus = 1; k < n_cpml_zn; k++, iplus++) {
            T rho_e = (n_cpml_zn - k - 0.75) / (T) n_cpml_zn;
            //T rho_e = (n_cpml_zn - k) /(T)n_cpml_zn;
            T rho_m = (n_cpml_zn - k - 0.25) / (T) n_cpml_zn;
            //T rho_m = (n_cpml_zn - k) /(T)n_cpml_zn;
            T rho_e_pmlOrder = pow(fabs(rho_e), pmlOrder);
            T rho_m_pmlOrder = pow(fabs(rho_m), pmlOrder);
            T sigma_pez = sigmaOpt*rho_e_pmlOrder;
            T sigma_pmz = sigmaOpt * rho_m_pmlOrder;
            T kappa_ez = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            T kappa_mz = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            T alpha_ez = alphaMax * fabs(rho_e);
            T alpha_mz = alphaMax * fabs(rho_m);
            cpml_b_ez_zn.p[k] = exp((-dt / eps_0) * (sigma_pez / kappa_ez + alpha_ez));
            cpml_b_mz_zn.p[k] = exp((-dt / eps_0) * (sigma_pmz / kappa_mz + alpha_mz));
            if (rho_e != 0) {
                cpml_a_ez_zn.p[k] = 1 / dz * (cpml_b_ez_zn.p[k] - 1.0) * sigma_pez / (kappa_ez * (sigma_pez + kappa_ez * alpha_ez));
            } else {
                cpml_a_ez_zn.p[k] = 0;
            }
            if (rho_m != 0) {
                cpml_a_mz_zn.p[k] = 1 / dz * (cpml_b_mz_zn.p[k] - 1.0) * sigma_pmz / (kappa_mz * (sigma_pmz + kappa_mz * alpha_mz));
            } else {
                cpml_a_mz_zn.p[k] = 0;
            }

            for (unsigned j = 0; j < Psi_exz_zn.ny; j++) {
                for (unsigned i = 0; i < Psi_exz_zn.nx; i++) {
                    CPsi_exz_zn.p[i][j][k] = Cexhy.p[i][j][iplus] * dz;
                    Cexhy.p[i][j][iplus] = Cexhy.p[i][j][iplus] / kappa_ez;
                }
            }
            for (unsigned j = 0; j < Psi_eyz_zn.nx; j++) {
                for (unsigned i = 0; i < Psi_eyz_zn.ny; i++) {
                    CPsi_eyz_zn.p[j][i][k] = Ceyhx.p[j][i][iplus] * dz;
                    Ceyhx.p[j][i][iplus] = Ceyhx.p[j][i][iplus] / kappa_ez;
                }
            }
            for (unsigned j = 0; j < Psi_hxz_zn.ny; j++) {
                for (unsigned i = 0; i < Psi_hxz_zn.nx; i++) {
                    CPsi_hxz_zn.p[i][j][k] = Chxey.p[i][j][k] * dz;
                    Chxey.p[i][j][k] = Chxey.p[i][j][k] / kappa_mz;
                }
            }
            for (unsigned j = 0; j < Psi_hyz_zn.nx; j++) {
                for (unsigned i = 0; i < Psi_hyz_zn.ny; i++) {
                    CPsi_hyz_zn.p[j][i][k] = Chyex.p[j][i][k] * dz;
                    Chyex.p[j][i][k] = Chyex.p[j][i][k] / kappa_mz;
                }
            }
            //T cev = 1 / (eps_0 / dt + sigma_pez / 2) / dz;
            //T chv = 1 / (mu_0 / dt + sigma_pmz / 2) / dz;
            //for (unsigned j = 0; j < Psi_exz_zn.ny; j++) {
            //    for (unsigned i = 0; i < Psi_exz_zn.nx; i++) {
            //        Cexhy.p[i][j][iplus] = -cev;
            //        CPsi_exz_zn.p[i][j][k] = Cexhy.p[i][j][iplus] * dz;
            //        Cexhy.p[i][j][iplus] = Cexhy.p[i][j][iplus] / kappa_ez;
            //    }
            //}
            //for (unsigned j = 0; j < Psi_eyz_zn.nx; j++) {
            //    for (unsigned i = 0; i < Psi_eyz_zn.ny; i++) {
            //        Ceyhx.p[j][i][iplus] = cev;
            //        CPsi_eyz_zn.p[j][i][k] = Ceyhx.p[j][i][iplus] * dz;
            //        Ceyhx.p[j][i][iplus] = Ceyhx.p[j][i][iplus] / kappa_ez;
            //    }
            //}
            //for (unsigned j = 0; j < Psi_hxz_zn.ny; j++) {
            //    for (unsigned i = 0; i < Psi_hxz_zn.nx; i++) {
            //        Chxey.p[i][j][k] = chv;
            //        CPsi_hxz_zn.p[i][j][k] = Chxey.p[i][j][k] * dz;
            //        Chxey.p[i][j][k] = Chxey.p[i][j][k] / kappa_mz;
            //    }
            //}
            //for (unsigned j = 0; j < Psi_hyz_zn.nx; j++) {
            //    for (unsigned i = 0; i < Psi_hyz_zn.ny; i++) {
            //        Chyex.p[j][i][k] = -chv;
            //        CPsi_hyz_zn.p[j][i][k] = Chyex.p[j][i][k] * dz;
            //        Chyex.p[j][i][k] = Chyex.p[j][i][k] / kappa_mz;
            //    }
            //}
        }
    }
}

template<class T>
void cpml<T>::initCoefficientArraysZP(short pmlOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dz,
        data3d<T>&Ceyhx, data3d<T>&Cexhy, data3d<T>&Chyex, data3d<T>&Chxey) {
    if (is_cpml_zp) {
        T sigmaOpt = sigmaMax * (pmlOrder + 1) / (150 * M_PI * dz * sqrt(epsR));
        unsigned iez = Ceyhx.nz - n_cpml_zp - 1;
        unsigned ihz = Chyex.nz - n_cpml_zp;
        for (unsigned k = 0; k < n_cpml_zp; k++) {
            T rho_e = (k + 0.25) / (T) n_cpml_zp;
            //T rho_e = (k) /(T)n_cpml_zp;
            T rho_m = (k + 0.75) / (T) n_cpml_zp;
            //T rho_m = (k) /(T)n_cpml_zp;
            T rho_e_pmlOrder = pow(fabs(rho_e), pmlOrder);
            T rho_m_pmlOrder = pow(fabs(rho_m), pmlOrder);
            T sigma_pez = sigmaOpt*rho_e_pmlOrder;
            T sigma_pmz = sigmaOpt * rho_m_pmlOrder;
            T kappa_ez = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            T kappa_mz = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            T alpha_ez = alphaMax * fabs(rho_e);
            T alpha_mz = alphaMax * fabs(rho_m);

            cpml_b_ez_zp.p[k] = exp((-dt / eps_0) * (sigma_pez / kappa_ez + alpha_ez));
            cpml_b_mz_zp.p[k] = exp((-dt / eps_0) * (sigma_pmz / kappa_mz + alpha_mz));
            if (rho_e != 0) {
                cpml_a_ez_zp.p[k] = 1 / dz * (cpml_b_ez_zp.p[k] - 1.0) * sigma_pez / (kappa_ez * (sigma_pez + kappa_ez * alpha_ez));
            } else {
                cpml_a_ez_zp.p[k] = 0;
            }
            if (rho_m != 0) {
                cpml_a_mz_zp.p[k] = 1 / dz * (cpml_b_mz_zp.p[k] - 1.0) * sigma_pmz / (kappa_mz * (sigma_pmz + kappa_mz * alpha_mz));
            } else {
                cpml_a_mz_zp.p[k] = 0;
            }
            for (unsigned j = 0; j < Psi_eyz_zp.ny; j++) {
                for (unsigned i = 0; i < Psi_eyz_zp.nx; i++) {
                    CPsi_eyz_zp.p[i][j][k] = Ceyhx.p[i][j][iez] * dz;
                    Ceyhx.p[i][j][iez] = Ceyhx.p[i][j][iez] / kappa_ez;
                }
            }
            for (unsigned j = 0; j < Psi_exz_zp.nx; j++) {
                for (unsigned i = 0; i < Psi_exz_zp.ny; i++) {
                    CPsi_exz_zp.p[j][i][k] = Cexhy.p[j][i][iez] * dz;
                    Cexhy.p[j][i][iez] = Cexhy.p[j][i][iez] / kappa_ez;
                }
            }
            for (unsigned j = 0; j < Psi_hyz_zp.ny; j++) {
                for (unsigned i = 0; i < Psi_hyz_zp.nx; i++) {
                    CPsi_hyz_zp.p[i][j][k] = Chyex.p[i][j][ihz] * dz;
                    Chyex.p[i][j][ihz] = Chyex.p[i][j][ihz] / kappa_mz;
                }
            }
            for (unsigned j = 0; j < Psi_hxz_zp.nx; j++) {
                for (unsigned i = 0; i < Psi_hxz_zp.ny; i++) {
                    CPsi_hxz_zp.p[j][i][k] = Chxey.p[j][i][ihz] * dz;
                    Chxey.p[j][i][ihz] = Chxey.p[j][i][ihz] / kappa_mz;
                }
            }
            //T cev = 1 / (eps_0 / dt + sigma_pez / 2) / dz;
            //T chv = 1 / (mu_0 / dt + sigma_pmz / 2) / dz;
            //for (unsigned j = 0; j < Psi_eyz_zp.ny; j++) {
            //    for (unsigned i = 0; i < Psi_eyz_zp.nx; i++) {
            //        Ceyhx.p[i][j][iez] = cev;
            //        CPsi_eyz_zp.p[i][j][k] = Ceyhx.p[i][j][iez] * dz;
            //        Ceyhx.p[i][j][iez] = Ceyhx.p[i][j][iez] / kappa_ez;
            //    }
            //}
            //for (unsigned j = 0; j < Psi_exz_zp.nx; j++) {
            //    for (unsigned i = 0; i < Psi_exz_zp.ny; i++) {
            //        Cexhy.p[j][i][iez] = -cev;
            //        CPsi_exz_zp.p[j][i][k] = Cexhy.p[j][i][iez] * dz;
            //        Cexhy.p[j][i][iez] = Cexhy.p[j][i][iez] / kappa_ez;
            //    }
            //}
            //for (unsigned j = 0; j < Psi_hyz_zp.ny; j++) {
            //    for (unsigned i = 0; i < Psi_hyz_zp.nx; i++) {
            //        Chyex.p[i][j][ihz] = -chv;
            //        CPsi_hyz_zp.p[i][j][k] = Chyex.p[i][j][ihz] * dz;
            //        Chyex.p[i][j][ihz] = Chyex.p[i][j][ihz] / kappa_mz;
            //    }
            //}
            //for (unsigned j = 0; j < Psi_hxz_zp.nx; j++) {
            //    for (unsigned i = 0; i < Psi_hxz_zp.ny; i++) {
            //        Chxey.p[j][i][ihz] = chv;
            //        CPsi_hxz_zp.p[j][i][k] = Chxey.p[j][i][ihz] * dz;
            //        Chxey.p[j][i][ihz] = Chxey.p[j][i][ihz] / kappa_mz;
            //    }
            //}
            ihz++;
            iez++;
        }
    }
}

template<class T>
void cpml<T>::updateCPML_E_Fields(data3d<T>& Ex, data3d<T>& Ey, data3d<T>& Ez,
        const data3d<T>& Hx, const data3d<T>& Hy, const data3d<T>& Hz) {
    updatePsiForEFields(Hx, Hy, Hz);
    updateEFieldCPML_x(Ey, Ez);
    updateEFieldCPML_y(Ex, Ez);
    updateEFieldCPML_z(Ex, Ey);
}

template<class T>
void cpml<T>::updateCPML_M_Fields(data3d<T>& Hx, data3d<T>& Hy, data3d<T>& Hz,
        const data3d<T>& Ex, const data3d<T>& Ey, const data3d<T>& Ez) {
    updatePsiForMFields(Ex, Ey, Ez);
    updateMFieldCPML_x(Hy, Hz);
    updateMFieldCPML_y(Hx, Hz);
    updateMFieldCPML_z(Hx, Hy);
}

template<class T>
void cpml<T>::updatePsiForEFields(const data3d<T>& Hx, const data3d<T>& Hy, const data3d<T>& Hz) {
    if (is_cpml_xn) {
        updatePsi_eyz_zn(Hx);
        updatePsi_ezy_yn(Hx);
    }
    if (is_cpml_xp) {
        updatePsi_eyz_zp(Hx);
        updatePsi_ezy_yp(Hx);
    }
    if (is_cpml_yn) {
        updatePsi_exz_zn(Hy);
        updatePsi_ezx_xn(Hy);
    }
    if (is_cpml_yp) {
        updatePsi_exz_zp(Hy);
        updatePsi_ezx_xp(Hy);
    }
    if (is_cpml_zn) {
        updatePsi_eyx_xn(Hz);
        updatePsi_exy_yn(Hz);
    }
    if (is_cpml_zp) {
        updatePsi_eyx_xp(Hz);
        updatePsi_exy_yp(Hz);
    }
}

template<class T>
void cpml<T>::updatePsiForMFields(const data3d<T>& Ex, const data3d<T>& Ey, const data3d<T>& Ez) {
    if (is_cpml_xn) {
        updatePsi_hyz_zn(Ex);
        updatePsi_hzy_yn(Ex);
    }
    if (is_cpml_xp) {
        updatePsi_hyz_zp(Ex);
        updatePsi_hzy_yp(Ex);
    }
    if (is_cpml_yn) {
        updatePsi_hxz_zn(Ey);
        updatePsi_hzx_xn(Ey);
    }
    if (is_cpml_yp) {
        updatePsi_hxz_zp(Ey);
        updatePsi_hzx_xp(Ey);
    }
    if (is_cpml_zn) {
        updatePsi_hyx_xn(Ez);
        updatePsi_hxy_yn(Ez);
    }
    if (is_cpml_zp) {
        updatePsi_hyx_xp(Ez);
        updatePsi_hxy_yp(Ez);
    }
}

template<class T>
void cpml<T>::updateEFieldCPML_x(data3d<T>& Ey, data3d<T>& Ez) {
    // x negetive region
    if (is_cpml_xn) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Ey.ny; j++) {
            for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
                for (unsigned k = 0; k < Ey.nz; k++) {
                    Ey.p[iplus][j][k] += CPsi_eyx_xn.p[i][j][k] * Psi_eyx_xn.p[i][j][k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Ez.ny; j++) {
            for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ez.p[iplus][j][k] += CPsi_ezx_xn.p[i][j][k] * Psi_ezx_xn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_xp) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Ey.ny; j++) {
            for (unsigned i = 0, ihy = Ey.nx - n_cpml_xp - 1; i < n_cpml_xp; i++, ihy++) {
                for (unsigned k = 0; k < Ey.nz; k++) {
                    Ey.p[ihy][j][k] += CPsi_eyx_xp.p[i][j][k] * Psi_eyx_xp.p[i][j][k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Ez.ny; j++) {
            for (unsigned i = 0, ihz = Ez.nx - n_cpml_xp - 1; i < n_cpml_xp; i++, ihz++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ez.p[ihz][j][k] += CPsi_ezx_xp.p[i][j][k] * Psi_ezx_xp.p[i][j][k];
                }
            }
        }
    }
}

template<class T>
void cpml<T>::updateEFieldCPML_y(data3d<T>& Ex, data3d<T>& Ez) {
    // x negetive region
    if (is_cpml_yn) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned i = 0; i < Ex.nx; i++) {
            for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
                for (unsigned k = 0; k < Ex.nz; k++) {
                    Ex.p[i][jplus][k] += CPsi_exy_yn.p[i][j][k] * Psi_exy_yn.p[i][j][k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned i = 0; i < Ez.nx; i++) {
            for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ez.p[i][jplus][k] += CPsi_ezy_yn.p[i][j][k] * Psi_ezy_yn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_yp) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned i = 0; i < Ex.nx; i++) {
            for (unsigned j = 0, jhx = Ex.ny - n_cpml_yp - 1; j < n_cpml_yp; j++, jhx++) {
                for (unsigned k = 0; k < Ex.nz; k++) {
                    Ex.p[i][jhx][k] += CPsi_exy_yp.p[i][j][k] * Psi_exy_yp.p[i][j][k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned i = 0; i < Ez.nx; i++) {
            for (unsigned j = 0, jhz = Ez.ny - n_cpml_yp - 1; j < n_cpml_yp; j++, jhz++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ez.p[i][jhz][k] += CPsi_ezy_yp.p[i][j][k] * Psi_ezy_yp.p[i][j][k];
                }
            }
        }
    }
}

template<class T>
void cpml<T>::updateEFieldCPML_z(data3d<T>& Ex, data3d<T>& Ey) {
    // x negetive region
    if (is_cpml_zn) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Ey.ny; j++) {
            for (unsigned k = 0, kplus = 1; k < n_cpml_zn; k++, kplus++) {

                for (unsigned i = 0; i < Ey.nx; i++) {
                    Ey.p[i][j][kplus] += CPsi_eyz_zn.p[i][j][k] * Psi_eyz_zn.p[i][j][k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Ex.ny; j++) {
            for (unsigned k = 0, kplus = 1; k < n_cpml_zn; k++, kplus++) {
                for (unsigned i = 0; i < Ex.nx; i++) {
                    Ex.p[i][j][kplus] += CPsi_exz_zn.p[i][j][k] * Psi_exz_zn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_zp) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Ey.ny; j++) {
            for (unsigned k = 0, khy = Ey.nz - n_cpml_zp - 1; k < n_cpml_zp; k++, khy++) {
                for (unsigned i = 0; i < Ey.nx; i++) {
                    Ey.p[i][j][khy] += CPsi_eyz_zp.p[i][j][k] * Psi_eyz_zp.p[i][j][k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Ex.ny; j++) {
            for (unsigned k = 0, khx = Ex.nz - n_cpml_zp - 1; k < n_cpml_zp; k++, khx++) {
                for (unsigned i = 0; i < Ex.nx; i++) {
                    Ex.p[i][j][khx] += CPsi_exz_zp.p[i][j][k] * Psi_exz_zp.p[i][j][k];
                }
            }
        }
    }
}

template<class T>
void cpml<T>::updateMFieldCPML_x(data3d<T>& Hy, data3d<T>& Hz) {
    // x negetive region
    if (is_cpml_xn) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Hy.ny; j++) {
            for (unsigned i = 0; i < n_cpml_xn; i++) {
                for (unsigned k = 0; k < Hy.nz; k++) {
                    Hy.p[i][j][k] += CPsi_hyx_xn.p[i][j][k] * Psi_hyx_xn.p[i][j][k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Hz.ny; j++) {
            for (unsigned i = 0; i < n_cpml_xn; i++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hz.p[i][j][k] += CPsi_hzx_xn.p[i][j][k] * Psi_hzx_xn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_xp) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Hy.ny; j++) {
            for (unsigned i = 0, ihy = Hy.nx - n_cpml_xp; i < n_cpml_xp; i++, ihy++) {
                for (unsigned k = 0; k < Hy.nz; k++) {
                    Hy.p[ihy][j][k] += CPsi_hyx_xp.p[i][j][k] * Psi_hyx_xp.p[i][j][k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Hz.ny; j++) {
            for (unsigned i = 0, ihz = Hz.nx - n_cpml_xp; i < n_cpml_xp; i++, ihz++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hz.p[ihz][j][k] += CPsi_hzx_xp.p[i][j][k] * Psi_hzx_xp.p[i][j][k];
                }
            }
        }
    }
}

template<class T>
void cpml<T>::updateMFieldCPML_y(data3d<T>& Hx, data3d<T>& Hz) {
    // x negetive region
    if (is_cpml_yn) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned i = 0; i < Hx.nx; i++) {
            for (unsigned j = 0; j < n_cpml_yn; j++) {
                for (unsigned k = 0; k < Hx.nz; k++) {
                    Hx.p[i][j][k] += CPsi_hxy_yn.p[i][j][k] * Psi_hxy_yn.p[i][j][k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned i = 0; i < Hz.nx; i++) {
            for (unsigned j = 0; j < n_cpml_yn; j++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hz.p[i][j][k] += CPsi_hzy_yn.p[i][j][k] * Psi_hzy_yn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_yp) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned i = 0; i < Hx.nx; i++) {
            for (unsigned j = 0, ihx = Hx.ny - n_cpml_yp; j < n_cpml_yp; j++, ihx++) {
                for (unsigned k = 0; k < Hx.nz; k++) {
                    Hx.p[i][ihx][k] += CPsi_hxy_yp.p[i][j][k] * Psi_hxy_yp.p[i][j][k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned i = 0; i < Hz.nx; i++) {
            for (unsigned j = 0, ihz = Hz.ny - n_cpml_yp; j < n_cpml_yp; j++, ihz++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hz.p[i][ihz][k] += CPsi_hzy_yp.p[i][j][k] * Psi_hzy_yp.p[i][j][k];
                }
            }
        }
    }
}

template<class T>
void cpml<T>::updateMFieldCPML_z(data3d<T>& Hx, data3d<T>& Hy) {
    //TODO add update Hz in pml region
    // x negetive region
    if (is_cpml_zn) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Hy.ny; j++) {
            for (unsigned k = 0; k < n_cpml_zn; k++) {
                for (unsigned i = 0; i < Hy.nx; i++) {
                    Hy.p[i][j][k] += CPsi_hyz_zn.p[i][j][k] * Psi_hyz_zn.p[i][j][k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned j = 0; j < Hx.ny; j++) {
            for (unsigned k = 0; k < n_cpml_zn; k++) {
                for (unsigned i = 0; i < Hx.nx; i++) {
                    Hx.p[i][j][k] += CPsi_hxz_zn.p[i][j][k] * Psi_hxz_zn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_zp) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned i = 0; i < Hy.nx; i++) {
            for (unsigned k = 0, khy = Hy.nz - n_cpml_zp; k < n_cpml_zp; k++, khy++) {
                for (unsigned j = 0; j < Hy.ny; j++) {
                    Hy.p[i][j][khy] += CPsi_hyz_zp.p[i][j][k] * Psi_hyz_zp.p[i][j][k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
        for (unsigned i = 0; i < Hx.nx; i++) {
            for (unsigned k = 0, khx = Hx.nz - n_cpml_zp; k < n_cpml_zp; k++, khx++) {
                for (unsigned j = 0; j < Hx.ny; j++) {
                    Hx.p[i][j][khx] += CPsi_hxz_zp.p[i][j][k] * Psi_hxz_zp.p[i][j][k];
                }
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_eyz_zp(const data3d<T>& Hx) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_eyz_zp.nx; i++) {
        for (unsigned k = 0, ikz = Hx.nz - n_cpml_zp; k < n_cpml_zp; k++, ikz++) {
            for (unsigned j = 0; j < Psi_eyz_zp.ny; j++) {
                Psi_eyz_zp.p[i][j][k] = Psi_eyz_zp.p[i][j][k] * cpml_b_ez_zp.p[k] + cpml_a_ez_zp.p[k]*(Hx.p[i][j][ikz] - Hx.p[i][j][ikz - 1]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_ezy_yp(const data3d<T>& Hx) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_ezy_yp.nx; i++) {
        for (unsigned j = 0, jhy = Hx.ny - n_cpml_yp; j < n_cpml_yp; j++, jhy++) {
            for (unsigned k = 0; k < Psi_ezy_yp.nz; k++) {
                Psi_ezy_yp.p[i][j][k] = Psi_ezy_yp.p[i][j][k] * cpml_b_ey_yp.p[j] + cpml_a_ey_yp.p[j]*(Hx.p[i][jhy][k] - Hx.p[i][jhy - 1][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_exz_zp(const data3d<T>& Hy) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_exz_zp.nx; i++) {
        for (unsigned k = 0, iez = Hy.nz - n_cpml_zp; k < n_cpml_zp; k++, iez++) {
            for (unsigned j = 0; j < Psi_exz_zp.ny; j++) {
                Psi_exz_zp.p[i][j][k] = cpml_b_ez_zp.p[k] * Psi_exz_zp.p[i][j][k] + cpml_a_ez_zp.p[k]*(Hy.p[i][j][iez ] - Hy.p[i][j][iez - 1]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_ezx_xp(const data3d<T>& Hy) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned j = 0; j < Psi_ezx_xp.ny; j++) {
        for (unsigned i = 0, iex = Hy.nx - n_cpml_xp; i < n_cpml_xp; i++, iex++) {
            for (unsigned k = 0; k < Psi_ezx_xp.nz; k++) {
                Psi_ezx_xp.p[i][j][k] = cpml_b_ex_xp.p[i] * Psi_ezx_xp.p[i][j][k] + cpml_a_ex_xp.p[i]*(Hy.p[iex][j][k] - Hy.p[iex - 1][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_eyx_xp(const data3d<T>& Hz) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned j = 0; j < Psi_eyx_xp.ny; j++) {
        for (unsigned i = 0, ihz = Hz.nx - n_cpml_xp; i < n_cpml_xp; i++, ihz++) {
            for (unsigned k = 0; k < Psi_eyx_xp.nz; k++) {
                Psi_eyx_xp.p[i][j][k] = cpml_b_ex_xp.p[i] * Psi_eyx_xp.p[i][j][k] + cpml_a_ex_xp.p[i]*(Hz.p[ihz][j][k] - Hz.p[ihz - 1][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_exy_yp(const data3d<T>& Hz) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_exy_yp.nx; i++) {
        for (unsigned j = 0, ihy = Hz.ny - n_cpml_yp; j < n_cpml_yp; j++, ihy++) {
            for (unsigned k = 0; k < Psi_exy_yp.nz; k++) {
                Psi_exy_yp.p[i][j][k] = Psi_exy_yp.p[i][j][k] * cpml_b_ey_yp.p[j] + cpml_a_ey_yp.p[j]*(Hz.p[i][ihy][k] - Hz.p[i][ihy - 1][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_hyz_zp(const data3d<T>& Ex) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_hyz_zp.nx; i++) {
        for (unsigned k = 0, ikz = Ex.nz - n_cpml_zp; k < n_cpml_zp; k++, ikz++) {
            for (unsigned j = 0; j < Psi_hyz_zp.ny; j++) {
                Psi_hyz_zp.p[i][j][k] = Psi_hyz_zp.p[i][j][k] * cpml_b_mz_zp.p[k] + cpml_a_mz_zp.p[k]*(Ex.p[i][j][ikz] - Ex.p[i][j][ikz - 1]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_hzy_yp(const data3d<T>& Ex) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_hzy_yp.nx; i++) {
        for (unsigned j = 0, jhy = Ex.ny - n_cpml_yp; j < n_cpml_yp; j++, jhy++) {
            for (unsigned k = 0; k < Psi_hzy_yp.nz; k++) {
                Psi_hzy_yp.p[i][j][k] = Psi_hzy_yp.p[i][j][k] * cpml_b_my_yp.p[j] + cpml_a_my_yp.p[j]*(Ex.p[i][jhy][k] - Ex.p[i][jhy - 1 ][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_hxz_zp(const data3d<T>& Ey) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned j = 0; j < Psi_hxz_zp.ny; j++) {
        for (unsigned k = 0, khz = Ey.nz - n_cpml_zp; k < n_cpml_zp; k++, khz++) {
            for (unsigned i = 0; i < Psi_hxz_zp.nx; i++) {
                Psi_hxz_zp.p[i][j][k] = Psi_hxz_zp.p[i][j][k] * cpml_b_mz_zp.p[k] + cpml_a_mz_zp.p[k]*(Ey.p[i][j][khz] - Ey.p[i][j][khz - 1]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_hzx_xp(const data3d<T>& Ey) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned j = 0; j < Psi_hzx_xp.ny; j++) {
        for (unsigned i = 0, ihx = Ey.nx - n_cpml_xp; i < n_cpml_xp; i++, ihx++) {
            for (unsigned k = 0; k < Psi_hzx_xp.nz; k++) {
                Psi_hzx_xp.p[i][j][k] = Psi_hzx_xp.p[i][j][k] * cpml_b_mx_xp.p[i] + cpml_a_mx_xp.p[i]*(Ey.p[ihx][j][k] - Ey.p[ihx - 1][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_hyx_xp(const data3d<T>& Ez) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned j = 0; j < Psi_hyx_xp.ny; j++) {
        for (unsigned i = 0, ihx = Ez.nx - n_cpml_xp; i < n_cpml_xp; i++, ihx++) {
            for (unsigned k = 0; k < Psi_hyx_xp.nz; k++) {
                Psi_hyx_xp.p[i][j][k] = Psi_hyx_xp.p[i][j][k] * cpml_b_mx_xp.p[i] + cpml_a_mx_xp.p[i]*(Ez.p[ihx][j][k] - Ez.p[ihx - 1 ][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_hxy_yp(const data3d<T>& Ez) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_hxy_yp.nx; i++) {
        for (unsigned j = 0, jhy = Ez.ny - n_cpml_yp; j < n_cpml_yp; j++, jhy++) {
            for (unsigned k = 0; k < Psi_hxy_yp.nz; k++) {
                Psi_hxy_yp.p[i][j][k] = Psi_hxy_yp.p[i][j][k] * cpml_b_my_yp.p[j] + cpml_a_my_yp.p[j]*(Ez.p[i][jhy][k] - Ez.p[i][jhy - 1][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_eyz_zn(const data3d<T>& Hx) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_eyz_zn.nx; i++) {
        for (unsigned k = 0, kplus = 1; k < n_cpml_zn; k++, kplus++) {
            for (unsigned j = 0; j < Psi_eyz_zn.ny; j++) {
                Psi_eyz_zn.p[i][j][k] = cpml_b_ez_zn.p[k] * Psi_eyz_zn.p[i][j][k] + cpml_a_ez_zn.p[k]*(Hx.p[i][j][kplus] - Hx.p[i][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_ezy_yn(const data3d<T>& Hx) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_ezy_yn.nx; i++) {
        for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
            for (unsigned k = 0; k < Psi_ezy_yn.nz; k++) {
                Psi_ezy_yn.p[i][j][k] = cpml_b_ey_yn.p[j] * Psi_ezy_yn.p[i][j][k] + cpml_a_ey_yn.p[j]*(Hx.p[i][jplus][k] - Hx.p[i][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_exz_zn(const data3d<T>& Hy) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_exz_zn.nx; i++) {
        for (unsigned k = 0, kplus = 1; k < n_cpml_zn; k++, kplus++) {
            for (unsigned j = 0; j < Psi_exz_zn.ny; j++) {
                Psi_exz_zn.p[i][j][k] = cpml_b_ez_zn.p[k] * Psi_exz_zn.p[i][j][k] + cpml_a_ez_zn.p[k]*(Hy.p[i][j][kplus] - Hy.p[i][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_ezx_xn(const data3d<T>& Hy) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned j = 0; j < Psi_ezx_xn.ny; j++) {
        for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
            for (unsigned k = 0; k < Psi_ezx_xn.nz; k++) {
                Psi_ezx_xn.p[i][j][k] = cpml_b_ex_xn.p[i] * Psi_ezx_xn.p[i][j][k] + cpml_a_ex_xn.p[i]*(Hy.p[iplus][j][k] - Hy.p[i][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_eyx_xn(const data3d<T>& Hz) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned j = 0; j < Psi_eyx_xn.ny; j++) {
        for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
            for (unsigned k = 0; k < Psi_eyx_xn.nz; k++) {
                Psi_eyx_xn.p[i][j][k] = cpml_b_ex_xn.p[i] * Psi_eyx_xn.p[i][j][k] + cpml_a_ex_xn.p[i]*(Hz.p[iplus][j][k] - Hz.p[i][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_exy_yn(const data3d<T>& Hz) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned k = 0; k < Psi_exy_yn.nz; k++) {
        for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
            for (unsigned i = 0; i < Psi_exy_yn.ny; i++) {
                Psi_exy_yn.p[i][j][k] = cpml_b_ey_yn.p[j] * Psi_exy_yn.p[i][j][k] + cpml_a_ey_yn.p[j]*(Hz.p[i][jplus][k] - Hz.p[i][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_hyz_zn(const data3d<T>& Ex) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned j = 0; j < Psi_hyz_zn.ny; j++) {
        for (unsigned k = 0, kplus = 1; k < n_cpml_zn; k++, kplus++) {
            for (unsigned i = 0; i < Psi_hyz_zn.nx; i++) {
                Psi_hyz_zn.p[i][j][k] = cpml_b_mz_zn.p[k] * Psi_hyz_zn.p[i][j][k] + cpml_a_mz_zn.p[k]*(Ex.p[i][j][kplus] - Ex.p[i][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_hzy_yn(const data3d<T>& Ex) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_hzy_yn.nx; i++) {
        for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
            for (unsigned k = 0; k < Psi_hzy_yn.nz; k++) {
                Psi_hzy_yn.p[i][j][k] = cpml_b_my_yn.p[j] * Psi_hzy_yn.p[i][j][k] + cpml_a_my_yn.p[j]*(Ex.p[i][jplus][k] - Ex.p[i][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_hxz_zn(const data3d<T>& Ey) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_hxz_zn.nx; i++) {
        for (unsigned k = 0, kplus = 1; k < n_cpml_zn; k++, kplus++) {
            for (unsigned j = 0; j < Psi_hxz_zn.ny; j++) {
                Psi_hxz_zn.p[i][j][k] = cpml_b_mz_zn.p[k] * Psi_hxz_zn.p[i][j][k] + cpml_a_mz_zn.p[k]*(Ey.p[i][j][kplus] - Ey.p[i][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_hzx_xn(const data3d<T>& Ey) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned j = 0; j < Psi_hzx_xn.ny; j++) {
        for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
            for (unsigned k = 0; k < Psi_hzx_xn.nz; k++) {
                Psi_hzx_xn.p[i][j][k] = Psi_hzx_xn.p[i][j][k] * cpml_b_mx_xn.p[i] + cpml_a_mx_xn.p[i]*(Ey.p[iplus][j][k] - Ey.p[i][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_hyx_xn(const data3d<T>& Ez) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned j = 0; j < Psi_hyx_xn.ny; j++) {
        for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
            for (unsigned k = 0; k < Psi_hyx_xn.nz; k++) {
                Psi_hyx_xn.p[i][j][k] = Psi_hyx_xn.p[i][j][k] * cpml_b_mx_xn.p[i] + cpml_a_mx_xn.p[i]*(Ez.p[iplus][j][k] - Ez.p[i][j][k]);
            }
        }
    }
}

template<class T>
void cpml<T>::updatePsi_hxy_yn(const data3d<T>& Ez) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) //shared(Hz,Ey,Ex,pml,DA,DB,dx,dy)
#endif
    for (unsigned i = 0; i < Psi_hxy_yn.nx; i++) {
        for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
            for (unsigned k = 0; k < Psi_hxy_yn.nz; k++) {
                Psi_hxy_yn.p[i][j][k] = Psi_hxy_yn.p[i][j][k] * cpml_b_my_yn.p[j] + cpml_a_my_yn.p[j]*(Ez.p[i][jplus][k] - Ez.p[i][j][k]);
            }
        }
    }
}

//#include "cpml.cpp"
#endif	/* CPML_H */


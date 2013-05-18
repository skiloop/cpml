/* 
 * File:   cpml.h
 * Author: skiloop
 *
 * Created on 2013年5月14日, 上午9:58
 */

#ifndef CPML_H
#define	CPML_H
#include <math.h>
#include "datastruct.h"

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
template<class type1>
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
     * @param sigmaRatio 
     * @param kappaMax
     * @param alphaMax
     */
    cpml(unsigned short width_xn, unsigned short width_xp, unsigned short width_yn,
            unsigned short width_yp, unsigned short width_zn, unsigned short width_zp,
            unsigned imax, unsigned jmax, unsigned kmax,
            unsigned short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax);

    /**
     * constructor
     * @param cpmlWidth width of cpml layer on all sides
     * @param imax max number of cells in x direction
     * @param jmax max number of cells in y direction
     * @param kmax max number of cells in z direction
     * @param pmlOrder pml order
     * @param sigmaRatio 
     * @param kappaMax
     * @param alphaMax
     */
    cpml(unsigned short cpmlWidth, unsigned imax, unsigned jmax, unsigned kmax,
            unsigned short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax);

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
    unsigned short n_cpml_xp;
    unsigned short n_cpml_xn;
    unsigned short n_cpml_yp;
    unsigned short n_cpml_yn;
    unsigned short n_cpml_zp;
    unsigned short n_cpml_zn;
    // cpml flags
    bool is_cpml_xp;
    bool is_cpml_xn;
    bool is_cpml_yp;
    bool is_cpml_yn;
    bool is_cpml_zp;
    bool is_cpml_zn;
    //============================================
    // cpml arrays
    //============================================
    // xn arrays
    data3d<type1> Psi_eyx_xn;
    data3d<type1> Psi_ezx_xn;
    data3d<type1> Psi_hyx_xn;
    data3d<type1> Psi_hzx_xn;
    // xp arrays
    data3d<type1> Psi_eyx_xp;
    data3d<type1> Psi_ezx_xp;
    data3d<type1> Psi_hyx_xp;
    data3d<type1> Psi_hzx_xp;
    // yn arrays
    data3d<type1> Psi_exy_yn;
    data3d<type1> Psi_ezy_yn;
    data3d<type1> Psi_hxy_yn;
    data3d<type1> Psi_hzy_yn;
    // yp arrays
    data3d<type1> Psi_exy_yp;
    data3d<type1> Psi_ezy_yp;
    data3d<type1> Psi_hxy_yp;
    data3d<type1> Psi_hzy_yp;
    // zn arrays
    data3d<type1> Psi_exz_zn;
    data3d<type1> Psi_eyz_zn;
    data3d<type1> Psi_hxz_zn;
    data3d<type1> Psi_hyz_zn;
    // zp arrays
    data3d<type1> Psi_exz_zp;
    data3d<type1> Psi_eyz_zp;
    data3d<type1> Psi_hxz_zp;
    data3d<type1> Psi_hyz_zp;
    //============================================
    // cpml coefficient arrays
    //============================================
    // xn arrays
    data3d<type1> CPsi_eyx_xn;
    data3d<type1> CPsi_ezx_xn;
    data3d<type1> CPsi_hyx_xn;
    data3d<type1> CPsi_hzx_xn;
    // xp arrays
    data3d<type1> CPsi_eyx_xp;
    data3d<type1> CPsi_ezx_xp;
    data3d<type1> CPsi_hyx_xp;
    data3d<type1> CPsi_hzx_xp;
    // yn arrays
    data3d<type1> CPsi_exy_yn;
    data3d<type1> CPsi_ezy_yn;
    data3d<type1> CPsi_hxy_yn;
    data3d<type1> CPsi_hzy_yn;
    // yp arrays
    data3d<type1> CPsi_exy_yp;
    data3d<type1> CPsi_ezy_yp;
    data3d<type1> CPsi_hxy_yp;
    data3d<type1> CPsi_hzy_yp;
    // zn arrays
    data3d<type1> CPsi_exz_zn;
    data3d<type1> CPsi_eyz_zn;
    data3d<type1> CPsi_hxz_zn;
    data3d<type1> CPsi_hyz_zn;
    // zp arrays
    data3d<type1> CPsi_exz_zp;
    data3d<type1> CPsi_eyz_zp;
    data3d<type1> CPsi_hxz_zp;
    data3d<type1> CPsi_hyz_zp;

    // a and b for cpml to update Psi
    // x direction
    data1d<type1> cpml_a_ex_xn;
    data1d<type1> cpml_b_ex_xn;
    data1d<type1> cpml_a_mx_xn;
    data1d<type1> cpml_b_mx_xn;
    data1d<type1> cpml_a_ex_xp;
    data1d<type1> cpml_b_ex_xp;
    data1d<type1> cpml_a_mx_xp;
    data1d<type1> cpml_b_mx_xp;
    // y direction
    data1d<type1> cpml_a_ey_yn;
    data1d<type1> cpml_b_ey_yn;
    data1d<type1> cpml_a_my_yn;
    data1d<type1> cpml_b_my_yn;
    data1d<type1> cpml_a_ey_yp;
    data1d<type1> cpml_b_ey_yp;
    data1d<type1> cpml_a_my_yp;
    data1d<type1> cpml_b_my_yp;
    // z direction
    data1d<type1> cpml_a_ez_zn;
    data1d<type1> cpml_b_ez_zn;
    data1d<type1> cpml_a_mz_zn;
    data1d<type1> cpml_b_mz_zn;
    data1d<type1> cpml_a_ez_zp;
    data1d<type1> cpml_b_ez_zp;
    data1d<type1> cpml_a_mz_zp;
    data1d<type1> cpml_b_mz_zp;

    //=================================================
    // PUBLIC interface 
    //=================================================
public:
    /**
     * update electric fields in CPML region
     * @param Ex
     * @param Ey
     * @param Ez
     * @param Hx
     * @param Hy
     * @param Hz
     */
    void updateCPML_E_Fields(data3d<type1> &Ex, data3d<type1>& Ey, data3d<type1> &Ez,
            const data3d<type1> & Hx, const data3d<type1> & Hy, const data3d<type1> & Hz);
    /**
     * update magnetic fields in CPML region
     * @param Hx
     * @param Hy
     * @param Hz
     * @param Ex
     * @param Ey
     * @param Ez
     */
    void updateCPML_M_Fields(data3d<type1> &Hx, data3d<type1>& Hy, data3d<type1> &Hz,
            const data3d<type1> & Ex, const data3d<type1> & Ey, const data3d<type1> & Ez);

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
     * @param sigmaRatio
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
    void initCoefficientArrays(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dx, type1 dy, type1 dz,
            data3d<type1>&Ceyhz, data3d<type1>&Cezhy, data3d<type1>&Chyez, data3d<type1>&Chzey,
            data3d<type1>&Cexhz, data3d<type1>&Cezhx, data3d<type1>&Chxez, data3d<type1>&Chzex,
            data3d<type1>&Ceyhx, data3d<type1>&Cexhy, data3d<type1>&Chyex, data3d<type1>&Chxey);
    //=======================================================
    // private functions
    //=======================================================
private:
    /**
     * 
     * @param pmlOrder
     * @param sigmaRatio
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dx
     * @param Ceyhz
     * @param Cezhy
     * @param Chyez
     * @param Chzey
     */
    void initCoefficientArraysXN(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dx, data3d<type1>&Ceyhz, data3d<type1>&Cezhy, data3d<type1>&Chyez, data3d<type1>&Chzey);

    /**
     * 
     * @param pmlOrder
     * @param sigmaRatio
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dx
     * @param Ceyhz
     * @param Cezhy
     * @param Chyez
     * @param Chzey
     */
    void initCoefficientArraysXP(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dx, data3d<type1>&Ceyhz, data3d<type1>&Cezhy, data3d<type1>&Chyez, data3d<type1>&Chzey);

    /**
     * 
     * @param pmlOrder
     * @param sigmaRatio
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dy
     * @param Cexhz
     * @param Cezhx
     * @param Chxez
     * @param Chzex
     */
    void initCoefficientArraysYN(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dy, data3d<type1>&Cexhz, data3d<type1>&Cezhx, data3d<type1>&Chxez, data3d<type1>&Chzex);
    /**
     * 
     * @param pmlOrder
     * @param sigmaRatio
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dy
     * @param Cexhz
     * @param Cezhx
     * @param Chxez
     * @param Chzex
     */
    void initCoefficientArraysYP(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dy, data3d<type1>&Cexhz, data3d<type1>&Cezhx, data3d<type1>&Chxez, data3d<type1>&Chzex);
    /**
     * 
     * @param pmlOrder
     * @param sigmaRatio
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dz
     * @param Ceyhx
     * @param Cexhy
     * @param Chyex
     * @param Chxey
     */
    void initCoefficientArraysZN(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dz, data3d<type1>&Ceyhx, data3d<type1>&Cexhy, data3d<type1>&Chyex, data3d<type1>&Chxey);
    /**
     * 
     * @param pmlOrder
     * @param sigmaRatio
     * @param kappaMax
     * @param alphaMax
     * @param dt
     * @param dz
     * @param Ceyhx
     * @param Cexhy
     * @param Chyex
     * @param Chxey
     */
    void initCoefficientArraysZP(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dz, data3d<type1>&Ceyhx, data3d<type1>&Cexhy, data3d<type1>&Chyex, data3d<type1>&Chxey);

    /**
     * 
     * @param Hx
     * @param Hy
     * @param Hz
     */
    void updatePsiForEFields(const data3d<type1>& Hx, const data3d<type1>& Hy, const data3d<type1>& Hz);
    /**
     * 
     * @param Ex
     * @param Ey
     * @param Ez
     */
    void updatePsiForMFields(const data3d<type1>& Ex, const data3d<type1>& Ey, const data3d<type1>& Ez);

    void updatePsi_eyz_zn(const data3d<type1>& Hx);
    void updatePsi_ezy_yn(const data3d<type1>& Hx);
    void updatePsi_exz_zn(const data3d<type1>& Hy);
    void updatePsi_ezx_xn(const data3d<type1>& Hy);
    void updatePsi_eyx_xn(const data3d<type1>& Hz);
    void updatePsi_exy_yn(const data3d<type1>& Hz);

    void updatePsi_hyz_zn(const data3d<type1>& Ex);
    void updatePsi_hzy_yn(const data3d<type1>& Ex);
    void updatePsi_hxz_zn(const data3d<type1>& Ey);
    void updatePsi_hzx_xn(const data3d<type1>& Ey);
    void updatePsi_hyx_xn(const data3d<type1>& Ez);
    void updatePsi_hxy_yn(const data3d<type1>& Ez);

    void updatePsi_eyz_zp(const data3d<type1>& Hx);
    void updatePsi_ezy_yp(const data3d<type1>& Hx);
    void updatePsi_exz_zp(const data3d<type1>& Hy);
    void updatePsi_ezx_xp(const data3d<type1>& Hy);
    void updatePsi_eyx_xp(const data3d<type1>& Hz);
    void updatePsi_exy_yp(const data3d<type1>& Hz);

    void updatePsi_hyz_zp(const data3d<type1>& Ex);
    void updatePsi_hzy_yp(const data3d<type1>& Ex);
    void updatePsi_hxz_zp(const data3d<type1>& Ey);
    void updatePsi_hzx_xp(const data3d<type1>& Ey);
    void updatePsi_hyx_xp(const data3d<type1>& Ez);
    void updatePsi_hxy_yp(const data3d<type1>& Ez);

    void updateEFieldCPML_x(data3d<type1>&Ey, data3d<type1>&Ez);
    void updateEFieldCPML_y(data3d<type1>&Ex, data3d<type1>&Ez);
    void updateEFieldCPML_z(data3d<type1>&Ex, data3d<type1>&Ey);

    void updateMFieldCPML_x(data3d<type1>&Hy, data3d<type1>&Hz);
    void updateMFieldCPML_y(data3d<type1>&Hx, data3d<type1>&Hz);
    void updateMFieldCPML_z(data3d<type1>&Hx, data3d<type1>&Hy);
};


template<class type1>
cpml<type1>::cpml()
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

template<class type1>
cpml<type1>::cpml(const cpml& orig)
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

template<class type1>
cpml<type1>::~cpml() {
}

template<class type1>
void cpml<type1>::createCPMLArrays(unsigned nx, unsigned ny, unsigned nz) {
    //============================================
    // cpml arrays
    //============================================
    // xn arrays
    unsigned nxp1 = nx + 1;
    unsigned nyp1 = ny + 1;
    unsigned nzp1 = nz + 1;
    if (is_cpml_xn) {
        // x direction
        cpml_a_ex_xn.CreateStruct(n_cpml_xn);
        cpml_b_ex_xn.CreateStruct(n_cpml_xn);
        cpml_a_mx_xn.CreateStruct(n_cpml_xn);
        cpml_b_mx_xn.CreateStruct(n_cpml_xn);
        Psi_eyx_xn.CreateStruct(n_cpml_xn, ny, nzp1);
        Psi_ezx_xn.CreateStruct(n_cpml_xn, nyp1, nz);
        Psi_hyx_xn.CreateStruct(n_cpml_xn, nyp1, nz);
        Psi_hzx_xn.CreateStruct(n_cpml_xn, ny, nzp1);
        CPsi_eyx_xn.CreateStruct(n_cpml_xn, ny, nzp1);
        CPsi_ezx_xn.CreateStruct(n_cpml_xn, nyp1, nz);
        CPsi_hyx_xn.CreateStruct(n_cpml_xn, nyp1, nz);
        CPsi_hzx_xn.CreateStruct(n_cpml_xn, ny, nzp1);
    }
    // xp arrays
    if (is_cpml_xp) {
        cpml_a_ex_xp.CreateStruct(n_cpml_xp);
        cpml_b_ex_xp.CreateStruct(n_cpml_xp);
        cpml_a_mx_xp.CreateStruct(n_cpml_xp);
        cpml_b_mx_xp.CreateStruct(n_cpml_xp);
        Psi_eyx_xp.CreateStruct(n_cpml_xp, ny, nzp1);
        Psi_ezx_xp.CreateStruct(n_cpml_xp, nyp1, nz);
        Psi_hyx_xp.CreateStruct(n_cpml_xp, nyp1, nz);
        Psi_hzx_xp.CreateStruct(n_cpml_xp, ny, nzp1);
        CPsi_eyx_xp.CreateStruct(n_cpml_xp, ny, nzp1);
        CPsi_ezx_xp.CreateStruct(n_cpml_xp, nyp1, nz);
        CPsi_hyx_xp.CreateStruct(n_cpml_xp, nyp1, nz);
        CPsi_hzx_xp.CreateStruct(n_cpml_xp, ny, nzp1);
    }
    // yn arrays
    if (is_cpml_yn) {
        // y direction
        cpml_a_ey_yn.CreateStruct(n_cpml_yn);
        cpml_b_ey_yn.CreateStruct(n_cpml_yn);
        cpml_a_my_yn.CreateStruct(n_cpml_yn);
        cpml_b_my_yn.CreateStruct(n_cpml_yn);
        Psi_exy_yn.CreateStruct(nx, n_cpml_yn, nzp1);
        Psi_ezy_yn.CreateStruct(nxp1, n_cpml_yn, nz);
        Psi_hxy_yn.CreateStruct(nxp1, n_cpml_yn, nz);
        Psi_hzy_yn.CreateStruct(nx, n_cpml_yn, nzp1);
        CPsi_exy_yn.CreateStruct(nx, n_cpml_yn, nzp1);
        CPsi_ezy_yn.CreateStruct(nxp1, n_cpml_yn, nz);
        CPsi_hxy_yn.CreateStruct(nxp1, n_cpml_yn, nz);
        CPsi_hzy_yn.CreateStruct(nx, n_cpml_yn, nzp1);
    }
    // yp arrays
    if (is_cpml_yp) {
        cpml_a_ey_yp.CreateStruct(n_cpml_yp);
        cpml_b_ey_yp.CreateStruct(n_cpml_yp);
        cpml_a_my_yp.CreateStruct(n_cpml_yp);
        cpml_b_my_yp.CreateStruct(n_cpml_yp);
        Psi_exy_yp.CreateStruct(nx, n_cpml_yp, nzp1);
        Psi_ezy_yp.CreateStruct(nxp1, n_cpml_yp, nz);
        Psi_hxy_yp.CreateStruct(nxp1, n_cpml_yp, nz);
        Psi_hzy_yp.CreateStruct(nx, n_cpml_yp, nzp1);
        CPsi_exy_yp.CreateStruct(nx, n_cpml_yp, nzp1);
        CPsi_ezy_yp.CreateStruct(nxp1, n_cpml_yp, nz);
        CPsi_hxy_yp.CreateStruct(nxp1, n_cpml_yp, nz);
        CPsi_hzy_yp.CreateStruct(nx, n_cpml_yp, nzp1);
    }
    // zn arrays
    if (is_cpml_zn) {
        // z direction
        cpml_a_ez_zn.CreateStruct(n_cpml_zn);
        cpml_b_ez_zn.CreateStruct(n_cpml_zn);
        cpml_a_mz_zn.CreateStruct(n_cpml_zn);
        cpml_b_mz_zn.CreateStruct(n_cpml_zn);
        Psi_exz_zn.CreateStruct(nx, nyp1, n_cpml_zn);
        Psi_eyz_zn.CreateStruct(nxp1, ny, n_cpml_zn);
        Psi_hxz_zn.CreateStruct(nxp1, ny, n_cpml_zn);
        Psi_hyz_zn.CreateStruct(nx, nyp1, n_cpml_zn);
        CPsi_exz_zn.CreateStruct(nx, nyp1, n_cpml_zn);
        CPsi_eyz_zn.CreateStruct(nxp1, ny, n_cpml_zn);
        CPsi_hxz_zn.CreateStruct(nxp1, ny, n_cpml_zn);
        CPsi_hyz_zn.CreateStruct(nx, nyp1, n_cpml_zn);
    }
    // zp arrays
    if (is_cpml_zp) {
        cpml_a_ez_zp.CreateStruct(n_cpml_zp);
        cpml_b_ez_zp.CreateStruct(n_cpml_zp);
        cpml_a_mz_zp.CreateStruct(n_cpml_zp);
        cpml_b_mz_zp.CreateStruct(n_cpml_zp);
        Psi_exz_zp.CreateStruct(nx, nyp1, n_cpml_zp);
        Psi_eyz_zp.CreateStruct(nxp1, ny, n_cpml_zp);
        Psi_hxz_zp.CreateStruct(nxp1, ny, n_cpml_zp);
        Psi_hyz_zp.CreateStruct(nx, nyp1, n_cpml_zp);
        CPsi_exz_zp.CreateStruct(nx, nyp1, n_cpml_zp);
        CPsi_eyz_zp.CreateStruct(nxp1, ny, n_cpml_zp);
        CPsi_hxz_zp.CreateStruct(nxp1, ny, n_cpml_zp);
        CPsi_hyz_zp.CreateStruct(nx, nyp1, n_cpml_zp);
    }
}

template<class type1>
void cpml<type1>::setCPMLRegion(short pmlWidth) {
    setCPMLRegion(pmlWidth, pmlWidth, pmlWidth, pmlWidth, pmlWidth, pmlWidth);
}

template<class type1>
void cpml<type1>::setCPMLRegion(short width_xn, short width_xp, short width_yn, short width_yp, short width_zn, short width_zp) {
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

template<class type1>
void cpml<type1>::initCoefficientArrays(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dx, type1 dy, type1 dz,
        data3d<type1>&Ceyhz, data3d<type1>&Cezhy, data3d<type1>&Chyez, data3d<type1>&Chzey,
        data3d<type1>&Cexhz, data3d<type1>&Cezhx, data3d<type1>&Chxez, data3d<type1>&Chzex,
        data3d<type1>&Ceyhx, data3d<type1>&Cexhy, data3d<type1>&Chyex, data3d<type1>&Chxey) {
    initCoefficientArraysXN(pmlOrder, sigmaRatio, kappaMax, alphaMax, dt, dx, Ceyhz, Cezhy, Chyez, Chzey);
    initCoefficientArraysXP(pmlOrder, sigmaRatio, kappaMax, alphaMax, dt, dx, Ceyhz, Cezhy, Chyez, Chzey);
    initCoefficientArraysYN(pmlOrder, sigmaRatio, kappaMax, alphaMax, dt, dy, Cexhz, Cezhx, Chxez, Chzex);
    initCoefficientArraysYP(pmlOrder, sigmaRatio, kappaMax, alphaMax, dt, dy, Cexhz, Cezhx, Chxez, Chzex);
    initCoefficientArraysZN(pmlOrder, sigmaRatio, kappaMax, alphaMax, dt, dz, Ceyhx, Cexhy, Chyex, Chxey);
    initCoefficientArraysZP(pmlOrder, sigmaRatio, kappaMax, alphaMax, dt, dz, Ceyhx, Cexhy, Chyex, Chxey);
}

template<class type1>
void cpml<type1>::initCoefficientArraysXN(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dx,
        data3d<type1>&Ceyhz, data3d<type1>&Cezhy, data3d<type1>&Chyez, data3d<type1>&Chzey) {
    if (is_cpml_xn) {
        type1 sigmaMax = sigmaRatio * (pmlOrder + 1) / (150 * M_PI * dx);

        for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
            type1 rho_e = (n_cpml_xn - i - 0.75) / n_cpml_xn;
            type1 rho_m = (n_cpml_xn - i - 0.25) / n_cpml_xn;
            type1 rho_e_pmlOrder = pow(rho_e, pmlOrder);
            type1 rho_m_pmlOrder = pow(rho_m, pmlOrder);
            type1 sigma_pex = sigmaMax*rho_e_pmlOrder;
            type1 sigma_pmx = sigmaMax * rho_m_pmlOrder*Mu0DivEps0;
            type1 kappa_ex = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            type1 kappa_mx = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            type1 alpha_ex = alphaMax*rho_e_pmlOrder;
            type1 alpha_mx = alphaMax*rho_m_pmlOrder;
            cpml_b_ex_xn.p[i] = exp((-dt / eps_0) * sigma_pex / kappa_ex + alpha_ex);
            cpml_b_mx_xn.p[i] = exp((-dt / mu_0) * sigma_pmx / kappa_mx + alpha_mx);
            cpml_a_ex_xn.p[i] = 1 / dx * (cpml_b_ex_xn.p[i] - 1.0) * sigma_pex / (kappa_ex * (sigma_pex + kappa_ex * alpha_ex));
            cpml_a_mx_xn.p[i] = 1 / dx * (cpml_b_mx_xn.p[i] - 1.0) * sigma_pmx / (kappa_mx * (sigma_pmx + kappa_mx * alpha_mx));
            for (unsigned j = 0; j < Psi_eyx_xn.ny; j++) {
                for (unsigned k = 0; k < Psi_eyx_xn.nz; k++) {
                    CPsi_eyx_xn.p[i][j][k] = Ceyhz.p[iplus][j][k] * dx;
                    CPsi_ezx_xn.p[i][k][j] = Cezhy.p[iplus][k][j] * dx;
                    Ceyhz.p[iplus][j][k] = Ceyhz.p[iplus][j][k] / kappa_ex;
                    Cezhy.p[iplus][k][j] = Cezhy.p[iplus][k][j] / kappa_ex;
                }
            }
            for (unsigned j = 0; j < Psi_hyx_xn.ny; j++) {
                for (unsigned k = 0; k < Psi_hyx_xn.nz; k++) {
                    CPsi_hyx_xn.p[i][j][k] = Chyez.p[i][j][k] * dx;
                    CPsi_hzx_xn.p[i][k][j] = Chzey.p[i][k][j] * dx;
                    Chyez.p[i][j][k] = Chyez.p[i][j][k] / kappa_mx;
                    Chzey.p[i][k][j] = Chzey.p[i][k][j] / kappa_mx;
                }
            }
        }
    }
}

template<class type1>
void cpml<type1>::initCoefficientArraysXP(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dx,
        data3d<type1>&Ceyhz, data3d<type1>&Cezhy, data3d<type1>&Chyez, data3d<type1>&Chzey) {
    if (is_cpml_xp) {
        type1 sigmaMax = sigmaRatio * (pmlOrder + 1) / (150 * M_PI * dx);
        unsigned iex = Ceyhz.nx - n_cpml_xp - 1;
        unsigned ihx = Chyez.nx - n_cpml_xp;
        for (unsigned i = 0; i < n_cpml_xp; i++) {
            type1 rho_e = (i + 0.25) / n_cpml_xp;
            type1 rho_m = (i + 0.75) / n_cpml_xp;
            type1 rho_e_pmlOrder = pow(rho_e, pmlOrder);
            type1 rho_m_pmlOrder = pow(rho_m, pmlOrder);
            type1 sigma_pex = sigmaMax*rho_e_pmlOrder;
            type1 sigma_pmx = sigmaMax * rho_m_pmlOrder*Mu0DivEps0;
            type1 kappa_ex = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            type1 kappa_mx = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            type1 alpha_ex = alphaMax*rho_e_pmlOrder;
            type1 alpha_mx = alphaMax*rho_m_pmlOrder;
            cpml_b_ex_xp.p[i] = exp((-dt / eps_0) * sigma_pex / kappa_ex + alpha_ex);
            cpml_b_mx_xp.p[i] = exp((-dt / mu_0) * sigma_pmx / kappa_mx + alpha_mx);
            cpml_a_ex_xp.p[i] = 1 / dx * (cpml_b_ex_xp.p[i] - 1.0) * sigma_pex / (kappa_ex * (sigma_pex + kappa_ex * alpha_ex));
            cpml_a_mx_xp.p[i] = 1 / dx * (cpml_b_mx_xp.p[i] - 1.0) * sigma_pmx / (kappa_mx * (sigma_pmx + kappa_mx * alpha_mx));

            for (unsigned j = 0; j < Psi_eyx_xp.ny; j++) {
                for (unsigned k = 0; k < Psi_eyx_xp.nz; k++) {
                    CPsi_eyx_xp.p[i][j][k] = Ceyhz.p[iex][j][k] * dx;
                    CPsi_ezx_xp.p[i][k][j] = Cezhy.p[iex][k][j] * dx;
                    Ceyhz.p[iex][j][k] = Ceyhz.p[iex][j][k] / kappa_ex;
                    Cezhy.p[iex][k][j] = Cezhy.p[iex][k][j] / kappa_ex;
                }
            }
            for (unsigned j = 0; j < Psi_hyx_xp.ny; j++) {
                for (unsigned k = 0; k < Psi_hyx_xp.nz; k++) {
                    CPsi_hyx_xp.p[i][j][k] = Chyez.p[ihx][j][k] * dx;
                    CPsi_hzx_xp.p[i][k][j] = Chzey.p[ihx][k][j] * dx;
                    Chyez.p[ihx][j][k] = Chyez.p[ihx][j][k] / kappa_mx;
                    Chzey.p[ihx][k][j] = Chzey.p[ihx][k][j] / kappa_mx;
                }
            }
            ihx++;
            iex++;
        }
    }
}

template<class type1>
void cpml<type1>::initCoefficientArraysYN(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dy,
        data3d<type1>&Cexhz, data3d<type1>&Cezhx, data3d<type1>&Chxez, data3d<type1>&Chzex) {
    if (is_cpml_yn) {
        type1 sigmaMax = sigmaRatio * (pmlOrder + 1) / (150 * M_PI * dy);
        for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
            type1 rho_e = (n_cpml_yn - j - 0.75) / n_cpml_yn;
            type1 rho_m = (n_cpml_yn - j - 0.25) / n_cpml_yn;
            type1 rho_e_pmlOrder = pow(rho_e, pmlOrder);
            type1 rho_m_pmlOrder = pow(rho_m, pmlOrder);
            type1 sigma_pey = sigmaMax*rho_e_pmlOrder;
            type1 sigma_pmy = sigmaMax * rho_m_pmlOrder*Mu0DivEps0;
            type1 kappa_ey = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            type1 kappa_my = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            type1 alpha_ey = alphaMax*rho_e_pmlOrder;
            type1 alpha_my = alphaMax*rho_m_pmlOrder;
            cpml_b_ey_yn.p[j] = exp((-dt / eps_0) * sigma_pey / kappa_ey + alpha_ey);
            cpml_b_my_yn.p[j] = exp((-dt / mu_0) * sigma_pmy / kappa_my + alpha_my);
            cpml_a_ey_yn.p[j] = 1 / dy * (cpml_b_ey_yn.p[j] - 1.0) * sigma_pey / (kappa_ey * (sigma_pey + kappa_ey * alpha_ey));
            cpml_a_my_yn.p[j] = 1 / dy * (cpml_b_my_yn.p[j] - 1.0) * sigma_pmy / (kappa_my * (sigma_pmy + kappa_my * alpha_my));

            for (unsigned i = 0; i < Psi_exy_yn.nx; i++) {
                for (unsigned k = 0; k < Psi_exy_yn.nz; k++) {
                    CPsi_exy_yn.p[i][j][k] = Cexhz.p[i][jplus][k] * dy;
                    CPsi_ezy_yn.p[j][i][k] = Cezhx.p[k][jplus][i] * dy;
                    Cexhz.p[i][jplus][k] = Cexhz.p[i][jplus][k] / kappa_ey;
                    Cezhx.p[k][jplus][i] = Cezhx.p[k][jplus][i] / kappa_ey;
                }
            }
            for (unsigned i = 0; i < Psi_hxy_yn.nx; i++) {
                for (unsigned k = 0; k < Psi_hxy_yn.nz; k++) {
                    CPsi_hxy_yn.p[i][j][k] = Chxez.p[i][j][k] * dy;
                    CPsi_hzy_yn.p[k][j][i] = Chzex.p[k][j][i] * dy;
                    Chxez.p[i][j][k] = Chxez.p[i][j][k] / kappa_my;
                    Chzex.p[k][j][i] = Chzex.p[k][j][i] / kappa_my;
                }
            }
        }
    }
}

template<class type1>
void cpml<type1>::initCoefficientArraysYP(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dy,
        data3d<type1>&Cexhz, data3d<type1>&Cezhx, data3d<type1>&Chxez, data3d<type1>&Chzex) {
    if (is_cpml_yp) {
        type1 sigmaMax = sigmaRatio * (pmlOrder + 1) / (150 * M_PI * dy);
        unsigned iex = Psi_exy_yp.ny - n_cpml_yp - 1;
        unsigned ihx = Psi_hxy_yp.ny - n_cpml_yp - 1;
        for (unsigned j = 0; j < n_cpml_yp; j++) {
            type1 rho_e = (j + 0.25) / n_cpml_yp;
            type1 rho_m = (j + 0.75) / n_cpml_yp;
            type1 rho_e_pmlOrder = pow(rho_e, pmlOrder);
            type1 rho_m_pmlOrder = pow(rho_m, pmlOrder);
            type1 sigma_pey = sigmaMax*rho_e_pmlOrder;
            type1 sigma_pmy = sigmaMax * rho_m_pmlOrder*Mu0DivEps0;
            type1 kappa_ey = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            type1 kappa_my = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            type1 alpha_ey = alphaMax*rho_e_pmlOrder;
            type1 alpha_my = alphaMax*rho_m_pmlOrder;
            cpml_b_ey_yp.p[j] = exp((-dt / eps_0) * sigma_pey / kappa_ey + alpha_ey);
            cpml_b_my_yp.p[j] = exp((-dt / mu_0) * sigma_pmy / kappa_my + alpha_my);
            cpml_a_ey_yp.p[j] = 1 / dy * (cpml_b_ey_yp.p[j] - 1.0) * sigma_pey / (kappa_ey * (sigma_pey + kappa_ey * alpha_ey));
            cpml_a_my_yp.p[j] = 1 / dy * (cpml_b_my_yp.p[j] - 1.0) * sigma_pmy / (kappa_my * (sigma_pmy + kappa_my * alpha_my));

            for (unsigned i = 0; i < Psi_exy_yp.nx; i++) {
                for (unsigned k = 0; k < Psi_exy_yp.nz; k++) {
                    // TODO fix region size
                    CPsi_exy_yp.p[i][j][k] = Cexhz.p[i][iex][k] * dy;
                    CPsi_ezy_yp.p[k][j][i] = Cezhx.p[k][iex][i] * dy;
                    Cexhz.p[i][iex][k] = Cexhz.p[i][iex][k] / kappa_ey;
                    Cezhx.p[k][iex][i] = Cezhx.p[k][iex][i] / kappa_ey;
                }
            }
            for (unsigned i = 0; i < Psi_hxy_yp.nx; i++) {
                for (unsigned k = 0; k < Psi_hxy_yp.nz; k++) {
                    CPsi_hxy_yp.p[i][j][k] = Chxez.p[i][ihx][k] * dy;
                    CPsi_hzy_yp.p[k][j][i] = Chzex.p[k][ihx][i] * dy;
                    Chxez.p[i][ihx][k] = Chxez.p[i][ihx][k] / kappa_my;
                    Chzex.p[k][ihx][i] = Chzex.p[k][ihx][i] / kappa_my;
                }
            }
            ihx++;
            iex++;
        }
    }
}

template<class type1>
void cpml<type1>::initCoefficientArraysZN(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dz,
        data3d<type1>&Ceyhx, data3d<type1>&Cexhy, data3d<type1>&Chyex, data3d<type1>&Chxey) {
    if (is_cpml_zn) {
        type1 sigmaMax = sigmaRatio * (pmlOrder + 1) / (150 * M_PI * dz);
        for (unsigned k = 0, iplus = 1; k < n_cpml_zn; k++, iplus++) {
            type1 rho_e = (n_cpml_zn - k - 0.75) / n_cpml_zn;
            type1 rho_m = (n_cpml_zn - k - 0.25) / n_cpml_zn;
            type1 rho_e_pmlOrder = pow(rho_e, pmlOrder);
            type1 rho_m_pmlOrder = pow(rho_m, pmlOrder);
            type1 sigma_pez = sigmaMax*rho_e_pmlOrder;
            type1 sigma_pmz = sigmaMax * rho_m_pmlOrder*Mu0DivEps0;
            type1 kappa_ez = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            type1 kappa_mz = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            type1 alpha_ez = alphaMax*rho_e_pmlOrder;
            type1 alpha_mz = alphaMax*rho_m_pmlOrder;
            cpml_b_ez_zn.p[k] = exp((-dt / eps_0) * sigma_pez / kappa_ez + alpha_ez);
            cpml_b_mz_zn.p[k] = exp((-dt / mu_0) * sigma_pmz / kappa_mz + alpha_mz);
            cpml_a_ez_zn.p[k] = 1 / dz * (cpml_b_ez_zn.p[k] - 1.0) * sigma_pez / (kappa_ez * (sigma_pez + kappa_ez * alpha_ez));
            cpml_a_mz_zn.p[k] = 1 / dz * (cpml_b_mz_zn.p[k] - 1.0) * sigma_pmz / (kappa_mz * (sigma_pmz + kappa_mz * alpha_mz));

            for (unsigned j = 0; j < Psi_exz_zn.ny; j++) {
                for (unsigned i = 0; i < Psi_exz_zn.nx; i++) {
                    // TODO fix region size
                    CPsi_exz_zn.p[i][j][k] = Cexhy.p[i][j][iplus] * dz;
                    CPsi_eyz_zn.p[j][i][k] = Ceyhx.p[j][i][iplus] * dz;
                    Cexhy.p[i][j][iplus] = Cexhy.p[i][j][iplus] / kappa_ez;
                    Ceyhx.p[j][i][iplus] = Ceyhx.p[j][i][iplus] / kappa_ez;
                }
            }
            for (unsigned j = 0; j < Psi_hxz_zn.ny; j++) {
                for (unsigned i = 0; i < Psi_hxz_zn.nx; i++) {
                    CPsi_hxz_zn.p[i][j][k] = Chxey.p[i][j][k] * dz;
                    CPsi_hyz_zn.p[j][i][k] = Chyex.p[j][i][k] * dz;
                    Chxey.p[i][j][k] = Chxey.p[i][j][k] / kappa_mz;
                    Chyex.p[j][i][k] = Chyex.p[j][i][k] / kappa_mz;
                }
            }
        }
    }
}

template<class type1>
void cpml<type1>::initCoefficientArraysZP(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dz,
        data3d<type1>&Ceyhx, data3d<type1>&Cexhy, data3d<type1>&Chyex, data3d<type1>&Chxey) {
    if (is_cpml_zp) {
        type1 sigmaMax = sigmaRatio * (pmlOrder + 1) / (150 * M_PI * dz);
        unsigned iez = Psi_eyz_zp.nz - n_cpml_zp - 1;
        unsigned ihz = Psi_hyz_zp.nz - n_cpml_zp;
        for (unsigned k = 0; k < n_cpml_zp; k++) {
            type1 rho_e = (k + 0.25) / n_cpml_zp;
            type1 rho_m = (k + 0.75) / n_cpml_zp;
            type1 rho_e_pmlOrder = pow(rho_e, pmlOrder);
            type1 rho_m_pmlOrder = pow(rho_m, pmlOrder);
            type1 sigma_pez = sigmaMax*rho_e_pmlOrder;
            type1 sigma_pmz = sigmaMax * rho_m_pmlOrder*Mu0DivEps0;
            type1 kappa_ez = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            type1 kappa_mz = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            type1 alpha_ez = alphaMax*rho_e_pmlOrder;
            type1 alpha_mz = alphaMax*rho_m_pmlOrder;
            cpml_b_ez_zp.p[k] = exp((-dt / eps_0) * sigma_pez / kappa_ez + alpha_ez);
            cpml_b_mz_zp.p[k] = exp((-dt / mu_0) * sigma_pmz / kappa_mz + alpha_mz);
            cpml_a_ez_zp.p[k] = 1 / dz * (cpml_b_ez_zp.p[k] - 1.0) * sigma_pez / (kappa_ez * (sigma_pez + kappa_ez * alpha_ez));
            cpml_a_mz_zp.p[k] = 1 / dz * (cpml_b_mz_zp.p[k] - 1.0) * sigma_pmz / (kappa_mz * (sigma_pmz + kappa_mz * alpha_mz));

            for (unsigned j = 0; j < Psi_eyz_zp.ny; j++) {
                for (unsigned i = 0; i < Psi_eyz_zp.nx; i++) {
                    // TODO fix region size
                    CPsi_eyz_zp.p[i][j][k] = Ceyhx.p[i][j][iez] * dz;
                    CPsi_exz_zp.p[j][i][k] = Cexhy.p[j][i][iez] * dz;
                    Ceyhx.p[i][j][iez] = Ceyhx.p[i][j][iez] / kappa_ez;
                    Cexhy.p[j][i][iez] = Cexhy.p[j][i][iez] / kappa_ez;
                }
            }
            for (unsigned j = 0; j < Psi_hyz_zp.ny; j++) {
                for (unsigned i = 0; i < Psi_hyz_zp.nx; i++) {
                    CPsi_hyz_zp.p[i][j][k] = Chyex.p[i][j][ihz] * dz;
                    CPsi_hxz_zp.p[j][i][k] = Chxey.p[j][i][ihz] * dz;
                    Chyex.p[i][j][ihz] = Chyex.p[i][j][ihz] / kappa_mz;
                    Chxey.p[j][i][ihz] = Chxey.p[j][i][ihz] / kappa_mz;
                }
            }
            ihz++;
            iez++;
        }
    }
}

template<class type1>
void cpml<type1>::updateCPML_E_Fields(data3d<type1>& Ex, data3d<type1>& Ey, data3d<type1>& Ez,
        const data3d<type1>& Hx, const data3d<type1>& Hy, const data3d<type1>& Hz) {
    updatePsiForEFields(Hx, Hy, Hz);
    updateEFieldCPML_x(Ey, Ez);
    updateEFieldCPML_y(Ex, Ez);
    updateEFieldCPML_z(Ex, Ey);
}

template<class type1>
void cpml<type1>::updateCPML_M_Fields(data3d<type1>& Hx, data3d<type1>& Hy, data3d<type1>& Hz,
        const data3d<type1>& Ex, const data3d<type1>& Ey, const data3d<type1>& Ez) {
    updatePsiForMFields(Ex, Ey, Ez);
    updateMFieldCPML_x(Hy, Hz);
    updateMFieldCPML_y(Hx, Hz);
    updateMFieldCPML_z(Hx, Hy);
}

template<class type1>
void cpml<type1>::updatePsiForEFields(const data3d<type1>& Hx, const data3d<type1>& Hy, const data3d<type1>& Hz) {
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

template<class type1>
void cpml<type1>::updatePsiForMFields(const data3d<type1>& Ex, const data3d<type1>& Ey, const data3d<type1>& Ez) {
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

template<class type1>
void cpml<type1>::updateEFieldCPML_x(data3d<type1>& Ey, data3d<type1>& Ez) {
    // x negetive region
    if (is_cpml_xn) {
        for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
            for (unsigned j = 0; j < Ey.ny; j++) {
                for (unsigned k = 0; k < Ey.nz; k++) {
                    Ey.p[iplus][j][k] += CPsi_hyx_xn.p[i][j][k] * Psi_hyx_xn.p[i][j][k];
                }
            }
            for (unsigned j = 0; j < Ez.ny; j++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ez.p[iplus][j][k] += CPsi_hzx_xn.p[i][j][k] * Psi_hzx_xn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_xp) {
        unsigned ihz = Ez.nx - n_cpml_xp - 1;
        unsigned ihy = Ey.nx - n_cpml_xp - 1;
        for (unsigned i = 0; i < n_cpml_xp; i++) {
            for (unsigned j = 0; j < Ey.ny; j++) {
                for (unsigned k = 0; k < Ey.nz; k++) {
                    Ey.p[ihy][j][k] += CPsi_hyx_xp.p[i][j][k] * Psi_hyx_xp.p[i][j][k];
                }
            }
            for (unsigned j = 0; j < Ez.ny; j++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ez.p[ihz][j][k] += CPsi_hzx_xp.p[i][j][k] * Psi_hzx_xp.p[i][j][k];
                }
            }
            ihy++;
            ihz++;
        }
    }
}

template<class type1>
void cpml<type1>::updateEFieldCPML_y(data3d<type1>& Ex, data3d<type1>& Ez) {
    // x negetive region
    if (is_cpml_yn) {
        for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
            for (unsigned i = 0; i < Ex.nx; i++) {
                for (unsigned k = 0; k < Ex.nz; k++) {
                    Ex.p[i][jplus][k] += CPsi_hxy_yn.p[i][j][k] * Psi_hxy_yn.p[i][j][k];
                }
            }
            for (unsigned i = 0; i < Ez.nx; i++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ez.p[i][jplus][k] += CPsi_hzy_yn.p[i][j][k] * Psi_hzy_yn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_yp) {
        unsigned jhz = Ez.ny - n_cpml_yp - 1;
        unsigned jhx = Ex.ny - n_cpml_yp - 1;
        for (unsigned j = 0; j < n_cpml_yp; j++) {
            for (unsigned i = 0; i < Ex.nx; i++) {
                for (unsigned k = 0; k < Ex.nz; k++) {
                    Ex.p[i][jhx][k] += CPsi_hxy_yp.p[i][j][k] * Psi_hxy_yp.p[i][j][k];
                }
            }
            for (unsigned i = 0; i < Ez.nx; i++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ez.p[i][jhz][k] += CPsi_hzy_yp.p[i][j][k] * Psi_hzy_yp.p[i][j][k];
                }
            }
            jhz++;
            jhx++;
        }
    }
}

template<class type1>
void cpml<type1>::updateEFieldCPML_z(data3d<type1>& Ex, data3d<type1>& Ey) {
    // x negetive region
    if (is_cpml_zn) {
        for (unsigned k = 0, kplus = 1; k < n_cpml_zn; k++, kplus++) {
            for (unsigned j = 0; j < Ey.ny; j++) {
                for (unsigned i = 0; i < Ey.nx; i++) {
                    Ey.p[i][j][kplus] += CPsi_hyz_zn.p[i][j][k] * Psi_hyz_zn.p[i][j][k];
                }
            }
            for (unsigned j = 0; j < Ex.ny; j++) {
                for (unsigned i = 0; i < Ex.nx; i++) {
                    Ex.p[i][j][kplus] += CPsi_hxz_zn.p[i][j][k] * Psi_hxz_zn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_zp) {
        unsigned khy = Ey.nz - n_cpml_zp - 1;
        unsigned khx = Ex.nz - n_cpml_zp - 1;
        for (unsigned k = 0; k < n_cpml_zn; k++) {
            for (unsigned j = 0; j < Ey.ny; j++) {
                for (unsigned i = 0; i < Ey.nx; i++) {
                    Ey.p[i][j][khy] += CPsi_hyz_zn.p[i][j][k] * Psi_hyz_zn.p[i][j][k];
                }
            }
            for (unsigned j = 0; j < Ex.ny; j++) {
                for (unsigned i = 0; i < Ex.nx; i++) {
                    Ex.p[i][j][khx] += CPsi_hxz_zn.p[i][j][k] * Psi_hxz_zn.p[i][j][k];
                }
            }
            khx++;
            khy++;
        }
    }
}

template<class type1>
void cpml<type1>::updateMFieldCPML_x(data3d<type1>& Hy, data3d<type1>& Hz) {
    // x negetive region
    if (is_cpml_xn) {
        for (unsigned i = 0; i < n_cpml_xn; i++) {
            for (unsigned j = 0; j < Hy.ny; j++) {
                for (unsigned k = 0; k < Hy.nz; k++) {
                    Hy.p[i][j][k] += CPsi_hyx_xn.p[i][j][k] * Psi_hyx_xn.p[i][j][k];
                }
            }
            for (unsigned j = 0; j < Hz.ny; j++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hz.p[i][j][k] += CPsi_hzx_xn.p[i][j][k] * Psi_hzx_xn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_xp) {
        unsigned ihz = Hz.nx - n_cpml_xp;
        unsigned ihy = Hy.nx - n_cpml_xp;
        for (unsigned i = 0; i < n_cpml_xp; i++) {
            for (unsigned j = 0; j < Hy.ny; j++) {
                for (unsigned k = 0; k < Hy.nz; k++) {
                    Hy.p[ihy][j][k] += CPsi_hyx_xp.p[i][j][k] * Psi_hyx_xp.p[i][j][k];
                }
            }
            for (unsigned j = 0; j < Hz.ny; j++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hz.p[ihz][j][k] += CPsi_hzx_xp.p[i][j][k] * Psi_hzx_xp.p[i][j][k];
                }
            }
            ihy++;
            ihz++;
        }
    }
}

template<class type1>
void cpml<type1>::updateMFieldCPML_y(data3d<type1>& Hx, data3d<type1>& Hz) {
    // x negetive region
    if (is_cpml_yn) {
        for (unsigned j = 0; j < n_cpml_yn; j++) {
            for (unsigned i = 0; i < Hx.nx; i++) {
                for (unsigned k = 0; k < Hx.nz; k++) {
                    Hx.p[i][j][k] += CPsi_hxy_yn.p[i][j][k] * Psi_hxy_yn.p[i][j][k];
                }
            }
            for (unsigned i = 0; i < Hz.nx; i++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hz.p[i][j][k] += CPsi_hzy_yn.p[i][j][k] * Psi_hzy_yn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_yp) {
        unsigned ihx = Hx.ny - n_cpml_yp;
        unsigned ihz = Hz.ny - n_cpml_yp;
        for (unsigned j = 0; j < n_cpml_yp; j++) {
            for (unsigned i = 0; i < Hx.nx; i++) {
                for (unsigned k = 0; k < Hx.nz; k++) {
                    Hx.p[i][ihx][k] += CPsi_hxy_yp.p[i][j][k] * Psi_hxy_yp.p[i][j][k];
                }
            }
            for (unsigned i = 0; i < Hz.nx; i++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hz.p[i][ihz][k] += CPsi_hzy_yp.p[i][j][k] * Psi_hzy_yp.p[i][j][k];
                }
            }
            ihx++;
            ihz++;
        }
    }
}

template<class type1>
void cpml<type1>::updateMFieldCPML_z(data3d<type1>& Hx, data3d<type1>& Hy) {
    //TODO add update Hz in pml region
    // x negetive region
    if (is_cpml_zn) {
        for (unsigned k = 0; k < n_cpml_zn; k++) {
            for (unsigned j = 0; j < Hy.ny; j++) {
                for (unsigned i = 0; i < Hy.nx; i++) {
                    Hy.p[i][j][k] += CPsi_hyz_zn.p[i][j][k] * Psi_hyz_zn.p[i][j][k];
                }
            }
            for (unsigned j = 0; j < Hx.ny; j++) {
                for (unsigned i = 0; i < Hx.nx; i++) {
                    Hx.p[i][j][k] += CPsi_hxz_zn.p[i][j][k] * Psi_hxz_zn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_zp) {
        unsigned khx = Hx.nz - n_cpml_zp;
        unsigned khy = Hy.nz - n_cpml_zp;
        for (unsigned k = 0; k < n_cpml_zp; k++) {
            for (unsigned i = 0; i < Hy.nx; i++) {
                for (unsigned j = 0; j < Hy.ny; j++) {
                    Hy.p[i][j][khy] += CPsi_hyz_zn.p[i][j][k] * Psi_hyz_zn.p[i][j][k];
                }
            }
            for (unsigned i = 0; i < Hx.nx; i++) {
                for (unsigned j = 0; j < Hx.ny; j++) {
                    Hx.p[i][j][khx] += CPsi_hxz_zn.p[i][j][k] * Psi_hxz_zn.p[i][j][k];
                }
            }
            khx++;
            khy++;
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_eyz_zp(const data3d<type1>& Hx) {
    for (unsigned k = 0, ikz = Hx.nz - n_cpml_zp; k < n_cpml_zp; k++, ikz++) {
        for (unsigned i = 0; i < Psi_eyz_zp.nx; i++) {
            for (unsigned j = 0; j < Psi_eyz_zp.ny; j++) {
                Psi_eyz_zp.p[i][j][k] = Psi_eyz_zp.p[i][j][k] * cpml_b_ez_zp.p[k] + cpml_a_ez_zp.p[k]*(Hx.p[i][j][ikz] - Hx.p[i][j][ikz - 1]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_ezy_yp(const data3d<type1>& Hx) {
    for (unsigned j = 0, jhy = Hx.ny - n_cpml_yp; j < n_cpml_yp; j++, jhy++) {
        for (unsigned i = 0; i < Psi_ezy_yp.nx; i++) {
            for (unsigned k = 0; k < Psi_ezy_yp.nz; k++) {
                Psi_ezy_yp.p[i][j][k] = Psi_ezy_yp.p[i][j][k] * cpml_b_ey_yp.p[j] + cpml_a_ey_yp.p[j]*(Hx.p[i][jhy][k] - Hx.p[i][jhy - 1][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_exz_zp(const data3d<type1>& Hy) {
    for (unsigned k = 0, iez = Hy.nz - n_cpml_zp; k < n_cpml_zp; k++, iez++) {
        for (unsigned i = 0; i < Psi_exz_zp.nx; i++) {
            for (unsigned j = 0; j < Psi_exz_zp.ny; j++) {
                Psi_exz_zp.p[i][j][k] = cpml_b_ez_zp.p[k] * Psi_exz_zp.p[i][j][k] + cpml_a_ez_zp.p[k]*(Hy.p[i][j][iez ] - Hy.p[i][j][iez - 1]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_ezx_xp(const data3d<type1>& Hy) {
    for (unsigned i = 0, iex = Hy.nx - n_cpml_xp; i < n_cpml_xp; i++, iex++) {
        for (unsigned j = 0; j < Psi_ezx_xp.ny; j++) {
            for (unsigned k = 0; k < Psi_ezx_xp.nz; k++) {
                Psi_ezx_xp.p[i][j][k] = cpml_b_ex_xp.p[i] * Psi_ezx_xp.p[i][j][k] + cpml_a_ex_xp.p[i]*(Hy.p[iex][j][k] - Hy.p[iex - 1][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_eyx_xp(const data3d<type1>& Hz) {
    for (unsigned i = 0, ihz = Hz.nx - n_cpml_xp; i < n_cpml_xp; i++, ihz++) {
        for (unsigned j = 0; j < Psi_eyx_xp.ny; j++) {
            for (unsigned k = 0; k < Psi_eyx_xp.nz; k++) {
                Psi_eyx_xp.p[i][j][k] = cpml_b_ex_xp.p[i] * Psi_eyx_xp.p[i][j][k] + cpml_a_ex_xp.p[i]*(Hz.p[ihz][j][k] - Hz.p[ihz - 1][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_exy_yp(const data3d<type1>& Hz) {
    for (unsigned j = 0, imy = Hz.ny - n_cpml_yp; j < n_cpml_yp; j++, imy++) {
        for (unsigned i = 0; i < Psi_exy_yp.nx; i++) {
            for (unsigned k = 0; k < Psi_exy_yp.nz; k++) {
                Psi_exy_yp.p[i][j][k] = Psi_exy_yp.p[i][j][k] * cpml_b_ey_yp.p[j] + cpml_a_ey_yp.p[j]*(Hz.p[i][imy][k] - Hz.p[i][imy - 1][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_hyz_zp(const data3d<type1>& Ex) {
    for (unsigned k = 0, ikz = Ex.nz - n_cpml_zp; k < n_cpml_zp; k++, ikz++) {
        for (unsigned i = 0; i < Psi_hyz_zp.nx; i++) {
            for (unsigned j = 0; j < Psi_hyz_zp.ny; j++) {
                Psi_hyz_zp.p[i][j][k] = Psi_hyz_zp.p[i][j][k] * cpml_b_mz_zp.p[k] + cpml_a_mz_zp.p[k]*(Ex.p[i][j][ikz] - Ex.p[i][j][ikz - 1]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_hzy_yp(const data3d<type1>& Ex) {
    for (unsigned j = 0, jhy = Ex.ny - n_cpml_yp; j < n_cpml_yp; j++, jhy++) {
        for (unsigned i = 0; i < Psi_hzy_yp.nx; i++) {
            for (unsigned k = 0; k < Psi_hzy_yp.nz; k++) {
                Psi_hzy_yp.p[i][j][k] = Psi_hzy_yp.p[i][j][k] * cpml_b_my_yp.p[j] + cpml_a_my_yp.p[j]*(Ex.p[i][jhy][k] - Ex.p[i][jhy - 1 ][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_hxz_zp(const data3d<type1>& Ey) {
    for (unsigned k = 0, khz = Ey.nz - n_cpml_zp; k < n_cpml_zp; k++, khz++) {
        for (unsigned j = 0; j < Psi_hxz_zp.ny; j++) {
            for (unsigned i = 0; i < Psi_hxz_zp.nx; i++) {
                Psi_hxz_zp.p[i][j][k] = Psi_hxz_zp.p[i][j][k] * cpml_b_mz_zp.p[k] + cpml_a_mz_zp.p[k]*(Ey.p[i][j][khz] - Ey.p[i][j][khz - 1]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_hzx_xp(const data3d<type1>& Ey) {
    for (unsigned i = 0, ihx = Ey.nx - n_cpml_xp; i < n_cpml_xp; i++, ihx++) {
        for (unsigned j = 0; j < Psi_hzx_xp.ny; j++) {
            for (unsigned k = 0; k < Psi_hzx_xp.nz; k++) {
                Psi_hzx_xp.p[i][j][k] = Psi_hzx_xp.p[i][j][k] * cpml_b_mx_xp.p[i] + cpml_a_mx_xp.p[i]*(Ey.p[ihx][j][k] - Ey.p[ihx - 1][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_hyx_xp(const data3d<type1>& Ez) {
    for (unsigned i = 0, imx = Ez.nx - n_cpml_zp; i < n_cpml_zp; i++, imx++) {
        for (unsigned j = 0; j < Psi_hyx_xp.ny; j++) {
            for (unsigned k = 0; k < Psi_hyx_xp.nz; k++) {
                Psi_hyx_xp.p[i][j][k] = Psi_hyx_xp.p[i][j][k] * cpml_b_mx_xp.p[i] + cpml_a_mx_xp.p[i]*(Ez.p[imx][j][k] - Ez.p[imx - 1 ][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_hxy_yp(const data3d<type1>& Ez) {
    for (unsigned j = 0, jmy = Ez.ny - n_cpml_yp; j < n_cpml_yp; j++, jmy++) {
        for (unsigned i = 0; i < Psi_hxy_yp.nx; i++) {
            for (unsigned k = 0; k < Psi_hxy_yp.nz; k++) {
                Psi_hxy_yp.p[i][j][k] = Psi_hxy_yp.p[i][j][k] * cpml_b_my_yp.p[j] + cpml_a_my_yp.p[j]*(Ez.p[i][jmy][k] - Ez.p[i][jmy - 1][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_eyz_zn(const data3d<type1>& Hx) {
    for (unsigned k = 0, kplus = 1; k < n_cpml_zn; k++, kplus++) {
        for (unsigned i = 0; i < Psi_eyz_zn.nx; i++) {
            for (unsigned j = 0; Psi_eyz_zn.ny; j++) {
                Psi_eyz_zn.p[i][j][k] = cpml_b_ez_zn.p[k] * Psi_eyz_zn.p[i][j][k] + cpml_a_ez_zn.p[k]*(Hx.p[i][j][kplus] - Hx.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_ezy_yn(const data3d<type1>& Hx) {
    for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
        for (unsigned i = 0; i < Psi_ezy_yn.nx; i++) {
            for (unsigned k = 0; k < Psi_ezy_yn.nz; k++) {
                Psi_ezy_yn.p[i][j][k] = cpml_b_ey_yn.p[j] * Psi_ezy_yn.p[i][j][k] + cpml_a_ey_yn.p[j]*(Hx.p[i][jplus][k] - Hx.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_exz_zn(const data3d<type1>& Hy) {
    for (unsigned k = 0, kplus = 1; k < n_cpml_zn; k++, kplus++) {
        for (unsigned i = 0; i < Psi_exz_zn.nx; i++) {
            for (unsigned j = 0; j < Psi_exz_zn.ny; j++) {
                Psi_exz_zn.p[i][j][k] = cpml_b_ez_zn.p[k] * Psi_exz_zn.p[i][j][k] + cpml_a_ez_zn.p[k]*(Hy.p[i][j][kplus] - Hy.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_ezx_xn(const data3d<type1>& Hy) {
    for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
        for (unsigned j = 0; j < Psi_ezx_xn.ny; j++) {
            for (unsigned k = 0; k < Psi_ezx_xn.nz; k++) {
                Psi_ezx_xn.p[i][j][k] = cpml_b_ex_xn.p[i] * Psi_ezx_xn.p[i][j][k] + cpml_a_ex_xn.p[i]*(Hy.p[iplus][j][k] - Hy.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_eyx_xn(const data3d<type1>& Hz) {
    for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
        for (unsigned j = 0; j < Psi_eyx_xn.ny; j++) {
            for (unsigned k = 0; k < Psi_eyx_xn.nz; k++) {
                Psi_eyx_xn.p[i][j][k] = cpml_b_ex_xn.p[i] * Psi_eyx_xn.p[i][j][k] + cpml_a_ex_xn.p[i]*(Hz.p[iplus][j][k] - Hz.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_exy_yn(const data3d<type1>& Hz) {
    for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
        for (unsigned k = 0; k < Psi_exy_yn.nz; k++) {
            for (unsigned i = 0; i < Psi_exy_yn.ny; i++) {
                Psi_exy_yn.p[i][j][k] = cpml_b_ey_yn.p[j] * Psi_exy_yn.p[i][j][k] + cpml_a_ey_yn.p[j]*(Hz.p[i][jplus][k] - Hz.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_hyz_zn(const data3d<type1>& Ex) {
    for (unsigned k = 0, kplus = 1; k < n_cpml_zn; k++, kplus++) {
        for (unsigned j = 0; j < Psi_hyz_zn.ny; j++) {
            for (unsigned i = 0; i < Psi_hyz_zn.nx; i++) {
                Psi_hyz_zn.p[i][j][k] = cpml_b_mz_zn.p[k] * Psi_hyz_zn.p[i][j][k] + cpml_a_mz_zn.p[k]*(Ex.p[i][j][kplus] - Ex.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_hzy_yn(const data3d<type1>& Ex) {
    for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
        for (unsigned i = 0; i < Psi_hzy_yn.nx; i++) {
            for (unsigned k = 0; k < Psi_hzy_yn.nz; k++) {
                Psi_hzy_yn.p[i][j][k] = cpml_b_my_yn.p[j] * Psi_hzy_yn.p[i][j][k] + cpml_a_my_yn.p[j]*(Ex.p[i][jplus][k] - Ex.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_hxz_zn(const data3d<type1>& Ey) {
    for (unsigned k = 0, kplus = 1; k < n_cpml_zn; k++, kplus++) {
        for (unsigned i = 0; i < Psi_hxz_zn.nx; i++) {
            for (unsigned j = 0; j < Psi_hxz_zn.ny; j++) {
                Psi_hxz_zn.p[i][j][k] = cpml_b_mz_zn.p[k] * Psi_hxz_zn.p[i][j][k] + cpml_a_mz_zn.p[k]*(Ey.p[i][j][kplus] - Ey.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_hzx_xn(const data3d<type1>& Ey) {
    for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
        for (unsigned j = 0; j < Psi_hzx_xn.ny; j++) {
            for (unsigned k = 0; k < Psi_hzx_xn.nz; k++) {
                Psi_hzx_xn.p[i][j][k] = Psi_hzx_xn.p[i][j][k] * cpml_b_mx_xn.p[i] + cpml_a_mx_xn.p[i]*(Ey.p[iplus][j][k] - Ey.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_hyx_xn(const data3d<type1>& Ez) {
    for (unsigned i = 0, iplus = 1; i < n_cpml_xn; i++, iplus++) {
        for (unsigned j = 0; j < Psi_hyx_xn.ny; j++) {
            for (unsigned k = 0; k < Psi_hyx_xn.nz; k++) {
                Psi_hyx_xn.p[i][j][k] = Psi_hyx_xn.p[i][j][k] * cpml_b_mx_xn.p[i] + cpml_a_mx_xn.p[i]*(Ez.p[iplus][j][k] - Ez.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml<type1>::updatePsi_hxy_yn(const data3d<type1>& Ez) {
    for (unsigned j = 0, jplus = 1; j < n_cpml_yn; j++, jplus++) {
        for (unsigned i = 0; i < Psi_hxy_yn.nx; i++) {
            for (unsigned k = 0; k < Psi_hxy_yn.nz; k++) {
                Psi_hxy_yn.p[i][j][k] = Psi_hxy_yn.p[i][j][k] * cpml_b_my_yn.p[j] + cpml_a_my_yn.p[j]*(Ez.p[i][jplus][k] - Ez.p[i][j][k]);
            }
        }
    }
}

//#include "cpml.cpp"
#endif	/* CPML_H */


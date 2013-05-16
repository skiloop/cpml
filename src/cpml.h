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
// some common constants
const double C = 2.99792458E8; // speed of light in free space
const double me = 9.110e-31; // electricity mass
const double e = 1.602e-19; // electricity charge
const double mu_0 = 4.0 * M_PI * 1.0E-7;
const double eps_0 = 1.0 / (C * C * mu_0);
const double M_PI_TWO = M_PI * 2;
const double Mu0DivEps0 = mu_0 / eps_0;

/***************************************************************
 * FDTD region 
 * Field    |    x     |    y    |   z
 * ----------------------------------------------
 * Ez       |   I      |    J    |   K
 * Ex       |   I-1    |    J    |   K-1
 * Ey       |   I      |   J-1   |   K-1
 * Hx       |   I      |   J-1   |   K
 * Hy       |   I-1    |    J    |   K
 * Hz       |   I-1    |   J-1   |   K-1
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
     * @param Hz1
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

    //=======================================================
    // private functions
    //=======================================================
private:
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
     * @param imax
     * @param jmax
     * @param kmax
     */
    void createCPMLArrays(unsigned imax, unsigned jmax, unsigned kmax);
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
    
    void updatePsi_eyz_n(const data3d<type1>& Hz);
    void updatePsi_ezy_n(const data3d<type1>& Hy);
    void updatePsi_exz_n(const data3d<type1>& Hz);
    void updatePsi_ezx_n(const data3d<type1>& Hx);
    void updatePsi_eyx_n(const data3d<type1>& Hx);
    void updatePsi_exy_n(const data3d<type1>& Hy);
    
    void updatePsi_hyz_n(const data3d<type1>& Ez);
    void updatePsi_hzy_n(const data3d<type1>& Ey);
    void updatePsi_hxz_n(const data3d<type1>& Ez);
    void updatePsi_hzx_n(const data3d<type1>& Ex);
    void updatePsi_hyx_n(const data3d<type1>& Ex);
    void updatePsi_hxy_n(const data3d<type1>& Ey);
    
        void updatePsi_eyz_p(const data3d<type1>& Hz);
    void updatePsi_ezy_p(const data3d<type1>& Hy);
    void updatePsi_exz_p(const data3d<type1>& Hz);
    void updatePsi_ezx_p(const data3d<type1>& Hx);
    void updatePsi_eyx_p(const data3d<type1>& Hx);
    void updatePsi_exy_p(const data3d<type1>& Hy);
    
    void updatePsi_hyz_p(const data3d<type1>& Ez);
    void updatePsi_hzy_p(const data3d<type1>& Ey);
    void updatePsi_hxz_p(const data3d<type1>& Ez);
    void updatePsi_hzx_p(const data3d<type1>& Ex);
    void updatePsi_hyx_p(const data3d<type1>& Ex);
    void updatePsi_hxy_p(const data3d<type1>& Ey);
    
    void updateEFieldCPML_x(data3d<type1>&Ey,data3d<type1>&Ez);
    void updateEFieldCPML_y(data3d<type1>&Ex,data3d<type1>&Ez);
    void updateEFieldCPML_z(data3d<type1>&Ex,data3d<type1>&Ey);
    
    void updateMFieldCPML_x(data3d<type1>&Hy,data3d<type1>&Hz);
    void updateMFieldCPML_y(data3d<type1>&Hx,data3d<type1>&Hz);
    void updateMFieldCPML_z(data3d<type1>&Hx,data3d<type1>&Hy);
};
#include "cpml.cpp"
#endif	/* CPML_H */


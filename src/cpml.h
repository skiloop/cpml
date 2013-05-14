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

template<class type1, class type2>
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

    //=======================================================
    // private functions
    //=======================================================
private:
    void setCPMLRegion(short pmlWidth);
    void setCPMLRegion(short width_xn, short width_xp, short width_yn, short width_yp, short width_zn, short width_zp);
    void createCPMLArrays(unsigned imax, unsigned jmax, unsigned kmax);
    void initCoefficientArrays(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dx, type1 dy, type1 dz);
    void initCoefficientArraysXN(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dx);
    void initCoefficientArraysXP(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dx);
    void initCoefficientArraysYN(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dy);
    void initCoefficientArraysYP(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dy);
    void initCoefficientArraysZN(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dz);
    void initCoefficientArraysZP(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dz);
};
#include "cpml.cpp"
#endif	/* CPML_H */


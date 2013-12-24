/* 
* File:   cpml.h
* Author: skiloop
*
* Created on April 3, 2013, 8:41 AM
*/

#ifndef CPML_H
#define	CPML_H

#include "data1d.h"
#include "data3d.h"
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
***************************************************************/
template<class T>
class cpml {
public:
	/**
	* default constructor
	*/
	cpml();

	/**
	* constructor
	* @param width_xn width of cpml layer on x negative side
	* @param width_xp width of cpml layer on x positive side
	* @param width_yn width of cpml layer on y negative side
	* @param width_yp width of cpml layer on y positive side
	* @param width_zn width of cpml layer on z negative side
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
	void initCoefficientArrays(short pmlOrder, T alphaOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dx, T dy, T dz,
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
	void initCoefficientArraysXN(short pmlOrder, T alphaOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dx, data3d<T>&Ceyhz, data3d<T>&Cezhy, data3d<T>&Chyez, data3d<T>&Chzey);

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
	void initCoefficientArraysXP(short pmlOrder, T alphaOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dx, data3d<T>&Ceyhz, data3d<T>&Cezhy, data3d<T>&Chyez, data3d<T>&Chzey);

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
	void initCoefficientArraysYN(short pmlOrder, T alphaOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dy, data3d<T>&Cexhz, data3d<T>&Cezhx, data3d<T>&Chxez, data3d<T>&Chzex);
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
	void initCoefficientArraysYP(short pmlOrder, T alphaOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dy, data3d<T>&Cexhz, data3d<T>&Cezhx, data3d<T>&Chxez, data3d<T>&Chzex);
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
	void initCoefficientArraysZN(short pmlOrder, T alphaOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dz, data3d<T>&Ceyhx, data3d<T>&Cexhy, data3d<T>&Chyex, data3d<T>&Chxey);
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
	void initCoefficientArraysZP(short pmlOrder, T alphaOrder, T sigmaMax, T kappaMax, T alphaMax, T epsR, T dt, T dz, data3d<T>&Ceyhx, data3d<T>&Cexhy, data3d<T>&Chyex, data3d<T>&Chxey);

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
#include "cpml.hpp"

#endif	/* CPML_H */


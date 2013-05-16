/* 
 * File:   cpml.cpp
 * Author: skiloop
 * 
 * Created on 2013年5月14日, 上午9:58
 */

#include "cpml.h"


//#include "cpml.h"

template<class type1>
cpml::cpml()
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
cpml::cpml(const cpml& orig)
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
cpml::~cpml() {
}

template<class type1>
void cpml::createCPMLArrays(unsigned imax, unsigned jmax, unsigned kmax) {
    //============================================
    // cpml arrays
    //============================================
    // xn arrays
    if (is_cpml_xn) {
        // x direction
        cpml_a_ex_xn.CreateStruct(n_cpml_xn);
        cpml_b_ex_xn.CreateStruct(n_cpml_xn);
        cpml_a_mx_xn.CreateStruct(n_cpml_xn);
        cpml_b_mx_xn.CreateStruct(n_cpml_xn);
        Psi_eyx_xn.CreateStruct(n_cpml_xn, jmax, kmax);
        Psi_ezx_xn.CreateStruct(n_cpml_xn, jmax, kmax);
        Psi_hyx_xn.CreateStruct(n_cpml_xn, jmax, kmax);
        Psi_hzx_xn.CreateStruct(n_cpml_xn, jmax, kmax);
        CPsi_eyx_xn.CreateStruct(n_cpml_xn, jmax, kmax);
        CPsi_ezx_xn.CreateStruct(n_cpml_xn, jmax, kmax);
        CPsi_hyx_xn.CreateStruct(n_cpml_xn, jmax, kmax);
        CPsi_hzx_xn.CreateStruct(n_cpml_xn, jmax, kmax);
    }
    // xp arrays
    if (is_cpml_xp) {
        cpml_a_ex_xp.CreateStruct(n_cpml_xp);
        cpml_b_ex_xp.CreateStruct(n_cpml_xp);
        cpml_a_mx_xp.CreateStruct(n_cpml_xp);
        cpml_b_mx_xp.CreateStruct(n_cpml_xp);
        Psi_eyx_xp.CreateStruct(n_cpml_xp, jmax, kmax);
        Psi_ezx_xp.CreateStruct(n_cpml_xp, jmax, kmax);
        Psi_hyx_xp.CreateStruct(n_cpml_xp, jmax, kmax);
        Psi_hzx_xp.CreateStruct(n_cpml_xp, jmax, kmax);
        CPsi_eyx_xp.CreateStruct(n_cpml_xp, jmax, kmax);
        CPsi_ezx_xp.CreateStruct(n_cpml_xp, jmax, kmax);
        CPsi_hyx_xp.CreateStruct(n_cpml_xp, jmax, kmax);
        CPsi_hzx_xp.CreateStruct(n_cpml_xp, jmax, kmax);
    }
    // yn arrays
    if (is_cpml_yn) {
        // y direction
        cpml_a_ey_yn.CreateStruct(n_cpml_yn);
        cpml_b_ey_yn.CreateStruct(n_cpml_yn);
        cpml_a_my_yn.CreateStruct(n_cpml_yn);
        cpml_b_my_yn.CreateStruct(n_cpml_yn);
        Psi_exy_yn.CreateStruct(imax, n_cpml_yn, kmax);
        Psi_ezy_yn.CreateStruct(imax, n_cpml_yn, kmax);
        Psi_hxy_yn.CreateStruct(imax, n_cpml_yn, kmax);
        Psi_hzy_yn.CreateStruct(imax, n_cpml_yn, kmax);
        CPsi_exy_yn.CreateStruct(imax, n_cpml_yn, kmax);
        CPsi_ezy_yn.CreateStruct(imax, n_cpml_yn, kmax);
        CPsi_hxy_yn.CreateStruct(imax, n_cpml_yn, kmax);
        CPsi_hzy_yn.CreateStruct(imax, n_cpml_yn, kmax);
    }
    // yp arrays
    if (is_cpml_yp) {
        cpml_a_ey_yp.CreateStruct(n_cpml_yp);
        cpml_b_ey_yp.CreateStruct(n_cpml_yp);
        cpml_a_my_yp.CreateStruct(n_cpml_yp);
        cpml_b_my_yp.CreateStruct(n_cpml_yp);
        Psi_exy_yp.CreateStruct(imax, n_cpml_yp, kmax);
        Psi_ezy_yp.CreateStruct(imax, n_cpml_yp, kmax);
        Psi_hxy_yp.CreateStruct(imax, n_cpml_yp, kmax);
        Psi_hzy_yp.CreateStruct(imax, n_cpml_yp, kmax);
        CPsi_exy_yp.CreateStruct(imax, n_cpml_yp, kmax);
        CPsi_ezy_yp.CreateStruct(imax, n_cpml_yp, kmax);
        CPsi_hxy_yp.CreateStruct(imax, n_cpml_yp, kmax);
        CPsi_hzy_yp.CreateStruct(imax, n_cpml_yp, kmax);
    }
    // zn arrays
    if (is_cpml_zn) {
        // z direction
        cpml_a_ez_zn.CreateStruct(n_cpml_zn);
        cpml_b_ez_zn.CreateStruct(n_cpml_zn);
        cpml_a_mz_zn.CreateStruct(n_cpml_zn);
        cpml_b_mz_zn.CreateStruct(n_cpml_zn);
        Psi_exz_zn.CreateStruct(imax, jmax, n_cpml_zn);
        Psi_eyz_zn.CreateStruct(imax, jmax, n_cpml_zn);
        Psi_hxz_zn.CreateStruct(imax, jmax, n_cpml_zn);
        Psi_hyz_zn.CreateStruct(imax, jmax, n_cpml_zn);
        CPsi_exz_zn.CreateStruct(imax, jmax, n_cpml_zn);
        CPsi_eyz_zn.CreateStruct(imax, jmax, n_cpml_zn);
        CPsi_hxz_zn.CreateStruct(imax, jmax, n_cpml_zn);
        CPsi_hyz_zn.CreateStruct(imax, jmax, n_cpml_zn);
    }
    // zp arrays
    if (is_cpml_zp) {
        cpml_a_ez_zp.CreateStruct(n_cpml_zp);
        cpml_b_ez_zp.CreateStruct(n_cpml_zp);
        cpml_a_mz_zp.CreateStruct(n_cpml_zp);
        cpml_b_mz_zp.CreateStruct(n_cpml_zp);
        Psi_exz_zp.CreateStruct(imax, jmax, n_cpml_zp);
        Psi_eyz_zp.CreateStruct(imax, jmax, n_cpml_zp);
        Psi_hxz_zp.CreateStruct(imax, jmax, n_cpml_zp);
        Psi_hyz_zp.CreateStruct(imax, jmax, n_cpml_zp);
        CPsi_exz_zp.CreateStruct(imax, jmax, n_cpml_zp);
        CPsi_eyz_zp.CreateStruct(imax, jmax, n_cpml_zp);
        CPsi_hxz_zp.CreateStruct(imax, jmax, n_cpml_zp);
        CPsi_hyz_zp.CreateStruct(imax, jmax, n_cpml_zp);
    }
}

template<class type1>
void cpml::initCoefficientArrays(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dx, type1 dy, type1 dz,
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
void cpml::initCoefficientArraysXN(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dx,
        data3d<type1>&Ceyhz, data3d<type1>&Cezhy, data3d<type1>&Chyez, data3d<type1>&Chzey) {
    if (is_cpml_xn) {
        type1 sigmaMax = sigmaRatio * (pmlOrder + 1) / (150 * M_PI * dx);

        for (unsigned i = 0; i < n_cpml_xn; i++) {
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
                    // TODO fix region size
                    CPsi_eyx_xn.p[i][j][k] = Ceyhz.p[i + 1][j][k] * dx;
                    CPsi_ezx_xn.p[i][j][k] = Cezhy.p[i + 1][j][k] * dx;
                    Ceyhz.p[i + 1][j][k] = Ceyhz.p[i + 1][j][k] / kappa_ex;
                    Cezhy.p[i + 1][j][k] = Cezhy.p[i + 1][j][k] / kappa_ex;
                }
            }
            for (unsigned j = 0; j < Psi_hyx_xn.ny; j++) {
                for (unsigned k = 0; k < Psi_hyx_xn.nz; k++) {
                    CPsi_hyx_xn.p[i][j][k] = Chyez.p[i][j][k] * dx;
                    CPsi_hzx_xn.p[i][j][k] = Chzey.p[i][j][k] * dx;
                    Chyez.p[i][j][k] = Chyez.p[i][j][k] / kappa_mx;
                    Chzey.p[i][j][k] = Chzey.p[i][j][k] / kappa_mx;
                }
            }
        }
    }
}

template<class type1>
void cpml::initCoefficientArraysXP(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dx,
        data3d<type1>&Ceyhz, data3d<type1>&Cezhy, data3d<type1>&Chyez, data3d<type1>&Chzey) {
    if (is_cpml_xp) {
        type1 sigmaMax = sigmaRatio * (pmlOrder + 1) / (150 * M_PI * dx);
        unsigned iex = Psi_eyx_xp.nx - n_cpml_xp - 1;
        unsigned ihx = Psi_hyx_xp.nx - n_cpml_xp - 1;
        for (unsigned i = 0; i < n_cpml_xp; i++) {
            type1 rho_e = (n_cpml_xp - i - 0.75) / n_cpml_xp;
            type1 rho_m = (n_cpml_xp - i - 0.25) / n_cpml_xp;
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
                    // TODO fix region size
                    CPsi_eyx_xp.p[i][j][k] = Ceyhz.p[iex][j][k] * dx;
                    CPsi_ezx_xp.p[i][j][k] = Cezhy.p[iex][j][k] * dx;
                    Ceyhz.p[iex][j][k] = Ceyhz.p[iex][j][k] / kappa_ex;
                    Cezhy.p[iex][j][k] = Cezhy.p[iex][j][k] / kappa_ex;
                }
            }
            for (unsigned j = 0; j < Psi_hyx_xp.ny; j++) {
                for (unsigned k = 0; k < Psi_hyx_xp.nz; k++) {
                    CPsi_hyx_xp.p[i][j][k] = Chyez.p[ihx][j][k] * dx;
                    CPsi_hzx_xp.p[i][j][k] = Chzey.p[ihx][j][k] * dx;
                    Chyez.p[ihx][j][k] = Chyez.p[ihx][j][k] / kappa_mx;
                    Chzey.p[ihx][j][k] = Chzey.p[ihx][j][k] / kappa_mx;
                }
            }
            ihx++;
            iex++;
        }
    }
}

template<class type1>
void cpml::initCoefficientArraysYN(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dy,
        data3d<type1>&Cexhz, data3d<type1>&Cezhx, data3d<type1>&Chxez, data3d<type1>&Chzex) {
    if (is_cpml_yn) {
        type1 sigmaMax = sigmaRatio * (pmlOrder + 1) / (150 * M_PI * dy);
        for (unsigned i = 0; i < n_cpml_yn; i++) {
            type1 rho_e = (n_cpml_yn - i - 0.75) / n_cpml_yn;
            type1 rho_m = (n_cpml_yn - i - 0.25) / n_cpml_yn;
            type1 rho_e_pmlOrder = pow(rho_e, pmlOrder);
            type1 rho_m_pmlOrder = pow(rho_m, pmlOrder);
            type1 sigma_pey = sigmaMax*rho_e_pmlOrder;
            type1 sigma_pmy = sigmaMax * rho_m_pmlOrder*Mu0DivEps0;
            type1 kappa_ey = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            type1 kappa_my = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            type1 alpha_ey = alphaMax*rho_e_pmlOrder;
            type1 alpha_my = alphaMax*rho_m_pmlOrder;
            cpml_b_ey_yn.p[i] = exp((-dt / eps_0) * sigma_pey / kappa_ey + alpha_ey);
            cpml_b_my_yn.p[i] = exp((-dt / mu_0) * sigma_pmy / kappa_my + alpha_my);
            cpml_a_ey_yn.p[i] = 1 / dy * (cpml_b_ey_yn.p[i] - 1.0) * sigma_pey / (kappa_ey * (sigma_pey + kappa_ey * alpha_ey));
            cpml_a_my_yn.p[i] = 1 / dy * (cpml_b_my_yn.p[i] - 1.0) * sigma_pmy / (kappa_my * (sigma_pmy + kappa_my * alpha_my));

            for (unsigned j = 0; j < Psi_exy_yn.ny; j++) {
                for (unsigned k = 0; k < Psi_exy_yn.nz; k++) {
                    // TODO fix region size
                    CPsi_exy_yn.p[i][j][k] = Cexhz.p[i + 1][j][k] * dy;
                    CPsi_ezy_yn.p[i][j][k] = Cezhx.p[i + 1][j][k] * dy;
                    Cexhz.p[i + 1][j][k] = Cexhz.p[i + 1][j][k] / kappa_ey;
                    Cezhx.p[i + 1][j][k] = Cezhx.p[i + 1][j][k] / kappa_ey;
                }
            }
            for (unsigned j = 0; j < Psi_hxy_yn.ny; j++) {
                for (unsigned k = 0; k < Psi_hxy_yn.nz; k++) {
                    CPsi_hxy_yn.p[i][j][k] = Chxez.p[i][j][k] * dy;
                    CPsi_hzy_yn.p[i][j][k] = Chzex.p[i][j][k] * dy;
                    Chxez.p[i][j][k] = Chxez.p[i][j][k] / kappa_my;
                    Chzex.p[i][j][k] = Chzex.p[i][j][k] / kappa_my;
                }
            }
        }
    }
}

template<class type1>
void cpml::initCoefficientArraysYP(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dy,
        data3d<type1>&Cexhz, data3d<type1>&Cezhx, data3d<type1>&Chxez, data3d<type1>&Chzex) {
    if (is_cpml_yp) {
        type1 sigmaMax = sigmaRatio * (pmlOrder + 1) / (150 * M_PI * dy);
        unsigned iex = Psi_exy_yp.nx - n_cpml_yp - 1;
        unsigned ihx = Psi_hxy_yp.nx - n_cpml_yp - 1;
        for (unsigned i = 0; i < n_cpml_yp; i++) {
            type1 rho_e = (n_cpml_yp - i - 0.75) / n_cpml_yp;
            type1 rho_m = (n_cpml_yp - i - 0.25) / n_cpml_yp;
            type1 rho_e_pmlOrder = pow(rho_e, pmlOrder);
            type1 rho_m_pmlOrder = pow(rho_m, pmlOrder);
            type1 sigma_pey = sigmaMax*rho_e_pmlOrder;
            type1 sigma_pmy = sigmaMax * rho_m_pmlOrder*Mu0DivEps0;
            type1 kappa_ey = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            type1 kappa_my = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            type1 alpha_ey = alphaMax*rho_e_pmlOrder;
            type1 alpha_my = alphaMax*rho_m_pmlOrder;
            cpml_b_ey_yp.p[i] = exp((-dt / eps_0) * sigma_pey / kappa_ey + alpha_ey);
            cpml_b_my_yp.p[i] = exp((-dt / mu_0) * sigma_pmy / kappa_my + alpha_my);
            cpml_a_ey_yp.p[i] = 1 / dy * (cpml_b_ey_yp.p[i] - 1.0) * sigma_pey / (kappa_ey * (sigma_pey + kappa_ey * alpha_ey));
            cpml_a_my_yp.p[i] = 1 / dy * (cpml_b_my_yp.p[i] - 1.0) * sigma_pmy / (kappa_my * (sigma_pmy + kappa_my * alpha_my));

            for (unsigned j = 0; j < Psi_exy_yp.ny; j++) {
                for (unsigned k = 0; k < Psi_exy_yp.nz; k++) {
                    // TODO fix region size
                    CPsi_exy_yp.p[i][j][k] = Cexhz.p[iex][j][k] * dy;
                    CPsi_ezy_yp.p[i][j][k] = Cezhx.p[iex][j][k] * dy;
                    Cexhz.p[iex][j][k] = Cexhz.p[iex][j][k] / kappa_ey;
                    Cezhx.p[iex][j][k] = Cezhx.p[iex][j][k] / kappa_ey;
                }
            }
            for (unsigned j = 0; j < Psi_hxy_yp.ny; j++) {
                for (unsigned k = 0; k < Psi_hxy_yp.nz; k++) {
                    CPsi_hxy_yp.p[i][j][k] = Chxez.p[ihx][j][k] * dy;
                    CPsi_hzy_yp.p[i][j][k] = Chzex.p[ihx][j][k] * dy;
                    Chxez.p[ihx][j][k] = Chxez.p[ihx][j][k] / kappa_my;
                    Chzex.p[ihx][j][k] = Chzex.p[ihx][j][k] / kappa_my;
                }
            }
            ihx++;
            iex++;
        }
    }
}

template<class type1>
void cpml::initCoefficientArraysZN(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dz,
        data3d<type1>&Ceyhx, data3d<type1>&Cexhy, data3d<type1>&Chyex, data3d<type1>&Chxey) {
    if (is_cpml_zn) {
        type1 sigmaMax = sigmaRatio * (pmlOrder + 1) / (150 * M_PI * dz);
        for (unsigned i = 0; i < n_cpml_zn; i++) {
            type1 rho_e = (n_cpml_zn - i - 0.75) / n_cpml_zn;
            type1 rho_m = (n_cpml_zn - i - 0.25) / n_cpml_zn;
            type1 rho_e_pmlOrder = pow(rho_e, pmlOrder);
            type1 rho_m_pmlOrder = pow(rho_m, pmlOrder);
            type1 sigma_pez = sigmaMax*rho_e_pmlOrder;
            type1 sigma_pmz = sigmaMax * rho_m_pmlOrder*Mu0DivEps0;
            type1 kappa_ez = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            type1 kappa_mz = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            type1 alpha_ez = alphaMax*rho_e_pmlOrder;
            type1 alpha_mz = alphaMax*rho_m_pmlOrder;
            cpml_b_ez_zn.p[i] = exp((-dt / eps_0) * sigma_pez / kappa_ez + alpha_ez);
            cpml_b_mz_zn.p[i] = exp((-dt / mu_0) * sigma_pmz / kappa_mz + alpha_mz);
            cpml_a_ez_zn.p[i] = 1 / dz * (cpml_b_ez_zn.p[i] - 1.0) * sigma_pez / (kappa_ez * (sigma_pez + kappa_ez * alpha_ez));
            cpml_a_mz_zn.p[i] = 1 / dz * (cpml_b_mz_zn.p[i] - 1.0) * sigma_pmz / (kappa_mz * (sigma_pmz + kappa_mz * alpha_mz));

            for (unsigned j = 0; j < Psi_exz_zn.ny; j++) {
                for (unsigned k = 0; k < Psi_exz_zn.nz; k++) {
                    // TODO fix region size
                    CPsi_exz_zn.p[i][j][k] = Cexhy.p[i + 1][j][k] * dz;
                    CPsi_eyz_zn.p[i][j][k] = Ceyhx.p[i + 1][j][k] * dz;
                    Cexhy.p[i + 1][j][k] = Cexhy.p[i + 1][j][k] / kappa_ez;
                    Ceyhx.p[i + 1][j][k] = Ceyhx.p[i + 1][j][k] / kappa_ez;
                }
            }
            for (unsigned j = 0; j < Psi_hxz_zn.ny; j++) {
                for (unsigned k = 0; k < Psi_hxz_zn.nz; k++) {
                    CPsi_hxz_zn.p[i][j][k] = Chxey.p[i][j][k] * dz;
                    CPsi_hyz_zn.p[i][j][k] = Chyex.p[i][j][k] * dz;
                    Chxey.p[i][j][k] = Chxey.p[i][j][k] / kappa_mz;
                    Chyex.p[i][j][k] = Chyex.p[i][j][k] / kappa_mz;
                }
            }
        }
    }
}

template<class type1>
void cpml::initCoefficientArraysZP(short pmlOrder, type1 sigmaRatio, type1 kappaMax, type1 alphaMax, type1 dt, type1 dz,
        data3d<type1>&Ceyhx, data3d<type1>&Cexhy, data3d<type1>&Chyex, data3d<type1>&Chxey) {
    if (is_cpml_zp) {
        type1 sigmaMax = sigmaRatio * (pmlOrder + 1) / (150 * M_PI * dz);
        unsigned iez = Psi_eyz_zp.nx - n_cpml_zp - 1;
        unsigned ihz = Psi_hyz_zp.nx - n_cpml_zp - 1;
        for (unsigned i = 0; i < n_cpml_zp; i++) {
            type1 rho_e = (n_cpml_zp - i - 0.75) / n_cpml_zp;
            type1 rho_m = (n_cpml_zp - i - 0.25) / n_cpml_zp;
            type1 rho_e_pmlOrder = pow(rho_e, pmlOrder);
            type1 rho_m_pmlOrder = pow(rho_m, pmlOrder);
            type1 sigma_pez = sigmaMax*rho_e_pmlOrder;
            type1 sigma_pmz = sigmaMax * rho_m_pmlOrder*Mu0DivEps0;
            type1 kappa_ez = 1 + (kappaMax - 1) * rho_e_pmlOrder;
            type1 kappa_mz = 1 + (kappaMax - 1) * rho_m_pmlOrder;
            type1 alpha_ez = alphaMax*rho_e_pmlOrder;
            type1 alpha_mz = alphaMax*rho_m_pmlOrder;
            cpml_b_ez_zp.p[i] = exp((-dt / eps_0) * sigma_pez / kappa_ez + alpha_ez);
            cpml_b_mz_zp.p[i] = exp((-dt / mu_0) * sigma_pmz / kappa_mz + alpha_mz);
            cpml_a_ez_zp.p[i] = 1 / dz * (cpml_b_ez_zp.p[i] - 1.0) * sigma_pez / (kappa_ez * (sigma_pez + kappa_ez * alpha_ez));
            cpml_a_mz_zp.p[i] = 1 / dz * (cpml_b_mz_zp.p[i] - 1.0) * sigma_pmz / (kappa_mz * (sigma_pmz + kappa_mz * alpha_mz));

            for (unsigned j = 0; j < Psi_eyz_zp.ny; j++) {
                for (unsigned k = 0; k < Psi_eyz_zp.nz; k++) {
                    // TODO fix region size
                    CPsi_eyz_zp.p[i][j][k] = Ceyhx.p[iez][j][k] * dz;
                    CPsi_exz_zp.p[i][j][k] = Cexhy.p[iez][j][k] * dz;
                    Ceyhx.p[iez][j][k] = Ceyhx.p[iez][j][k] / kappa_ez;
                    Cexhy.p[iez][j][k] = Cexhy.p[iez][j][k] / kappa_ez;
                }
            }
            for (unsigned j = 0; j < Psi_hyz_zp.ny; j++) {
                for (unsigned k = 0; k < Psi_hyz_zp.nz; k++) {
                    CPsi_hyz_zp.p[i][j][k] = Chyex.p[ihz][j][k] * dz;
                    CPsi_hxz_zp.p[i][j][k] = Chxey.p[ihz][j][k] * dz;
                    Chyex.p[ihz][j][k] = Chyex.p[ihz][j][k] / kappa_mz;
                    Chxey.p[ihz][j][k] = Chxey.p[ihz][j][k] / kappa_mz;
                }
            }
            ihz++;
            iez++;
        }
    }
}

template<class type1>
void cpml::updateCPML_E_Fields(data3d<type1>& Ex, data3d<type1>& Ey, data3d<type1>& Ez,
        const data3d<type1>& Hx, const data3d<type1>& Hy, const data3d<type1>& Hz) {
    updatePsiForEFields(Hx, Hy, Hz);
    updateCPML_x(Ex);
    updateCPML_y(Ey);
    updateCPML_z(Ez);
}

template<class type1>
void cpml::updateCPML_M_Fields(data3d<type1>& Hx, data3d<type1>& Hy, data3d<type1>& Hz,
        const data3d<type1>& Ex, const data3d<type1>& Ey, const data3d<type1>& Ez) {
    updatePsiForMFields(Ex, Ey, Ez);
    updateCPML_x(Hx);
    updateCPML_y(Hy);
    updateCPML_z(Hz);
}

template<class type1>
void cpml::updatePsiForEFields(const data3d<type1>& Hx, const data3d<type1>& Hy, const data3d<type1>& Hz) {
    if (is_cpml_xn) {
        updatePsi_eyz_n(Hz);
        updatePsi_ezy_n(Hy);
    }
    if (is_cpml_xp) {
        updatePsi_eyz_p(Hz);
        updatePsi_ezy_p(Hy);
    }
    if (is_cpml_yn) {
        updatePsi_exz_n(Hz);
        updatePsi_ezx_n(Hx);
    }
    if (is_cpml_yp) {
        updatePsi_exz_p(Hz);
        updatePsi_ezx_p(Hx);
    }
    if (is_cpml_zn) {
        updatePsi_eyx_n(Hx);
        updatePsi_exy_n(Hy);
    }
    if (is_cpml_zp) {
        updatePsi_eyx_p(Hx);
        updatePsi_exy_p(Hy);
    }
}

template<class type1>
void cpml::updatePsiForMFields(const data3d<type1>& Ex, const data3d<type1>& Ey, const data3d<type1>& Ez) {
    if (is_cpml_xn) {
        updatePsi_hyz_n(Ez);
        updatePsi_hzy_n(Ey);
    }
    if (is_cpml_xp) {
        updatePsi_hyz_p(Ez);
        updatePsi_hzy_p(Ey);
    }
    if (is_cpml_yn) {
        updatePsi_hxz_n(Ez);
        updatePsi_hzx_n(Ex);
    }
    if (is_cpml_yp) {
        updatePsi_hxz_p(Ez);
        updatePsi_hzx_p(Ex);
    }
    if (is_cpml_zn) {
        updatePsi_hyx_n(Ex);
        updatePsi_hxy_n(Ey);
    }
    if (is_cpml_zp) {
        updatePsi_hyx_p(Ex);
        updatePsi_hxy_p(Ey);
    }
}

template<class type1>
void cpml::updateEFieldCPML_x(data3d<type1>& Ey, data3d<type1>& Ez) {
    // x negetive region
    if (is_cpml_xn) {
        for (unsigned i = 0; i < n_cpml_xn; i++) {
            for (unsigned j = 0; j < Ey.ny; j++) {
                for (unsigned k = 0; k < Ey.nz; k++) {
                    Ey.p[i][j][k] += CPsi_hyx_xn.p[i][j][k] * Psi_hyx_xn.p[i][j][k];
                }
            }

            for (unsigned j = 0; j < Ez.ny; j++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ez.p[i][j][k] += CPsi_hzx_xn.p[i][j][k] * Psi_hzx_xn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_xp) {
        for (unsigned i = 0; i < n_cpml_xp; i++) {
            unsigned ihx = Ey.nx - n_cpml_xp - 1;
            for (unsigned j = 0; j < Ey.ny; j++) {
                for (unsigned k = 0; k < Ey.nz; k++) {
                    Ey.p[ihx][j][k] += CPsi_hyx_xp.p[i][j][k] * Psi_hyx_xp.p[i][j][k];
                }
            }
            ihx = Ez.nx - n_cpml_xp - 1;
            for (unsigned j = 0; j < Ez.ny; j++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ez.p[ihx][j][k] += CPsi_hzx_xp.p[i][j][k] * Psi_hzx_xp.p[i][j][k];
                }
            }
        }
    }
}

template<class type1>
void cpml::updateEFieldCPML_y(data3d<type1>& Ex, data3d<type1>& Ez) {
    // x negetive region
    if (is_cpml_yn) {
        for (unsigned i = 0; i < n_cpml_xn; i++) {
            for (unsigned j = 0; j < Ex.ny; j++) {
                for (unsigned k = 0; k < Ex.nz; k++) {
                    Ex.p[i][j][k] += CPsi_hxy_yn.p[i][j][k] * Psi_hxy_yn.p[i][j][k];
                }
            }
            for (unsigned j = 0; j < Ez.ny; j++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ez.p[i][j][k] += CPsi_hzy_yn.p[i][j][k] * Psi_hzy_yn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_yp) {
        for (unsigned i = 0; i < n_cpml_yp; i++) {
            unsigned ihx = Ex.nx - n_cpml_yp - 1;
            for (unsigned j = 0; j < Ex.ny; j++) {
                for (unsigned k = 0; k < Ex.nz; k++) {
                    Ex.p[ihx][j][k] += CPsi_hxy_yp.p[i][j][k] * Psi_hxy_yp.p[i][j][k];
                }
            }
            ihx = Ez.nx - n_cpml_yp - 1;
            for (unsigned j = 0; j < Ez.ny; j++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ez.p[ihx][j][k] += CPsi_hzy_yp.p[i][j][k] * Psi_hzy_yp.p[i][j][k];
                }
            }
        }
    }
}

template<class type1>
void cpml::updateEFieldCPML_z(data3d<type1>& Ex, data3d<type1>& Ey) {
    // x negetive region
    if (is_cpml_zn) {
        for (unsigned i = 0; i < n_cpml_zn; i++) {
            for (unsigned j = 0; j < Ey.ny; j++) {
                for (unsigned k = 0; k < Ey.nz; k++) {
                    Ey.p[i][j][k] += CPsi_hyz_zn.p[i][j][k] * Psi_hyz_zn.p[i][j][k];
                }
            }

            for (unsigned j = 0; j < Ez.ny; j++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ex.p[i][j][k] += CPsi_hxz_zn.p[i][j][k] * Psi_hxz_zn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_zp) {
        for (unsigned i = 0; i < n_cpml_zp; i++) {
            unsigned ihx = Ey.nx - n_cpml_zp - 1;
            for (unsigned j = 0; j < Ey.ny; j++) {
                for (unsigned k = 0; k < Ey.nz; k++) {
                    Ey.p[ihx][j][k] += CPsi_hyz_zn.p[i][j][k] * Psi_hyz_zn.p[i][j][k];
                }
            }
            ihx = Ex.nx - n_cpml_zp - 1;
            for (unsigned j = 0; j < Ez.ny; j++) {
                for (unsigned k = 0; k < Ez.nz; k++) {
                    Ex.p[ihx][j][k] += CPsi_hxz_zn.p[i][j][k] * Psi_hxz_zn.p[i][j][k];
                }
            }
        }
    }
}

template<class type1>
void cpml::updateMFieldCPML_x(data3d<type1>& Hy, data3d<type1>& Hz) {
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
        for (unsigned i = 0; i < n_cpml_xp; i++) {
            unsigned ihx = Hy.nx - n_cpml_xp - 1;
            for (unsigned j = 0; j < Hy.ny; j++) {
                for (unsigned k = 0; k < Hy.nz; k++) {
                    Hy.p[ihx][j][k] += CPsi_hyx_xp.p[i][j][k] * Psi_hyx_xp.p[i][j][k];
                }
            }
            ihx = Hz.nx - n_cpml_xp - 1;
            for (unsigned j = 0; j < Hz.ny; j++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hz.p[ihx][j][k] += CPsi_hzx_xp.p[i][j][k] * Psi_hzx_xp.p[i][j][k];
                }
            }
        }
    }
}

template<class type1>
void cpml::updateMFieldCPML_y(data3d<type1>& Hx, data3d<type1>& Hz) {
    // x negetive region
    if (is_cpml_yn) {
        for (unsigned i = 0; i < n_cpml_xn; i++) {
            for (unsigned j = 0; j < Hx.ny; j++) {
                for (unsigned k = 0; k < Hx.nz; k++) {
                    Hx.p[i][j][k] += CPsi_hxy_yn.p[i][j][k] * Psi_hxy_yn.p[i][j][k];
                }
            }
            for (unsigned j = 0; j < Hz.ny; j++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hz.p[i][j][k] += CPsi_hzy_yn.p[i][j][k] * Psi_hzy_yn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_yp) {
        for (unsigned i = 0; i < n_cpml_yp; i++) {
            unsigned ihx = Hx.nx - n_cpml_yp - 1;
            for (unsigned j = 0; j < Hx.ny; j++) {
                for (unsigned k = 0; k < Hx.nz; k++) {
                    Hx.p[ihx][j][k] += CPsi_hxy_yp.p[i][j][k] * Psi_hxy_yp.p[i][j][k];
                }
            }
            ihx = Hz.nx - n_cpml_yp - 1;
            for (unsigned j = 0; j < Hz.ny; j++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hz.p[ihx][j][k] += CPsi_hzy_yp.p[i][j][k] * Psi_hzy_yp.p[i][j][k];
                }
            }
        }
    }
}

template<class type1>
void cpml::updateMFieldCPML_z(data3d<type1>& Hx, data3d<type1>& Hy) {
    //TODO add update Hz in pml region
    // x negetive region
    if (is_cpml_zn) {
        for (unsigned i = 0; i < n_cpml_zn; i++) {
            for (unsigned j = 0; j < Hy.ny; j++) {
                for (unsigned k = 0; k < Hy.nz; k++) {
                    Hy.p[i][j][k] += CPsi_hyz_zn.p[i][j][k] * Psi_hyz_zn.p[i][j][k];
                }
            }

            for (unsigned j = 0; j < Hz.ny; j++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hx.p[i][j][k] += CPsi_hxz_zn.p[i][j][k] * Psi_hxz_zn.p[i][j][k];
                }
            }
        }
    }
    // x positive region
    if (is_cpml_zp) {
        for (unsigned i = 0; i < n_cpml_zp; i++) {
            unsigned ihx = Hy.nx - n_cpml_zp - 1;
            for (unsigned j = 0; j < Hy.ny; j++) {
                for (unsigned k = 0; k < Hy.nz; k++) {
                    Hy.p[ihx][j][k] += CPsi_hyz_zn.p[i][j][k] * Psi_hyz_zn.p[i][j][k];
                }
            }
            ihx = Hx.nx - n_cpml_zp - 1;
            for (unsigned j = 0; j < Hz.ny; j++) {
                for (unsigned k = 0; k < Hz.nz; k++) {
                    Hx.p[ihx][j][k] += CPsi_hxz_zn.p[i][j][k] * Psi_hxz_zn.p[i][j][k];
                }
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_eyz_p(const data3d<type1>& Hx) {
    for (unsigned k = 0, ikz = Hx.nz - n_cpml_zp; k < n_cpml_zp; k++, ihz++) {
        for (unsigned i = 0; i < Psi_eyz_zp.nx; i++) {
            for (unsigned j = 0; j < Psi_eyz_zp.ny; j++) {
                Psi_eyz_zp.p[i][j][k] = Psi_eyz_zp.p[i][j][k] * cpml_b_ez_zp.p[k] + cpml_a_ez_zp.p[k]*(Hx.p[i][j][ikz] - Hx.p[i][j][ikz - 1]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_ezy_p(const data3d<type1>& Hx) {
    for (unsigned j = 0, jhy = Hx.ny - n_cpml_yp; j < n_cpml_yp; j++, jhy++) {
        for (unsigned i = 0; i < Psi_ezy_yp.nx; i++) {
            for (unsigned k = 0; k < Psi_ezy_yp.nz; k++) {
                Psi_ezy_yp.p[i][j][k] = Psi_ezy_yp.p[i][j][k] * cpml_b_ey_yp.p[j] + cpml_a_ey_yp.p[j]*(Hx.p[i][jhy][k] - Hx.p[i][jhy - 1][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_exz_p(const data3d<type1>& Hy) {
    for (unsigned k = 0, iez = Hy.nz - n_cpml_zp - 2; k < n_cpml_zp; k++, iez++) {
        for (unsigned i = 0; i < Psi_exz_zp.nx; i++) {
            for (unsigned j = 0; j < Psi_exz_zp.ny; j++) {
                Psi_exz_zp.[i][j][k] = cpml_b_ez_zp.p[k] * Psi_exz_zp.p[i][j][k] + cpml_a_ez_zp.p[k]*(Hy.p[i][j][iez + 1] - Hy.p[i][j][iez]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_ezx_p(const data3d<type1>& Hy) {
    for (unsigned i = 0, iex = Hz.nx - n_cpml_xp - 1; i < n_cpml_xp; i++, iex++) {
        for (unsigned j = 0; j < Psi_ezx_xp.ny; j++) {
            for (unsigned k = 0; k < Psi_ezx_xp.nz; k++) {
                Psi_ezx_xp.p[i][j][k] = cpml_b_ex_xp.p[i] * Psi_ezx_xp.p[i][j][k] + cpml_a_ex_xp.p[i]*(Hz.p[iex][j][k] - Hz.p[iex - 1][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_eyx_p(const data3d<type1>& Hz) {
    for (unsigned i = 0; i < n_cpml_xp; i++) {
        for (unsigned j = 0; j < Psi_eyx_xp.ny; j++) {
            for (unsigned k = 0; k < Psi_eyx_xp.nz; k++) {
                Psi_eyx_xp.p[i][j][k] = cpml_b_ex_xp.p[i] * Psi_eyx_xp.p[i][j][k] + cpml_a_ex_xp.p[i]*(Hz.p[i + 1][j][k] - Hz.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_exy_p(const data3d<type1>& Hz) {
    for (unsigned j = 0, imy = Hz.ny - n_cpml_yp; j < n_cpml_yp; j++, imy++) {
        for (unsigned i = 0; i < Psi_exy_yp.nx; i++) {
            for (unsigned k = 0; k < Psi_exy_yp.nz; k++) {
                Psi_exy_yp.p[i][j][k] = Psi_exy_yp.p[i][j][k] * cpml_b_ey_yp.p[j] + cpml_a_ey_yp.p[j]*(Ez.p[i][imy][k] - Ez.p[i][imy - 1][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_hyz_p(const data3d<type1>& Ex) {
    for (unsigned k = 0, ikz = Ex.nz - n_cpml_zp; k < n_cpml_zp; k++, ihz++) {
        for (unsigned i = 0; i < Psi_hyz_zp.nx; i++) {
            for (unsigned j = 0; j < Psi_hyz_zp.ny; j++) {
                Psi_hyz_zp.p[i][j][k] = Psi_hyz_zp.p[i][j][k] * cpml_b_mz_zp.p[k] + cpml_a_mz_zp.p[k]*(Ex.p[i][j][ikz] - Ex.p[i][j][ikz - 1]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_hzy_p(const data3d<type1>& Ex) {
    for (unsigned j = 0, jhy = Ex.ny - n_cpml_yp - 1; j < n_cpml_yp; j++, jhy++) {
        for (unsigned i = 0; i < Psi_hzy_yp.nx; i++) {
            for (unsigned k = 0; k < Psi_hzy_yp.nz; k++) {
                Psi_hzy_yp.p[i][j][k] = Psi_hzy_yp.p[i][j][k] * cpml_b_my_yp.p[j] + cpml_a_my_yp.p[j]*(Ex.p[i][jhy][k] - Ex.p[i][jhy - 1][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_hxz_p(const data3d<type1>& Ey) {
    for (unsigned k = 0, ihz = Ey.nz - n_cpml_zp - 1; k < n_cpml_zp; k++, ihz++) {
        for (unsigned j = 0; j < Psi_hxz_zp.ny; j++) {
            for (unsigned i = 0; i < Psi_hxz_zp.nx; i++) {
                Psi_hxz_zp.p[i][j][k] = Psi_hxz_zp.p[i][j][k] * cpml_b_mz_zp.p[k] + cpml_a_mz_zp.p[k]*(Ey.p[i][j][ihz] - Ey.p[i][j][ihz - 1])
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_hzx_p(const data3d<type1>& Ey) {
    for (unsigned i = 0, ihx = Ey.nx - n_cpml_xp - 1; i < n_cpml_xp; i++, ihx++) {
        for (unsigned j = 0; j < Psi_hzx_xp.ny; j++) {
            for (unsigned k = 0; k < Psi_hzx_xp.nz; k++) {
                Psi_hzx_xp.p[i][j][k] = Psi_hzx_xp.p[i][j][k] * cpml_b_mx_xp.p[i] + cpml_a_mx_xp.p[i]*(Ey.p[ihx + 1][j][k] - Ey.p[ihx][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_hyx_p(const data3d<type1>& Ez) {
    for (unsigned i = 0, imx = Ez.nx - n_cpml_zp - 1; i < n_cpml_zp; i++, imx++) {
        for (unsigned j = 0; j < Psi_hyx_xp.ny; j++) {
            for (unsigned k = 0; k < Psi_hyx_xp.nz; k++) {
                Psi_hyx_xp.p[i][j][k] = Psi_hyx_xp.p[i][j][k] * cpml_b_mx_xp.p[i] + cpml_a_mx_xp.p[i]*(Ez.p[i + 1][j][k] - Ez.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_hxy_p(const data3d<type1>& Ez) {
    for (unsigned j = 0, imy = Ez.ny - n_cpml_yp - 1; j < n_cpml_yp; j++, imy++) {
        for (unsigned i = 0; i < Psi_hxy_yp.nx; i++) {
            for (unsigned k = 0; k < Psi_hxy_yp.nz; k++) {
                Psi_hxy_yp.p[i][j][k] = Psi_hxy_yp.p[i][j][k] * cpml_b_my_yp.p[j] + cpml_a_my_yp.p[j]*(Ez.p[i][imy + 1][k] - Ez.p[i][imy][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_eyz_n(const data3d<type1>& Hx) {
    for (unsigned k = 0; k < n_cpml_zn; k++) {
        for (unsigned i = 0; i < Psi_eyz_zn.nx; i++) {
            for (unsigned j = 0; Psi_eyz_zn.ny; j++) {
                Psi_eyz_zn.p[i][j][k] = cpml_b_ez_zn.p[k] * Psi_eyz_zn.p[i][j][k] + cpml_a_ez_zn.p[k]*(Hx.p[i][j][k + 1] - Hx.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_ezy_n(const data3d<type1>& Hx) {
    for (unsigned j = 0; j < n_cpml_yn; j++) {
        for (unsigned i = 0; i < Psi_ezy_yn.nx; i++) {
            for (unsigned k = 0; k < Psi_ezy_yn.nz; k++) {
                Psi_ezy_yn.[i][j][k] = cpml_b_ey_yn.p[j] * Psi_ezy_yn.p[i][j][k] + cpml_a_ey_yn.p[j]*(Hx.p[i][j + 1][k] - Hx.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_exz_n(const data3d<type1>& Hy) {
    for (unsigned k = 0; k < n_cpml_zn; k++) {
        for (unsigned i = 0; i < Psi_exz_zn.nx; i++) {
            for (unsigned j = 0; j < Psi_exz_zn.ny; j++) {
                Psi_exz_zn.[i][j][k] = cpml_b_ez_zn.p[k] * Psi_exz_zn.p[i][j][k] + cpml_a_ez_zn.p[k]*(Hy.p[i][j][k + 1] - Hy.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_ezx_n(const data3d<type1>& Hy) {
    for (unsigned i = 0; i < n_cpml_xn; i++) {
        for (unsigned j = 0; j < Psi_ezx_xn.ny; j++) {
            for (unsigned k = 0; k < Psi_ezx_xn.nz; k++) {
                Psi_ezx_xn.p[i][j][k] = cpml_b_ex_xn.p[i] * Psi_ezx_xn.p[i][j][k] + cpml_a_ex_xn.p[i]*(Hy.p[i + 1][j][k] - Hy.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_eyx_n(const data3d<type1>& Hz) {
    for (unsigned i = 0; i < n_cpml_xn; i++) {
        for (unsigned j = 0; j < Psi_eyx_xn.ny; j++) {
            for (unsigned k = 0; k < Psi_eyx_xn.nz; k++) {
                Psi_eyx_xn.p[i][j][k] = cpml_b_ex_xn.p[i] * Psi_eyx_xn.p[i][j][k] + cpml_a_ex_xn.p[i]*(Hz.p[i + 1][j][k] - Hz.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_exy_n(const data3d<type1>& Hz) {
    for (unsigned j = 0; j < n_cpml_yn; j++) {
        for (unsigned k = 0; k < Psi_exy_yn.nz; k++) {
            for (unsigned i = 0; i < Psi_exy_yn.ny; i++) {
                Psi_exy_yn.p[i][j][k] = cpml_b_ey_yn.p[j] * Psi_exy_yn.p[i][j][k] + cpml_a_ey_yn.p[j]*(Hz.p[i][j + 1][k] - Hz.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_hyz_n(const data3d<type1>& Ex) {
    for (unsigned k = 0; k < n_cpml_zn; k++) {
        for (unsigned j = 0; j < Psi_hyz_zn.ny; j++) {
            for (unsigned i = 0; i < Psi_hyz_zn.nx; i++) {
                Psi_hyz_zn.p[i][j][k] = cpml_b_mz_zn.p[k] * Psi_hyz_zn.p[i][j][k] + cpml_a_mz_zn.p[k]*(Ex.p[i][j][k + 1] - Ex.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_hzy_n(const data3d<type1>& Ex) {
    for (unsigned j = 0; j < n_cpml_yn; j++) {
        for (unsigned i = 0; i < Psi_hzy_yn.nx; i++) {
            for (unsigned k = 0; k < Psi_hzy_yn.nz; k++) {
                Psi_hzy_yn.p[i][j][k] = cpml_b_my_yn.p[j] * Psi_hzy_yn.p[i][j][k] + cpml_a_my_yn.p[j]*(Ex.p[i][j + 1][k] - Ex.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_hxz_n(const data3d<type1>& Ey) {
    for (unsigned k = 0; k < n_cpml_zn; k++) {
        for (unsigned i = 0; i < Psi_hxz_zn.nx; i++) {
            for (unsigned j = 0; j < Psi_hxz_zn.ny; j++) {
                Psi_hxz_zn.p[i][j][k] = cpml_b_mz_zn.p[k] * Psi_hxz_zn.p[i][j][k] + cpml_a_mz_zn.p[k]*(Ey.p[i][j][k + 1] - Ey.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_hzx_n(const data3d<type1>& Ex) {
    for (unsigned i = 0; i < n_cpml_xn; i++) {
        for (unsigned j = 0; j < Psi_hzx_xn.ny; j++) {
            for (unsigned k = 0; k < Psi_hzx_xn.nz; k++) {
                Psi_hzx_xn.p[i][j][k] = Psi_hzx_xn.p[i][j][k] * cpml_b_mx_zn.p[i] + cpml_a_mx_zn.p[i]*(Ex.p[i + 1][j][k] - Ex.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_hyx_n(const data3d<type1>& Ex) {
    for (unsigned i = 0; i < n_cpml_xn; i++) {
        for (unsigned j = 0; j < Psi_hyz_zn.ny; j++) {
            for (unsigned k = 0; k < Psi_hyz_zn.nz; k++) {
                Psi_hyz_zn.p[i][j][k] = Psi_hyz_zn.p[i][j][k] * cpml_b_mz_zn.p[i] + cpml_a_mz_zn.p[i]*(Ex.p[i + 1][j][k] - Ex.p[i][j][k]);
            }
        }
    }
}

template<class type1>
void cpml::updatePsi_hxy_n(const data3d<type1>& Ez) {
    for (unsigned j = 0; j < n_cpml_yn; j++) {
        for (unsigned i = 0; i < Psi_hxy_yn.nx; i++) {
            for (unsigned k = 0; k < Psi_hxy_yn.nz; k++) {
                Psi_hxy_yn.p[i][j][k] = Psi_hxy_yn.p[i][j][k] * cpml_b_my_yn.p[j] + cpml_a_my_yn.p[j]*(Ez.p[i][j + 1][k] - Ez.p[i][j][k])
            }
        }
    }
}

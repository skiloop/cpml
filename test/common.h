
#ifndef COMMON_H
#define COMMON_H

#include <math.h>

typedef double MyDataF;
typedef MyDataF* pMyDataF;
typedef pMyDataF* ppMyDataF;

const MyDataF C = 2.99792458E8; // speed of light in free space
const MyDataF me = 9.110e-31; // electricity mass
const MyDataF e = 1.602e-19; // electricity charge
const MyDataF mu_0 = 4.0 * M_PI * 1.0E-7;
const MyDataF eps_0 = 1.0 / (C * C * mu_0);
const MyDataF M_PI_TWO = M_PI * 2;

#endif // COMMON_H


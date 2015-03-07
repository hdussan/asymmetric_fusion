/*
 *
 */
#ifndef potentials_h
#define potentials_h
#include <cmath>
#include "constants.h"

 using namespace std;

 double coulomb(double r, double Ro1, double Ro2, int Z1,int Z2);

 double V_l(double r, double reducedMass, int l);

#endif

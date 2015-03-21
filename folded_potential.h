/******************************************************************
 *                         folded_potential.h                     *
 *                         ASYMMETRIC FUSION                      *
 *  DEFINITION OF POTENTIALS FOR THE BARRIER POTENTIAL FOMALISM   *
 *  Vf CALCULATED FROM W-S DENSITIES                              *
 *  n IS THE SIZE OF THE VECTOR X.                                *
 *****************************************************************/
#ifndef folded_potential_h
#define folded_potential_h
#include "constants.h"
#include "Gauss_q.h"

using namespace std;

double R_o(int baryonNumber);

double diffuseness(int protonNumber);

double central_density(double Radius, double a, int baryonNumber, 
                       int protonNumber);

double ws_density(double Radius, int baryonNumber, int protonNumber, 
                  double a, double density0, double radialDistance);

double folded_potential(double r, double wsRadius1, double wsRadius2,
                        int protonNumber1,int protonNumber2,
                        int baryonNumber1,int baryonNumber2);


#endif



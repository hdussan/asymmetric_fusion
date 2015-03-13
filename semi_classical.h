/*
 *                       semi_classical.h
 *  FUSION ^Z1+^Z2 REACTIONS
 *  folded_potential files available to calculate
 *  the Effective potential and return points in the
 *  semi-classical approach
 */

#ifndef semi_classical_h
#define semi_classical_h
#include "constants.h"
#include "potentials.h"
#include "folded_potential.h"
#include <cmath>
#include <vector>
using namespace std;
/*
  Reduced mass
 */
double mReduce(int B1, int B2);

double menosEffectiveV(double r, double Energy, 
                        double wsRadius1, double wsRadius2,
			int protonNumber1, int protonNumber2,
			int baryonNumber1, int baryonNumber2,
		        double reducedMass, int orbitalAngularMomentum);

double PotentialMenusEnergy(double r, double Energy,
                             double wsRadius1, double wsRadius2,
			     int protonNumber1, int protonNumber2,
			     int baryonNumber1, int baryonNumber2,
			     double reducedMass, int orbitalAngularMomentum);


double middleGuess(double r_infinity,  double Energy,
		   double wsRadius1, double wsRadius2,
		   int protonNumber1, int protonNumber2,
		   int baryonNumber1, int baryonNumber2,
		   double reducedMass, int orbitalAngularMomentum);
/********************************************************
 *   FUNCTIONS	SELF-CONSISTENT VELOCITY CALCULATION	*
 * To be declared  Don't have to be templates    !!!!!  *
 ********************************************************/


//// Self-consistent calculation of the velocity
//// Sao Paulo potential 
template<class T>
inline T f_v2(T A1,T A2,T v2)
{ 
  return A1 - A2 * exp(-4 * v2);
}

template<class T>
T speedy2(T En,T A1,T A2)
{
 T vold2,vnew2,vsonew2,v2;
 v2=0.02;
 inicio:
        vold2=v2;
        vnew2=f_v2(A1,A2,vold2);
        vsonew2=0.9*vnew2+(1-0.9)*vold2;
	if(fabs(vsonew2-vold2)>1.0e-15)
	  { 
           v2=vsonew2;
	   goto inicio;
          }
return v2;
}

#endif



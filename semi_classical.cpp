/*
 *                       semi_classical.cpp
 *   FUSION ^Z1+^Z2 REACTIONS
 *
 ********************************************************/

#include "semi_classical.h"
using namespace std;

/*     
  Reduced mass
 */
double mReduce(int B1,  int B2)
{
  return B1 * B2 * Uma /(B1 + B2);
}

/*
 * Sao Paulo potential non-locality
 * Incorporated by puting an energy dependent 
 * term:
 *  Velocity is calculated self-consistently
 */
double velocityDependency( double coeff1, double coeff2, double velocity2)
{
  return coeff1 - coeff2 * exp(-4 * velocity2);
}

double selfconsistentV2(double Energy, double coeff1, double coeff2)
{
  double vold2, vnew2, vsonew2,v2;
  v2 = 0.02;
  double const weight = 0.9;
  inicio:
    vold2 = v2;
    vnew2 = velocityDependency( coeff1, coeff2, vold2);
    vsonew2 = weight * vnew2 + (1 - weight) * vold2;
    if(fabs( vsonew2 - vold2) > 1.0e-15)
    { 
      v2 = vsonew2;
      goto inicio;
    }
  return v2;
}
/*
 *   Effective potential: (it returns  - (Effective potential) )
 *                - Coulomb - Sao Paulo - centrifugal potential
 *   Sao Paulo = Folded Potential * exp( -4 v^2)
 */
 double menosEffectiveV(double r, double Energy, 
                        double wsRadius1, double wsRadius2,
			int protonNumber1, int protonNumber2,
			int baryonNumber1, int baryonNumber2,
                        double reducedMass, int orbitalAngularMomentum)
 {
   double Ro1 = wsRadius1;
   double Ro2 = wsRadius2;
   double mu = reducedMass;
   int l =  orbitalAngularMomentum;
   double Vc = coulomb(r, Ro1, Ro2, protonNumber1, protonNumber2);
   double Vl = V_l(r, mu, l);
   double Vfolded = folded_potential(r, Ro1, Ro2, protonNumber1, protonNumber2,
                                     baryonNumber1, baryonNumber2);
   double coeff1 = 2 * (Energy - Vc)/ mu; 
   double coeff2 = 2 * Vfolded / mu;
   double v2 = selfconsistentV2(Energy, coeff1, coeff2);
   
   return -Vc - Vfolded * exp(-4 * v2) - Vl;
 }

/*
 *   Effective Potential - Energy
 *  To find the return points   => 
 *  TODO: Need a new specialised bisection function (because of arguments)
 */

 double PotentialMenusEnergy(double r, double Energy,
                             double wsRadius1, double wsRadius2,
			     int protonNumber1, int protonNumber2,
			     int baryonNumber1, int baryonNumber2,
                             double reducedMass, int orbitalAngularMomentum)
{
  double Veff = -menosEffectiveV(r, Energy, wsRadius1, wsRadius2,
			         protonNumber1, protonNumber2, 
                                  baryonNumber1, baryonNumber2,
				   reducedMass, orbitalAngularMomentum); 
  return Veff - Energy;
}

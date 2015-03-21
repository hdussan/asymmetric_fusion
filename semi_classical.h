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
#include "minimal.h"
#include "Derivate.h"
#include <cmath>
#include <vector>
using namespace std;

double const hbarc2 = hbarc * hbarc;
/*
  Reduced mass
 */
double mReduce(int B1, int B2);

double menosEffectiveV(double r, double Energy, 
                        double wsRadius1, double wsRadius2,
			int protonNumber1, int protonNumber2,
			int baryonNumber1, int baryonNumber2,
		        double reducedMass, int orbitalAngularMomentum);

double potentialMenusEnergy(double r, double Energy,
                            double wsRadius1, double wsRadius2,
			    int protonNumber1, int protonNumber2,
			    int baryonNumber1, int baryonNumber2,
			    double reducedMass, int orbitalAngularMomentum);


double middleGuess(double r_infinity,  double Energy,
		   double wsRadius1, double wsRadius2,
		   int protonNumber1, int protonNumber2,
		   int baryonNumber1, int baryonNumber2,
		   double reducedMass, int orbitalAngularMomentum);

double closeToOriginGuess(double r_o, double Energy,
		          double wsRadius1, double wsRadius2,
		          int protonNumber1, int protonNumber2,
		          int baryonNumber1, int baryonNumber2,
			  double reducedMass, int orbitalAngularMomentum);


double TransmissionCoeff(double r2, double r1, double Energy,
		         double wsRadius1, double wsRadius2,
		         int protonNumber1, int protonNumber2,
		         int baryonNumber1, int baryonNumber2,
			 double reducedMass, int orbitalAngularMomentum);


void getSfactorAndCrossSection(double r, double rAtInfinity, double Energy,
		               double wsRadius1, double wsRadius2,
		               int protonNumber1, int protonNumber2,
		               int baryonNumber1, int baryonNumber2,
			       double reducedMass,
                               double &Sfactor, double &crossSection);


#endif



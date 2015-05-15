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
#include <fstream>
//#include <pthread.h>
#include <boost/thread.hpp>
#include <boost/bind.hpp>

using namespace std;
using namespace boost;
double const hbarc2 = hbarc * hbarc;

/*
  Reduced mass
 */
double mReduce(int B1, int B2);

double selfconsistentV2(double Energy, double coeff1, double coeff2);

/*********************************************************
 *********************************************************/
class fusionSemiClassical
{
  private:
    double wsRadius1;
    int protonNumber1;
    int baryonNumber1;
    double wsRadius2;
    int protonNumber2;
    int baryonNumber2;
    double reducedMass;
  public:
    fusionSemiClassical(double R1, int Z1, int A1, double R2, int Z2, int A2, double mu);

    double menosEffectiveV(double r, double Energy, int orbitalAngularMomentum);

    double potentialMenusEnergy(double r, double Energy, int orbitalAngularMomentum);

    double middleGuess(double r_infinity,  double Energy, int orbitalAngularMomentum);

    double findReturnPoint(double rMenus, double rPlus, double Energy,
                           int orbitalAngularMomentum);

    double differentiateV(double const &xo, double const &h,
                          double input4function1, int input4function2);

    double closeToOriginGuess(double r_o, double Energy, int orbitalAngularMomentum);

    double TransmissionCoeff(double r2, double r1, double Energy, int orbitalAngularMomentum);

    double TransmissionCoeffHW(double b_, double rAtInfinity, double Energy, 
                               int orbitalAngularMomentum); 
    /*
    void getSfactorAndCrossSection(double r, double rAtInfinity, double Energy,
                                   double &Sfactor, double &crossSection);
    */
    void getSfactorAndCrossSection(double r, double rAtInfinity, double Energy,
                                   vector<double> &Sfactor, vector<double> &crossSection);

    double brentMinimise(const double ax, const double bx, const double cx,
                         const double input4function1, const int input4function2, 
			 const double tol, double &xmin);

    void saveAstrophysicalSfactor(vector<double> &energies);

};

void solveFusion(fusionSemiClassical &XY, double minEnergy,double maxEnergy, int numberData,
                 vector<double> &energyValues,  vector<double> &Sfactor, vector<double> &crossSection);

#endif



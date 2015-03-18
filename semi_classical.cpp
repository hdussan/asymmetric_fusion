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

 double potentialMenusEnergy(double r, double Energy,
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

/*********************************************************
 *  midGuess & smallGuess				 *
 *CHECK FOR GOOD GUESSES FOR FINDING THE TURNING POINTS	 *
 *INPUTS r_inf  or r_o in fm				 *
 *********************************************************/

double middleGuess(double r_infinity,  double Energy,
		   double wsRadius1, double wsRadius2,
		   int protonNumber1, int protonNumber2,
		   int baryonNumber1, int baryonNumber2,
		   double reducedMass, int orbitalAngularMomentum)
{
  double rAtInfinity = r_infinity;
  double const step = 0.2;
  int l = orbitalAngularMomentum;
  int z1 = protonNumber1;
  int z2 = protonNumber2;
  int B1 = baryonNumber1;
  int B2 = baryonNumber2;
  double mu = reducedMass;
  double a = 0;

  double deltaEnergy = potentialMenusEnergy(rAtInfinity, Energy,
                                            wsRadius1, wsRadius2,
			                    z1, z2, B1, B2, mu, l);
  
  if(deltaEnergy > 0)
  {
    a = rAtInfinity;
  }
  else
  {
    here:
    rAtInfinity -= step;
    deltaEnergy = potentialMenusEnergy(rAtInfinity, Energy, 
                                       wsRadius1, wsRadius2,
			               z1, z2, B1, B2, mu, l);
    if(deltaEnergy > 0)
    {
      a = rAtInfinity;
    }
    else
    {
      goto here;
    }
  } 
 return a;
}
/*
 * To replace smallGuess
 */
double closeToOriginGuess(double r_o, double Energy,
		          double wsRadius1, double wsRadius2,
		          int protonNumber1, int protonNumber2,
		          int baryonNumber1, int baryonNumber2,
			  double reducedMass, int orbitalAngularMomentum)
{
  int const maxIterations = 50;
  double En = Energy;
  double R1 = wsRadius1;
  double R2 = wsRadius2; 
  int z1 = protonNumber1;
  int z2 = protonNumber2;
  int B1 = baryonNumber1;
  int B2 = baryonNumber2;
  double mu = reducedMass;
  int l = orbitalAngularMomentum;

  double c=0, cplus, cminos;

  Zero<double, potentialMenusEnergy>(r_o, 9.0, En, R1, R2, z1, z2, B1, B2, mu, l, c, maxIterations);
  cplus  = c + 2 * r_o;
  cminos = c - 2 * r_o;
  double effectiveE = potentialMenusEnergy(cplus, En, R1, R2, z1, z2, B1, B2, mu, l);
  return (effectiveE < 0)? cplus:fabs(cminos); 
}

/********************************************************
 *    PROBABILITY OF PENETRATING THE BARRIER POTENTIAL	*
 *                                                      *
 *GIVEN TURNING POINTS r1 AND r2 GETS Tl		*
 *Tl   TRANSMISSION COEF.				*
 *Sl   PENETRATION PROBABILITY				*

 *********************************************************/
// To replace Tl method template
double TransmissionCoeff(double r2, double r1, double Energy,
		          double wsRadius1, double wsRadius2,
		           int protonNumber1, int protonNumber2,
		            int baryonNumber1, int baryonNumber2,
			     double reducedMass, int orbitalAngularMomentum) 
{
  int l = orbitalAngularMomentum;
  int z1 = protonNumber1;
  int z2 = protonNumber2;
  int B1 = baryonNumber1;
  int B2 = baryonNumber2;
  double mu = reducedMass;

  double Sl=0.0;
  vector<double> X(200), W(200), fl(200);
  GausLeg(r2,r1,X,W);
  for(int i = 0; i < X.size(); i++)
  {
    double r_aux=X[i];
    double effectiveEnergy = potentialMenusEnergy(r_aux, Energy, 
                                                  wsRadius1, wsRadius2,
						  z1, z2, B1, B2, mu, l);
    fl[i] = sqrt(8 * mu * fabs(effectiveEnergy));
  }
  Sl = dot(W, fl) / hbarc;
  return 1 / (1 + exp(Sl)); 
}
/***
 * To replace Tl_HW
 * TransmissionCoeffHW:	TRANSMISSION COEF. WHEN E>Vmax (HILL-WHEELER)	*
***********/

double TransmissionCoeffHW(double b_, double rAtInfinity, double Energy,
		           double wsRadius1, double wsRadius2,
		           int protonNumber1, int protonNumber2,
		           int baryonNumber1, int baryonNumber2,
			   double reducedMass, int orbitalAngularMomentum)
{
  double En = Energy;
  int l = orbitalAngularMomentum;
  int z1 = protonNumber1;
  int z2 = protonNumber2;
  int B1 = baryonNumber1;
  int B2 = baryonNumber2;
  double mu = reducedMass;

  double  Rmax, Vmax, omegal;
  double  Sl;
  double smallGuess = 1.2; //fm
  double  middle = fabs(0.5 * rAtInfinity);

  Vmax = -BRENT(smallGuess, middle, rAtInfinity,
                En, wsRadius1, wsRadius2,
                z1, z2, B1, B2, mu, l,  
                potentialMenusEnergy, 1e-7, Rmax); 

  omegal = fabs(D2_f5(Rmax,0.001,En, wsRadius1, wsRadius2,
                       z1, z2, B1, B2, mu, l, potentialMenusEnergy));

  omegal = hbarc * sqrt(omegal / mu);

  Sl = 2 * pi * (Vmax - En)/omegal;

  return 1 / (1 + exp(Sl));
}




/*********************************************************
 *   S_fact	  ASTROPHYSICAL S-FACTOR	      	 *
 *       Outputs                                         *
 *   Sfactor : Astrophysical S-factor                    *
 *   crossSection: cross section                         * 
 *********************************************************/
 // To replace S_fact
void getSfactorAndCrossSection(double r, double rAtInfinity, double Energy,
		               double wsRadius1, double wsRadius2,
		               int protonNumber1, int protonNumber2,
		               int baryonNumber1, int baryonNumber2,
			       double reducedMass,
                               double &Sfactor, double &crossSection)
{
  double En = Energy;
  double R1 = wsRadius1;
  double R2 = wsRadius2;
  int z1 = protonNumber1;
  int z2 = protonNumber2;
  int B1 = baryonNumber1;
  int B2 = baryonNumber2;
  double mu = reducedMass;

  double a_, b_, r1, r2;
  int local;
  Sfactor = 0.0;
  vector<double> T_l(21), subl(21);
  cout<<"calculating fusion at E = "<<En<<"MeV\n";
  double gamow = z1 * z2 * e2 * sqrt( mu / ( 2 * En * hbarc2));
  local = (En > 7.2)? 1:0;
  for(int l = local; l < T_l.size(); l++)
  {
    a_ = closeToOriginGuess(0.04, En, R1, R2, z1, z2, B1, B2, mu, l);
    b_ = middleGuess(r, En, R1, R2, z1, z2, B1, B2, mu, l);

    Zero<double, potentialMenusEnergy>(a_, b_, En, R1, R2, 
                                        z1, z2, B1, B2, mu, l, r2,100);
    Zero<double, potentialMenusEnergy>(b_, r, En, R1, R2, 
                                         z1, z2, B1, B2, mu, l, r1, 100);

    double deltar=fabs(r1-r2);
    //L = l;
    
    if(deltar < 1e-5 )
    {
      T_l[l] = TransmissionCoeffHW(b_, rAtInfinity, En, 
                                     R1, R2, z1, z2, B1, B2, mu, l);
    } 
    else
    {
      T_l[l] = TransmissionCoeff(r2, r1, En, 
                                   R1, R2, z1, z2, B1, B2, mu, l);
    } 
   
    subl[l] = 2 * l + 1;

    cout.precision(8);
    cout<<l<<" \t"<<r2<<" \t"<<r1<<" \t"<<T_l[l]<<endl;
  }

  Sfactor = pi * hbarc2 * dot(subl, T_l) * exp(2 * pi * gamow) / (2 * mu * 100); 
  crossSection = pi * hbarc2 * dot(subl, T_l)* 10 / (2 * mu * En);
}

/*
 *                       semi_classical.cpp
 *   FUSION ^Z1+^Z2 REACTIONS
 *
 ********************************************************/

#include "semi_classical.h"
using namespace std;

 inline void shift3( double &a, double &b, double &c,const double d)  
  { 
    a = b; 
    b = c; 
    c = d;
  }

/*     Reduced mass     */
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


/******************************************************

 ******************************************************/
fusionSemiClassical::fusionSemiClassical(double R1, int Z1, int A1, 
                                         double R2, int Z2, int A2, 
                                         double mu)
{
   wsRadius1  = R1;
   protonNumber1 = Z1;
   baryonNumber1 = A1;
   wsRadius2  = R2;
   protonNumber2 = Z2;
   baryonNumber2 = A2;
   reducedMass = mu;

}

/*
 *   Effective potential: (it returns  - (Effective potential) )
 *                - Coulomb - Sao Paulo - centrifugal potential
 *   Sao Paulo = Folded Potential * exp( -4 v^2)
 */
 double fusionSemiClassical::menosEffectiveV(double r, double Energy, 
                                             int orbitalAngularMomentum)
 {
   double Ro1 = wsRadius1;
   double Ro2 = wsRadius2;
   double mu = reducedMass;

   int l =  orbitalAngularMomentum;
   double Vc = coulomb(r, Ro1, Ro2, protonNumber1, protonNumber2);
   double Vl = V_l(r, mu, l);
   double Vfolded = folded_potential(r, Ro1, Ro2, protonNumber1, 
                                     protonNumber2,
                                     baryonNumber1, baryonNumber2);
   double coeff1 = 2 * (Energy - Vc)/ mu; 
   double coeff2 = 2 * Vfolded / mu;
  
   double v2 = selfconsistentV2(Energy, coeff1, coeff2);

   return -Vc - Vfolded * exp(-4 * v2) - Vl;
 }

double fusionSemiClassical::potentialMenusEnergy(double r, double Energy, 
                                                 int orbitalAngularMomentum)
{
  double Veff = -menosEffectiveV(r, Energy, orbitalAngularMomentum); 
  return Veff - Energy;
}

/*********************************************************
 *CHECK FOR GOOD GUESSES FOR FINDING THE TURNING POINTS  *
 *INPUTS r_inf  or r_o in fm				 *
 *********************************************************/
double fusionSemiClassical::middleGuess(double r_infinity, double Energy, 
                                        int orbitalAngularMomentum)
{
  double rAtInfinity = r_infinity;
  double const step = 0.5;
  int l = orbitalAngularMomentum;

  double a = 0;
  double deltaEnergy = potentialMenusEnergy(rAtInfinity, Energy, l);
  
  if(deltaEnergy > 0)
  {
    a = rAtInfinity;
  }
  else
   {
    
    here:
    rAtInfinity -= step;
    deltaEnergy = potentialMenusEnergy(rAtInfinity, Energy, l);
   
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

double fusionSemiClassical::findReturnPoint(double rMenus, double rPlus, 
                                            double Energy,
                                            int orbitalAngularMomentum)
{
  int const maxIterations = 128;
  double E = Energy;
  int l = orbitalAngularMomentum;

  double X;
  double point, newpoint;
  double a = rMenus;
  double b = rPlus;
  double c = (a + b) / 2.;
  X = c;
  for(int j = 0; j <= maxIterations; j++)
  {
    point = c;
    (potentialMenusEnergy(c, E, l) * potentialMenusEnergy(b, E, l) > 0)? b = point: a = point;
    newpoint = (a + b) / 2.;
    if(fabs(potentialMenusEnergy(newpoint, E, l)) < 1e-10)
    {
      X=newpoint;
      break;
    }
    else  c = newpoint;
   }

  return X;
}
// 5 points derivative 
double fusionSemiClassical::differentiateV(double const &xo, double const &h,
                                           double input4function1, 
                                           int input4function2)
{
  double f1 = potentialMenusEnergy(xo - 2 * h, input4function1, input4function2);
  double f2 = potentialMenusEnergy(xo - h, input4function1, input4function2);
  double f0 = potentialMenusEnergy(xo, input4function1, input4function2);
  double f3 = potentialMenusEnergy(xo + h, input4function1, input4function2);
  double f4 = potentialMenusEnergy(xo + 2 * h, input4function1, input4function2);

  double delta2_f = -f1 + 16 * f2 - 30 * f0 + 16 * f3 - f4;

  delta2_f = delta2_f / 12.;
  return delta2_f / (h * h);
}

double fusionSemiClassical::closeToOriginGuess(double r_o, double Energy,
                                               int orbitalAngularMomentum)
{
  int const maxIterations = 64;// should be 100 for the return points
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
  double r0plus = 9.0;//fm
  c = findReturnPoint(r_o, r0plus, En, l);

  cplus  = c + 2 * r_o;
  cminos = c - 2 * r_o;
  double effectiveE = potentialMenusEnergy(cplus, En, l);
  return (effectiveE < 0)? cplus:fabs(cminos); 
}
/********************************************************
 *    PROBABILITY OF PENETRATING THE BARRIER POTENTIAL	*
 *                                                      *
 *GIVEN TURNING POINTS r1 AND r2 GETS Tl		*
 *Tl   TRANSMISSION COEF.				*
 *Sl   PENETRATION PROBABILITY				*

 *********************************************************/

double fusionSemiClassical::TransmissionCoeff(double r2, double r1, 
                                              double Energy, 
                                              int orbitalAngularMomentum)
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
    double effectiveEnergy = potentialMenusEnergy(r_aux, Energy, l);
    fl[i] = sqrt(8 * mu * fabs(effectiveEnergy));
  }
  Sl = dot(W, fl) / hbarc;
  return 1 / (1 + exp(Sl)); 
}

/*
 * TransmissionCoeffHW:	TRANSMISSION COEF. WHEN E>Vmax (HILL-WHEELER)	
 */
double fusionSemiClassical::TransmissionCoeffHW(double b_, double rAtInfinity, 
                                                double Energy,
                                                int orbitalAngularMomentum)
{
  double En = Energy;
  int l = orbitalAngularMomentum;
  double mu = reducedMass;
  double  Rmax, Vmax, omegal;
  double  Sl;
  double smallGuess = 1.2; //fm
  double  middle = fabs(0.5 * rAtInfinity);
  double const tolerance = 1.e-7;
  Vmax = -brentMinimise(smallGuess, middle, rAtInfinity,
                        En, l, tolerance, Rmax); 
  
  double const dr = 0.001;
  omegal = fabs(differentiateV(Rmax, dr, En, l));

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
 void fusionSemiClassical::getSfactorAndCrossSection(double r, 
                                                     double rAtInfinity,
                                                     double Energy,
                                                     double &Sfactor, 
                                                     double &crossSection)
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

  double gamow = z1 * z2 * e2 * sqrt( mu / ( 2 * En * hbarc2));
  local = (En > 7.2)? 1:0;

  cout<<"calculating fusion at E = "<<En<<"MeV\n";
  cout<<" l   \t      r_1 \t  r_2 \t    T_l\n";
  for(int l = local; l < T_l.size(); l++)
  {
     a_ = closeToOriginGuess(0.04, En, l);
     b_ = middleGuess(r, En, l);
     
     r1 = findReturnPoint(a_, b_, En, l);
     r2 = findReturnPoint(b_, r, En, l);
             
     if( fabs(r1-r2) < 1e-5 )
     {
       T_l[l] = TransmissionCoeffHW(b_, rAtInfinity, En, l);
     } 
     else
     {
       T_l[l] = TransmissionCoeff(r1, r2, En, l);
     } 
     
     subl[l] = 2 * l + 1;
     
     cout.precision(8);
     cout<<l<<" \t"<<r1<<" \t"<<r2<<" \t"<<T_l[l]<<endl;
  }

  Sfactor = pi * hbarc2 * dot(subl, T_l) * exp(2 * pi * gamow) 
                                                  / (2 * mu * 100); 

  crossSection = pi * hbarc2 * dot(subl, T_l)* 10 / (2 * mu * En);
}
/**************************************************************/
/**************************************************************/


/**************************************************************
    Brent Algorithm: to find minimum of the potentialMenusEnergy()
  TODO: should be a better way to do this, since it is the function
        to minimise is hardwired
 **************************************************************/
double fusionSemiClassical::brentMinimise(const double ax, const double bx, 
                                          const double cx,
                                          const double input4function1, 
                                          const int input4function2, 
                                          const double tol, double &xmin)
{
  const int ITMAX = 100;   
  const double CGOLD = 0.3819660; 
  const double ZEPS = numeric_limits<double>::epsilon()*1.0e-3;
 
  int iter;
  double a, b, d = 0.0, etemp, fu, fv, fw, fx;
  double p, q, r, tol1, tol2, u, v, w, x, xm;
  /********************************************
   *e:=distance moved on the step before last *
   *a,b:=must be in ascending order	      *
   ********************************************/
  double e = 0.0;
  a = (ax < cx? ax:cx);
  b = (ax > cx? ax:cx);
  x = w = v = bx;

  fx = potentialMenusEnergy(x, input4function1, input4function2);
  fw = fv = fx;
                     /*MAIN LOOP STARTS*/

  for(iter = 0; iter < ITMAX; iter++)
  {
    xm = 0.5 * (a+b);
    tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
    if(fabs(x-xm) <= (tol2 - 0.5 * (b-a)))
    {
      xmin = x; 
      return fx;
    }

    if(fabs(e) > tol1)
    {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if(q > 0.0) 
      {
         p = -p;
      }
      q = fabs(q);
      etemp = e;
      e = d;
      if(fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b-x))
      {
	d = CGOLD * (e =(x >= xm)? (a - x):(b - x));
      }
      else
      {
	d = p / q;		     //parabolic step.
	u = x + d;
	if( (u - a) < tol2 || (b - u) < tol2)
	{
	  d = SIGN(tol1,xm - x);
        }
      }
    }
    else
    {
      d = CGOLD * ( e = (x >= xm)? a - x : b - x );
    }
    u = (fabs(d) >= tol1)? x + d: x + SIGN(tol1,d);

    fu = potentialMenusEnergy(u, input4function1, input4function2);

    if(fu <= fx)           //decide what to do w evaluation.
    {
      if(u >= x)
      {
        a = x;
      }
      else 
      {
        b = x;
      }
      shift3(v,w,x,u);
      shift3(fv,fw,fx,fu);
    }
    else
    {
      if(u < x) a = u;
      else b = u;
      if(fu <= fw || w == x)
      {
	 v = w;
	 w = u;
	 fv = fw;
	 fw = fu;
      }
      else if(fu <= fv || v == x || v == w)
	{
	   v = u; 
	  fv = fu;
	}
     }  
   }                  
        
 xmin = x;
 return fx;
}

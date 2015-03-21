/****************************************************************
 *                         folded_potential.cpp                 *
 *                         ASYMMETRIC FUSION                    *
 *DEFINITION OF POTENTIALS FOR THE BARRIER POTENTIAL FOMALISM   *
 *Vf CALCULATED FROM W-S DENSITIES                              *
 *n IS THE SIZE OF THE VECTOR X.                                *
 ****************************************************************/
#include "folded_potential.h"

using namespace std;


/****************************************************************
 *  WOOD-SAXON PARAMETERS                                       *
 ****************************************************************/

double R_o(int baryonNumber)
{
  return 1.31 * pow(baryonNumber,1.0 / 3.0) - 0.84;
}

double diffuseness(int protonNumber)
{
  return (protonNumber == 6)? 0.56:0.58;
}
//This does not replace two_vec function
// This method does not work
double central_density(double Radius, double a, int baryonNumber, 
                       int protonNumber)
{
  double n0 = 0.0;
  double R0 = Radius;
  int A = baryonNumber;
  int Z = protonNumber;

  /* Gaussian integration */
  double ws_denominator = 0.0;
  vector<double> r(200), dr(200);
  GausLeg(0., 20., r, dr);
  for(int i = 0; i < r.size(); i++)
  {
    double local = (r[i] - R0)/a;
    ws_denominator += r[i] * r[i] * dr[i] / (1 + exp(local));
  }  
  ws_denominator *=  4. * pi;
  n0 =  A / ws_denominator;
  return n0;
}

/*
  Woods-Saxon density as a function of radial distance r
   //Should replace rhows function
*/
double ws_density(double Radius, int baryonNumber, int protonNumber, 
                  double a, double density0, double radialDistance)
{
  double rho = 0.0;
  double rho0 = density0;
  double x = (radialDistance - Radius) / a;
  rho = rho0 / (1 + exp(x));  
  return rho;
}

/*
   Folded Potential
  Shall replace V_f function
 */
double folded_potential(double r, double wsRadius1, double wsRadius2,
                        int protonNumber1,int protonNumber2,
                        int baryonNumber1,int baryonNumber2)
{
  double V_folded = 0.0;
  double Ro1 = wsRadius1;
  double Ro2 = wsRadius2;
  int Z1 = protonNumber1;
  int Z2 = protonNumber2;
  int B1 = baryonNumber1;
  int B2 = baryonNumber2;
  double a1 = diffuseness(Z1);
  double a2 = diffuseness(Z2);
  double rho0_1 = central_density(Ro1, a1, B1, Z1); 
  double rho0_2 = central_density(Ro2, a2, B2, Z2); 

  vector<double> cosi(100),w(100),XX1(200),W(200);
  double Vr, integr1, integr2, folded, x1, xx, R1, d_x1_r;
  const double Vo=-456;

  GausLeg(-1.0, 1.0, cosi, w);
  GausLeg(0.0, 30.0, XX1, W);
  integr1 = 0.0;
  integr2 = 0.0;

  Vr = 0;
  for(int j=0;j<XX1.size();j++)
  {
    folded = 0.0;
    x1 = XX1[j];
    for(int i=0;i<cosi.size();i++)
    {
      xx=cosi[i];
      d_x1_r=norma(x1, r, xx);
      folded += w[i] * ws_density(Ro2, B2, Z2, a2, rho0_2, d_x1_r);
    }

    Vr += W[j] * x1 * x1 * ws_density(Ro1, B1, Z1, a1, rho0_1, x1) 
               * folded;
  }

 V_folded = 2 * pi * Vo * Vr;
 return V_folded;
}


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
double central_density(double Radius, double a, int baryonNumber, 
                       int protonNumber)
{
  double R0 = Radius;
  int A = baryonNumber;
  int Z = protonNumber;

  /* Gaussian integration */
  double ws_denominator = 0.0;
  vector<double> r(200), dr(200);
  for(int i = 0; i < r.size(); i++)
  {
    double local = (r[i] - R0)/a;
    ws_denominator += r[i] * r[i] * dr[i] / (1 + exp(local));
  }  
  ws_denominator *=  4. * pi;
  
  return  A / ws_denominator;
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
  /* Should not need these 
  two_vec(Ro1, Z1, integr1);
  two_vec(Ro2, Z2, integr2);
  */
  Vr = 0;
  for(int j=0;j<XX1.size();j++)
  {
    folded = 0;
    x1 = XX1[j];
    for(int i=0;i<cosi.size();i++)
    {
      xx=cosi[i];
      d_x1_r=norma(x1,r,xx);
      folded+=w[i]*ws_density(Ro2, B2, Z2, a2, rho0_2, d_x1_r);
    }
    Vr += W[j] * x1 * x1 *ws_density(Ro1, B1, Z1, a1, rho0_1, x1) 
               * folded;
  }

 V_folded = 2 * pi * Vo * Vr;
 return V_folded;
}

/****************************************************************
****************************************************************

template <class T>
T V_f(T r,T Ro1,T Ro2,int Z1,int Z2,int B1,int B2)
{
vector<T> cosi(100),w(100),XX1(200),W(200);
T Vr,integr1,integr2,folded,x1,xx,R1,d_x1_r;
const double Vo=-456;

GausLeg(-1.0,1.0,cosi,w);
GausLeg(0.0,30.0,XX1,W);
 integr1=0.0;
 integr2=0.0;
 two_vec(Ro1,Z1,integr1);
 two_vec(Ro2,Z2,integr2);
  Vr=0;
    for(int j=0;j<XX1.size();j++)
      {folded=0;
       x1=XX1[j];
       for(int i=0;i<cosi.size();i++)
             {xx=cosi[i];
              d_x1_r=norma(x1,r,xx);
              folded+=w[i]*rhows(B2,Z2,d_x1_r,Ro2,integr2);
             }
       Vr+=W[j]*x1*x1*rhows(B1,Z1,x1,Ro1,integr1)*folded;
      }
 return 2*pi*Vo*Vr;
}
****************************************************************

template<class T,T G(T,int)>
T turnPoint(T a,T b,int L,T &X,int N)
{
 T point,newpoint;
 T c=0.5*(a+b);
 X=c;
 for(int j=0;j<=N;j++)
    {
      point=c;
      (G(c,L)*G(b,L)>0)? b=point:a=point;
      newpoint=0.5*(a+b);
      if(fabs((newpoint-point)/newpoint)<1.0e-12){X=newpoint; break;}
      // if(fabs(G(newpoint,L))<1.0e-10){X=newpoint; break;}
       else c=newpoint;
    }
}
*****************************************************************/


/****************************************************************/
/*                         A_Vfolded.h                          */
/*                         ASYMMETRIC FUSION                    */
/*DEFINITION OF POTENTIALS FOR THE BARRIER POTENTIAL FOMALISM   */
/*Vf CALCULATED FROM W-S DENSITIES                              */
/*n IS THE SIZE OF THE VECTOR X.                                */
/****************************************************************/
#ifndef A_Vfolded
#define A_Vfolded
#include "Gauss_q.h"
using namespace std;

template <class T>
inline T norma(T x1,T r,T cosi){return sqrt(x1*x1+r*r-2*x1*r*cosi);}
/****************************************************************/
template <class T>
//T coulomb(T r,T Ro1,int Z1,int Z2)
T coulomb(T r,T Ro1,T Ro2,int Z1,int Z2)
{const double e2=1.4427713;     //[MeV*fm]
 if(r>(Ro1+Ro2)) return e2*Z1*Z2/r;
 else    return e2*Z1*Z2*(1.5-0.5*r*r/(Ro1*Ro1))/Ro1;//check!!!
}

template <class T>
T V_l(T r,T mr,int l){return 0.5*l*(l+1)*hbarc*hbarc/(mr*r*r);}

/****************************************************************/
/*WOOD-SAXON PARAMETERS                                         */
/****************************************************************/
//template<class T>
double diffuse(int Z){return (Z==6)? 0.56:0.58;}

template <class T>
void two_vec(T Ro,int Z,T &integr)
{vector<T> oneover(200);
 vector<T> X(200),W(200);
 GausLeg(0.0,20.0,X,W);
 double a=diffuse(Z); 
for(int i=0;i<X.size();i++)
    {
      T loc=(X[i]-Ro)/a;
      oneover[i]=X[i]*X[i]/(1+exp(loc));
    }
integr=4*pi*dot(W,oneover);
}


template <class T>
T rhows(int B,int Z,T r,T Ro,T integr)
{
  // const double a=0.58;  //[fm]
 double a=diffuse(Z);
 T local=(r-Ro)/a; 
 T rho0=B/integr;
 return rho0/(1+exp(local));
}

/****************************************************************/

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
/*****************************************************************/
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
/*****************************************************************/

#endif
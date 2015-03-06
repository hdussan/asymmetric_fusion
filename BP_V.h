/****************************************************************/
/*                         BP_V.h                               */
/*DEFINITION OF POTENTIALS FOR THE BARRIER POTENTIAL FOMALISM   */
/*Vf CALCULATED FROM W-S DENSITIES                              */
/*n IS THE SIZE OF THE VECTOR X.                                */
/****************************************************************/
#ifndef BP_V
#define BP_V
#include "Gauss_q.h"

using namespace std;
template <class T>
inline T norma(T x1,T r,T cosi){return sqrt(x1*x1+r*r-2*x1*r*cosi);}

/****************************************************************/
template <class T>
T coulomb(T r,T R0,int Z)
{const double e2=1.4427713;     //[MeV*fm]
 if(r>R0) return e2*Z*Z/r;
 else    return e2*Z*Z*(1.5-0.5*r*r/(R0*R0))/R0;
}

template <class T>
T V_l(T r,T mr,int l){return 0.5*l*(l+1)*hbarc*hbarc/(mr*r*r);}

/****************************************************************/
/*WOOD-SAXON PARAMETERS                                         */
/****************************************************************/

inline double Ro(int B){return 1.31*pow(B,1.0/3.0)-0.84;}

template <class T>
void two_vec(T Ro,T &integr)
{vector<T> oneover(200);
 vector<T> X(200),W(200);
 GausLeg(0.0,20.0,X,W);
 const double a=0.58; 
for(int i=0;i<X.size();i++)
    {
      T loc=(X[i]-Ro)/a;
      oneover[i]=X[i]*X[i]/(1+exp(loc));
    }

 integr=4*pi*dot(W,oneover);
}


template <class T>
T rhows(int B,T r,T Ro,T integr)
{const double a=0.58;  //[fm]
 T local=(r-Ro)/a; 
 T rho0=B/integr;
return rho0/(1+exp(local));
}


/****************************************************************/
template <class T>
T V_f(T r,T R0,int B)
{
vector<T> cosi(200),w(200),XX1(200),W(200);
T Vr,integr,folded,x1,xx,R1;
const double Vo=-456;

GausLeg(-1.0,1.0,cosi,w);
GausLeg(0.0,30.0,XX1,W);
 integr=0;
 two_vec(R0,integr);
  Vr=0;
    for(int j=0;j<XX1.size();j++)
      {folded=0;
       x1=XX1[j];
       for(int i=0;i<cosi.size();i++)
             {xx=cosi[i];
              R1=norma(x1,r,xx);
              folded+=w[i]*rhows(B,R1,R0,integr);
             }
       Vr+=W[j]*x1*x1*rhows(B,x1,R0,integr)*folded;
      }
 return 2*pi*Vo*Vr;
}

/*****************************************************************/

template<class T,T G(T,int)>
T turnPoint(T a,T b,int L,T &X,int N)
{T point,newpoint;

 T c=0.5*(a+b);
 X=c;
 for(int j=0;j<=N;j++)
    {
      point=c;
      (G(c,L)*G(b,L)>0)? b=point:a=point;
       newpoint=0.5*(a+b);
       //       if(fabs(G(newpoint,L))<1.0e-10){X=newpoint; break;}
       if(fabs((newpoint-point)/newpoint)<1.0e-12)
	 {X=newpoint; break;}
       else c=newpoint;
    }
}
/*****************************************************************/

#endif

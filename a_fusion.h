/*
 *  main header file for asymmetric fusion
 */
#ifndef a_fusion_h
#define a_fusion_h
#include <iostream>
#include <fstream>
#include <cstring>
#include "semi_classical.h"
#include "folded_potential.h"
#include "potentials.h"
#include "Derivate.h"
#include "minimal.h"
#include "constants.h"
#include "Gauss_q.h"
using namespace std;

/*  
   Global variables to be removed also
   Everything should be local
 */
int B1, B2, Z1, Z2;
int L;
double Rnot1, Rnot2, r, rInfinity, massR, En, S_E, crossSec;

/*********************************************************
 *menos Effective V (-Veff)                              *
 * To be replaced by function method in semi_classical   *
 *  menosEffectiveV( )                                   *
 *********************************************************/
template<class T>
T _Veff(T r_)
{
 T Vc = coulomb(r_,Rnot1,Rnot2,Z1,Z2);
 T Vl = V_l(r_,massR,L);
 T Vfold = V_f(r_,Rnot1,Rnot2,Z1,Z2,B1,B2);
 T A1 = 2 * (En - Vc)/massR; 
 T A2 = 2 * Vfold / massR;
 T v_2 = speedy2(En, A1, A2);
 return -Vc - Vfold * exp(-4 * v_2) - Vl;
}

/*********************************************************
 *Effective V - En                                       *
 * To be replaced by function method                     * 
 *  PotentialMenusEnergy()                               *
 *   in semi_classical                                   *
 *********************************************************/
template<class T>
T E_eff(T r,int l)
{
 T Vc=coulomb(r,Rnot1,Rnot2,Z1,Z2);
 T Vl=V_l(r,massR,l);
 T Vfold=V_f(r,Rnot1,Rnot2,Z1,Z2,B1,B2);//correct??
 T A1=2*(En-Vc)/massR; 
 T A2=2*Vfold/massR;
 T v_2=speedy2(En, A1, A2);
 return Vc+Vfold*exp(-4*v_2)+Vl-En;
}



/********* TO BE REMOVED ******************/
/********* TO BE REMOVED ****************/
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


/********* TO BE REMOVED ****************/
/**/
template <class T>
void two_vec(T Ro,int Z,T &integr)
{vector<T> oneover(200);
 vector<T> X(200),W(200);
 GausLeg(0.0,20.0,X,W);
 double a = (Z==6)? 0.56:0.58; 
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
 double a = (Z==6)? 0.56:0.58;
 T local=(r-Ro)/a; 
 T rho0=B/integr;
 return rho0/(1+exp(local));
}

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

/********* TO BE REMOVED ***/





#endif

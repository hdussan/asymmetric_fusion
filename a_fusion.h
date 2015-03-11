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
double const hbarc2 = hbarc * hbarc;

/********* TO BE REMOVED ******************/
/********* TO BE REMOVED ****************/
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

/*
 *                       semi_classical.h
 *  FUSION ^Z1+^Z2 REACTIONS
 */

#ifndef semi_classical_h
#define semi_classical_h
#include "constants.h"
#include <cmath>
#include <vector>
using namespace std;
/*
  Reduced mass
 */
double mReduce(int B1, int B2);

/* 
  Find classical return points 
*/
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

/********************************************************
 *   FUNCTIONS	SELF-CONSISTENT VELOCITY CALCULATION	*
 * To be declared  Don't have to be templates    !!!!!  *
 ********************************************************/
//// Self-consistent calculation of the velocity
//// Sao Paulo potential 
template<class T>
inline T f_v2(T A1,T A2,T v2){return A1-A2*exp(-4*v2);}

template<class T>
T speedy2(T En,T A1,T A2)
{T vold2,vnew2,vsonew2,v2;
 v2=0.02;
 inicio:
        vold2=v2;
        vnew2=f_v2(A1,A2,vold2);
        vsonew2=0.9*vnew2+(1-0.9)*vold2;
	if(fabs(vsonew2-vold2)>1.0e-15)
	  { 
           v2=vsonew2;
	   goto inicio;
          }
return v2;
}

#endif



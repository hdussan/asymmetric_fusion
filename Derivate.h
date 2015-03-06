/***************************************************/
/*                  Derivate.h                     */
/*NUMERICAL DIFFERENTIATION                        */
/*GIVEN f(x),xo AND STEP h                         */
/*EVALUATES D_f(x=xo)                              */
/*INCLUDES 4 AND 5 POINT DERIVATIVES               */
/*TESTED                                           */
/*APRIL 18,08                                      */
/***************************************************/
#ifndef Derivate
#define Derivate
#include "methods.h"
using namespace std;
/***************************************************/
template<class T>
T D_f(const T &xo,const T &h,T f(const T))
{return (f(xo+h)-f(xo-h))/(2*h);}

template<class T>
T D2_f(const T &xo,const T &h,T f(const T))
{return (f(xo+h)-2*f(xo)+f(xo-h))/(h*h);}
/***************************************************/
template<class T>
T D_f4(const T &xo,const T &h,T f(const T))
{return (-2*f(xo-h)-3*f(xo)+6*f(xo+h)-f(xo+2*h))/(6*h);}

template<class T>
T D2_f5(const T &xo,const T &h,T f(const T))
{T delta2_f=-f(xo-2*h)+16*f(xo-h)-30*f(xo)+16*f(xo+h)-f(xo+2*h);
 delta2_f=delta2_f/12;
 return delta2_f/(h*h);
}
#endif

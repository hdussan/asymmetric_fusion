/***************************************************
 *                  Derivate.h                     *
 *NUMERICAL DIFFERENTIATION                        *
 *GIVEN f(x),xo AND STEP h                         *
 *EVALUATES D_f(x = xo)                            *
 *INCLUDES 4 AND 5 POINT DERIVATIVES               *
 *TESTED                                           *
 *APRIL 18,08                                      *
 ***************************************************/
#ifndef Derivate
#define Derivate
#include "methods.h"
using namespace std;
/***************************************************/
template<class T>
T D_f(const T &xo,const T &h,T f(const T))
{ 
  return (f(xo+h)-f(xo-h))/(2*h);
}

template<class T>
T D2_f(const T &xo,const T &h,T f(const T))
{
  return (f(xo+h)-2*f(xo)+f(xo-h))/(h*h);
}
/***************************************************/
template<class T>
T D_f4(const T &xo,const T &h,T f(const T))
{ 
  return (-2*f(xo-h)-3*f(xo)+6*f(xo+h)-f(xo+2*h))/(6*h);
}

template<class T>
T D2_f5(const T &xo,const T &h,T f(const T))
{
  T delta2_f = -f(xo-2*h)+16*f(xo-h)-30*f(xo)+16*f(xo+h)-f(xo+2*h);
  delta2_f=delta2_f/12;
  return delta2_f/(h*h);
}
/* Another specialisation to more arguments */
template<class T>
T D2_f5(const T &xo,const T &h,
        const T input4function1, const T input4function2, 
        const T input4function3, 
        const int input4function4, const int input4function5, 
        const int input4function6, const int input4function7,
        const T input4function8, const int input4function9,
        T f(const T, const T, const T, const T, 
            const int, const int, const int, const int,
	    const T, const int))
{
  T f1 = f(xo - 2 * h, input4function1, input4function2, input4function3, 
           input4function4, input4function5, input4function6, 
	   input4function7, input4function8,  input4function9);
  T f2 = f(xo - h, input4function1, input4function2, input4function3, 
           input4function4, input4function5, input4function6, 
	   input4function7, input4function8,  input4function9);
  T f0 = f(xo, input4function1, input4function2, input4function3, 
           input4function4, input4function5, input4function6, 
	   input4function7, input4function8,  input4function9);
  T f3 = f(xo + h, input4function1, input4function2, input4function3, 
           input4function4, input4function5, input4function6, 
	   input4function7, input4function8,  input4function9);
  T f4 = f(xo + 2 * h, input4function1, input4function2, input4function3, 
           input4function4, input4function5, input4function6, 
	   input4function7, input4function8,  input4function9);

  T delta2_f = -f1 + 16 * f2 - 30 * f0 + 16 * f3 - f4;
  delta2_f = delta2_f / 12;
  return delta2_f / (h * h);
}

#endif

/****************************************************************
 *                      minimal.h                               *
 *MINIMISES A FUNCTION f(x)                                     *
 *INPUT POINTS a,b,c                                            *
 *      AND    tol(TOLERANCE>SQRT(MACHINE FLOAT POINT PRECISION)*
 *OUTPUT  xmin AND f(xmin)                                      *
 *MINIMUM WILL BRACKETED BETWEEN a,b,c                          *
 *BRENT'S ALGORITHM                                             *
 *TESTED WORKS PERFECTLY                                        *
 *APRIL 18,08                                                   *
 ****************************************************************/
#ifndef minimal
#define minimal
#include <limits>
#include "methods.h"
using namespace std;

namespace
{
  template<class T>
  inline void shft3(T &a,T &b,T &c,const T d)  
  { 
    a = b; 
    b = c; 
    c = d;
  }
}

/****************************************************************
 *		  BRENT'S METHOD				*
 *ITMAX:=   Maximum allowed number of iterations		*
 *CGOLD:=   Golden ratio					*
 *ZEPS :=   Avoid fractional accuracy for a min at zero		*
 ****************************************************************/
template<class T>
T BRENT(const T ax,const T bx,const T cx,T f(const T),const T tol,T &xmin)
{
  const int ITMAX=100;   
  const T CGOLD=0.3819660; 
  const T ZEPS=numeric_limits<T>::epsilon()*1.0e-3;
 
  int iter;
  T a,b,d=0.0,etemp,fu,fv,fw,fx;
  T p,q,r,tol1,tol2,u,v,w,x,xm;
  /*******************************************/
  /*e:=distance moved on the step before last*/
  /*a,b:=must be in ascending order          */
  /*******************************************/
  T e=0.0;
  a=(ax<cx? ax:cx);
  b=(ax>cx? ax:cx);
  x=w=v=bx;
  fw=fv=fx=f(x);
                     /*MAIN LOOP STARTS*/
  for(iter=0;iter<ITMAX;iter++)
     {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
      if(fabs(x-xm)<=(tol2-0.5*(b-a)))
	{xmin=x; return fx;}

      if(fabs(e)>tol1)
	{r=(x-w)*(fx-fv);
         q=(x-v)*(fx-fw);
	 p=(x-v)*q-(x-w)*r;
	 q=2.0*(q-r);
	 if(q>0.0)p=-p;
	 q=fabs(q);
	 etemp=e;
	 e=d;
	 if(fabs(p)>=fabs(0.5*q*etemp)|| p<=q*(a-x) || p>=q*(b-x))
	    d=CGOLD*(e=(x>=xm? a-x:b-x));
         else{
	   d=p/q;               //parbolic step.
	      u=x+d;
	      if(u-a<tol2 || b-u<tol2)
		d=SIGN(tol1,xm-x);
	     }
        }
      else d=CGOLD*(e=(x>=xm? a-x:b-x));
    
      u=(fabs(d)>=tol1? x+d:x+SIGN(tol1,d));
      fu=f(u);                //one evaluation per iteration.
	 if(fu<=fx)           //decide what to do w evaluation.
            {if(u>=x)a=x;else b=x;
	     shft3(v,w,x,u);
	     shft3(fv,fw,fx,fu);
            }
         else
             {
	       if(u<x)a=u;else b=u;
	       if(fu<=fw || w==x)
	         {
		   v=w;
		   w=u;
		   fv=fw;
		   fw=fu;
		 }
	       else if(fu<=fv ||v==x ||v==w)
	              {v=u; fv=fu;}
             }  
     }                  //END MAIN LOOP.
                        //DONE WITH HOUSEKEEPING
                        //BACK ANOTHER ITERATION.
      
xmin=x;
return fx;
}
/*****************************************************************
 * Specialisation for more arguments of input function           *
 * TODO: Pass arguments in function as an object                 *
 *****************************************************************/
template<class T>
T BRENT(const T ax, const T bx, const T cx,
        const T input4function1, const T input4function2, 
        const T input4function3, 
        const int input4function4, const int input4function5, 
        const int input4function6, const int input4function7,
        const T input4function8, const int input4function9,
        T f(const T, const T, const T, const T, 
            const int, const int, const int, const int,
	    const T, const int),
        const T tol,T &xmin)
{
  const int ITMAX = 100;   
  const T CGOLD = 0.3819660; 
  const T ZEPS = numeric_limits<T>::epsilon()*1.0e-3;
 
  int iter;
  T a, b, d = 0.0, etemp, fu, fv, fw, fx;
  T p, q, r, tol1, tol2, u, v, w, x, xm;
  /********************************************
   *e:=distance moved on the step before last *
   *a,b:=must be in ascending order	      *
   ********************************************/
  T e = 0.0;
  a = (ax < cx? ax:cx);
  b = (ax > cx? ax:cx);
  x = w = v = bx;
  fx = f(x, input4function1, input4function2, input4function3, 
            input4function4, input4function5, input4function6, 
	    input4function7, input4function8,  input4function9);
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
    fu = f(u, input4function1, input4function2, input4function3, 
              input4function4, input4function5, input4function6, 
	      input4function7, input4function8,  input4function9);

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
      shft3(v,w,x,u);
      shft3(fv,fw,fx,fu);
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

#endif

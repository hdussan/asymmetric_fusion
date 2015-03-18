/****************************************************************
 *******************   Some standard methods  *******************
 ****************************************************************/
#ifndef methods
#define methods
#include<iostream>
#include<stdlib.h>
#include<iomanip>
#include<math.h>
#include<iterator>
#include<algorithm>
#include<cmath>
#include<vector> 
#include<errno.h>

using namespace std;

template<class T>
inline void input(char *S,T &x) 
{ 
  cout<<S; cin>>x; 
}

/****************************************************************
 ***************                                     ************
 ***************   useful operations                 ************
 ***************                                     ************
 ****************************************************************/

template<class T>
inline T sign(T x) 
{ 
  return x<0? -1:1;
}

template<class T>
inline const T SIGN(const T &a,const T &b)
{
  return b >= 0? ( a >= 0? a:-a):(a >= 0? -a:a);
}

inline double maxi(double a,double b){return a<b? b:a;}
inline double mini(double a,double b){return a<b? a:b;}

template <class T> 
inline T pyth(T x,T y)
{
  return sqrt(x * x + y * y);
}

template <class T>
inline T norma(T x1,T r,T cosi)
{
  return sqrt(x1 * x1 + r * r - 2 * x1 * r * cosi);
}

/***************************************************************
 *		     INTERPOLATIONS			       *  
 *Linear interpolation between two given points.	       *
 *LAGRANGE						       *
 *NEWTON						       *
 ***************************************************************/
template <class T>
T inter(T x,T x1,T x2,T y1,T y2)
{
  return ( x * (y2 - y1) + (y1 * x2 - y2 * x1)) / (x2 - x1);
}

template <class T>
T linear(const vector<T>& vx,const vector<T>& vy,T x)
{
 int n = vx.size() - 1;
 T b = fabs(x - vx[0]); 
 for(int i = 1; i <= n; i++)
 {
   b = min(b, fabs(x - vx[i]));                  
 }
 T y = 0; 
 for(int i = 0; i <= n; i++) 
 {
   if(b==fabs(x-vx[i])) 
   y=vy[i]*(x-vx[i+1])/(vx[i]-vx[i+1])+vy[i+1]*(x-vx[i])/(vx[i+1]-vx[i]);
 }
 return y;
}

template <class T>
T lagrange(int order,const vector<T>& vx,const vector<T>& vy,T x)
{
 int n=vx.size()-1;
 vector<T> auxX(order),auxY(order);
 T b=fabs(x-vx[0]);
 for(int i=1;i<=n;i++)b=min(b,fabs(x-vx[i]));
 T y=0;
 for(int i=0;i<=n;i++)
   {if(b==fabs(x-vx[i]))
      {for(int j=0;j<order;j++)
	  {if(i==0){auxX[j]=vx[i+j];  auxY[j]=vy[i+j];}
           if(i==n+1){auxX[j]=vx[i-4+j]; auxY[j]=vy[i-4+j];}
           else { auxX[j]=vx[i-1+j]; auxY[j]=vy[i-1+j];} 
          }
      }
   }

 for(int i=0;i<order;i++)    
    {T temp=1;
       for(int j=0;j<=order;j++)
          if(j!=i)temp *=(x-auxX[j])/(auxX[i]-auxX[j]);
       y+=temp*auxY[i];
    }
 return y;
}

template <class T>
T newton(const vector<T>& vx,const vector<T>& vy,T x)
{
 vector<T> b=vy;
 int n=vx.size()-1;
 //find coff. in Newton's form
 for(int j=1;j<=n;j++)
   for(int i=n;i>=j;i--)
     b[i]=(b[i]-b[i-1])/(vx[i]-vx[i-j]);
 //evaluate interpolation polynomial at x
 T u=b[n];
 for(int i=n-1;i>=0;i--)u=b[i]+(x-vx[i])*u;
 return u;
}

template <class T>
void polint(vector<T> &xa,vector<T> &ya,const T x,T &y,T &d_y)
{
 int i,m,ns=0;
 T den,dif,dift,ho,hp,w;
 int n=xa.size();
 vector<T> c(n),d(n);
 dif=fabs(x-xa[0]);
 for(i=0;i<n;i++){          /*find the closest index ns*/
   if((dift=fabs(x-xa[i]))<dif)
       {ns=i;
        dif=dift;
        }
      c[i]=ya[i];
      d[i]=ya[i];
    }

 y=ya[ns--];                 /*initial approx to y      */
 for(m=1;m<n;m++){
    for(i=0;i<(n-m);i++){
        ho=xa[i]-x;
	hp=xa[i+m]-x;
	w=c[i+1]-d[i];
	if((den=ho-hp)==0.0) break;//cerr<<"Error in polint\n";
	den=w/den;
	d[i]=hp*den;          /*the update of c's and d's*/ 
	c[i]=ho*den;
        }
    y+= (d_y=(2*(ns+1)<(n-m) ? c[ns+1] : d[ns--]));
    }
}
/****************************************************************/
/**********    Integration Methods                     **********/
/****************************************************************/
template<class Funct>
double euler(double a,double b,Funct f,int N)
{
  double h=(b-a)/N;
  double sum=f(a);//initialising
  for(int i=1;i<N;i++) 
  {
    sum+=f(a+i*h);
    sum+=f(b);   
  }    //adding last rectangle
  return sum*h;
};

template<class Funct>
double trapezoidal(double a,double b,Funct f,int N)
{
  double h=(b-a)/N;
  double sum=f(a)*0.5;//initialising
  for(int i=1;i<N;i++) {sum+=f(a+i*h);
                        sum+=f(b)*0.5;}    //adding last rectangle
  return sum*h;
};

template<class T,T F(T)>
T gaussian(int n)
{
   T result;
   if(n==2)
   {
     double w[]={1.0,1.0};
     double x[]={-0.57735,0.57735};
     for(int i=0;i<n;i++)result+=w[i]*F(x[i]);
     return result; 
   }
   if(n==5)
   {
     double w[]={0.2369268851,0.2369268851,0.4786286705,0.4786286705,0.5688888889};
     double x[]={-0.9061798459,0.9061798459,-0.5384693101,0.5384693101,0.0};
     for(int i=0;i<n;i++)result+=w[i]*F(x[i]);
     return result; 
   }
}

/****************************************************************
 *****************		*********************************
 ***************** Finding roots ********************************
 *****************		*********************************
 *****************************************************************/

template<class T,T F(T,T,T)>
  T Zero(T a,T b,T Y,T Z,T &X,int N)
    {T point,newpoint;
     T c=(a+b)/2;
     X=c;
     for(int j=0;j<=N;j++)
       {point=c;
       (F(Y,Z,c)*F(Y,Z,b)>0)? b=point:a=point;
       newpoint=(a+b)/2;
       if(fabs(F(Y,Z,newpoint))<1e-10){X=newpoint;break;}
	else  c=newpoint;
	}
    }

/* Specialised vesion of the bisection method */
template<class T,T F(T,int)>
T Zero(T a,T b,int l,T &X,int N)
{
  T point, newpoint;
  T c = (a + b) / 2;
  X = c;
  for(int j=0;j<=N;j++)
  {
    point=c;
    (F(c, l) * F(b, l) > 0)? b = point: a = point;
    newpoint = (a + b) / 2;
    if(fabs(F(newpoint, l)) < 1e-10)
    {
      X=newpoint;
      break;
    }
    else  c=newpoint;
   }
 }

/* 
    Even more specialised vesion of the bisection method 
    function with more arguments...
*/
template<class T,T F(T, T, T, T, int, int, int, int, T, int)>
  T Zero(T a,T b, T E, T R1, T R2, int Z1, int Z2, int B1, int B2, T Mu,
         int l,T &X,int N)
{
  T point, newpoint;
  T c = (a + b) / 2;
  X = c;
  for(int j = 0;j <= N; j++)
  {
    point = c;
    (F(c, E, R1, R2, Z1, Z2, B1, B2, Mu, l) 
        * F(b, E, R1, R2, Z1, Z2, B1, B2, Mu, l) > 0)? b = point: a = point;
    newpoint = (a + b) / 2;
    if(fabs(F(newpoint, E, R1, R2, Z1, Z2, B1, B2, Mu, l)) < 1e-10)
    {
      X=newpoint;
      break;
    }
    else  c = newpoint;
   }
 }


/****************************************************************
 **********		    VECTORS		       **********
 ****************************************************************/
template<class T>
T dot(const vector<T>& u,const vector<T>& v)
{
 T dotprod = 0;
 int n = u.size();
 for(int i=0; i<n; i++)
 {
   dotprod += u[i] * v[i];
 }
 return dotprod;
}

#endif

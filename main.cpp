/*
 *                        main.cpp
 *FUSION ^Z1+^Z2 REACTIONS
 *NEW CODE FOR ASYMMETRIC FUSION.
 *A_Vfolded.h GETS VFOLDED FROM W-S.
 *BUILDING NOV 15,07.
 */
/********************************************************/

#include "a_fusion.h"
using namespace std;

/*********************************************************
 *  midGuess & smallGuess				 *
 *CHECK FOR GOOD GUESSES FOR FINDING THE TURNING POINTS	 *
 *INPUTS r_inf $ r_o					 *
 * to be replaced by middleGuess()  in semiclassical     *
/*********************************************************/
template<class T>
T midGuess(T r_inf,int l)
{T a=0;
  if(E_eff(r_inf, l)>0) a=r_inf;
  else{
       here:
       r_inf-=0.2;
       if(E_eff(r_inf, l)>0)
	 a=r_inf;
       else goto here;
      } 
 return a;
}

/*
template<class T>
T smallGuess(T r_o,int l)
{T c=0;

 if(E_eff(r_o,l)<0) c=r_o;
 else{
      aqui:
      r_o+=0.02;
      if(E_eff(r_o,l)<0)
	 c=r_o;
       else goto aqui;
      } 
 return c;
}*/

/*****************************************************/
template<class T>
T smallGuess(T r_o,int l)
{
  int const maxIterations = 50;
  T c=0, cplus, cminos;
  Zero<T, E_eff>(r_o, 9.,l, c, maxIterations);
  cplus=c+2*r_o;
  cminos=c-2*r_o;
  return (E_eff(cplus, l)<0)? cplus:fabs(cminos); 
}

/********************************************************
 *    PROBABILITY OF PENETRATING THE BARRIER POTENTIAL	*
 *GIVEN TURNING POINTS r1 AND r2 GETS Tl		*
 *Tl   TRANSMISSION COEF.				*
 *Sl   PENETRATION PROBABILITY				*
 *Tl_HW	TRANSMISSION COEF. WHEN E>Vmax (HILL-WHEELER)	*
 *********************************************************/
template<class T>
T Tl(T r2,T r1,int l)
{
 T Sl=0.0;
 vector<T> X(200),W(200),fl(200);
 GausLeg(r2,r1,X,W);
 for(int i=0;i<X.size();i++)
   {T r_aux=X[i];
    fl[i]=sqrt(8*massR*fabs(E_eff(r_aux,l)));
   }
 Sl=dot(W,fl)/hbarc;
 return 1/(1+exp(Sl)); 
}

template<class T>
T Tl_HW(T b_)
{T Rmax,Vmax,omegal;
 T Sl;
 T middle=fabs(0.5*rInfinity);
 Vmax=-BRENT(1.2,middle,rInfinity,_Veff,1e-7,Rmax); 
 omegal=fabs(D2_f5(Rmax,0.001,_Veff));
 omegal=hbarc*sqrt(omegal/massR);
 Sl=2*pi*(Vmax-En)/omegal;
 return 1/(1+exp(Sl));
}

/*********************************************************/
/*S_fact   ASTROPHYSICAL S-FACTOR                        */
/*********************************************************/
template<class T>
void S_fact(T r,T En,T &S_E,T &crossSec)
{
 T a_,b_,r1,r2;
 int local;
 S_E=0.0;
 vector<T> T_l(21),subl(21);
 cout<<"calculating for E= "<<En<<"MeV\n";
 T gamow=Z1*Z2*e2*sqrt(massR/(2*En*hbarc2));
 (En>7.2)?local=1:local=0;
for(int l_=local;l_<T_l.size();l_++)
    {
      a_=smallGuess(0.04,l_);
      b_=midGuess(r,l_);
      /* VERY SLOW */
      Zero<T,E_eff>(a_,b_,l_,r2,100);
      Zero<T,E_eff>(b_,r,l_,r1,100);

      T deltar=fabs(r1-r2);
      L=l_;
      T_l[l_]=(deltar<1e-5? Tl_HW(b_):Tl(r2,r1,l_));      
      subl[l_]=2*l_+1;
      cout.precision(8);
      cout<<l_<<" \t"<<r2<<" \t"<<r1<<" \t"<<T_l[l_]<<endl;
    }
 S_E=pi*hbarc2*dot(subl,T_l)*exp(2*pi*gamow)/(2*massR*100); 
 crossSec=pi*hbarc2*dot(subl,T_l)*10/(2*massR*En);
}

/*********************************************************/
int main()
{
 ofstream Otro;
 char afuera[15];
 /** Commented for faster testing
 cout<<"  ASYMMETRIC FUSION CALCULATION\n";
 cout<<" A1=? and Z1=?\n";
 cin>>B1;    cin>>Z1; 
 cout<<" A2=? and Z2=?\n";
 cin>>B2;    cin>>Z2;  
 **/
 B1 = 18;
 Z1 = 8;
 B2 =18;
 Z2 = 8;
 Rnot1 = R_o(B1);
 Rnot2 = R_o(B2);
 massR = mReduce(B1, B2);

 cout<<"W-S parameter Ro1="<<Rnot1<<" fm\n";
 cout<<"W-S parameter Ro2="<<Rnot2<<" fm\n";
 cout<<"reduced mass mu="<<massR<<" MeV\n";
 /**
 cout<<"Enter guessed r infinity ="; 
 cin>>r; 
 cout<<"fm";
**/
 r = 50;
 rInfinity = 30;
 cout<<"\n";
 
 
 En =1.0;
 S_E=0.0; crossSec=0.0;
 sprintf(afuera,"%s%d%s%d%s","S_",int(B1),"_",int(B2),".dat");
  Otro.open(afuera,ios::out);
 for(int i=0;i<=12;i++)
   {
     S_fact(r,En,S_E,crossSec);
     Otro<<En<<"  \t"<<S_E<<" \t"<<crossSec<<endl;
     cout<<"energy \t"<<"s-factor \t"<<"cross section\n";
     cout<<En<<"  \t"<<S_E<<" \t"<<crossSec<<endl;
     En+=0.5;
   } 
 Otro.close();
 
}

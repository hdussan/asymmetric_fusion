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

int main()
{
 ofstream Otro;
 char afuera[15];
 double energyStep = 0.5;
 double r, rInfinity;
 int Z1, Z2, B1, B2;

 cout<<"  ASYMMETRIC FUSION CALCULATION\n";
 
 cout<<" A1=? and Z1=?\n";
 cin>>B1;    cin>>Z1; 
 cout<<" A2=? and Z2=?\n";
 cin>>B2;    cin>>Z2;  
 
 double Rnot1 = R_o(B1);
 double Rnot2 = R_o(B2);
 double massR = mReduce(B1, B2);

 cout<<"W-S parameter Ro1="<<Rnot1<<" fm\n";
 cout<<"W-S parameter Ro2="<<Rnot2<<" fm\n";
 cout<<"reduced mass mu="<<massR<<" MeV\n";
 
 fusionSemiClassical afusion(Rnot1, Z1, B1, Rnot2,  Z2, B2, massR);

 r = 50;
 rInfinity = 40;
 cout<<"\n";
 
 
 double En =1.0;
 double S_E=0.0, crossSec=0.0;
 sprintf(afuera,"%s%d%s%d%s","S_",int(B1),"_",int(B2),"old.dat");
 Otro.open(afuera,ios::out);
 
 for(int i = 0; i <= 10 ; i++)
   {
     afusion.getSfactorAndCrossSection(r, rInfinity, En,
                                       S_E, crossSec);

     Otro<<En<<"  \t"<<S_E<<" \t"<<crossSec<<endl;
     cout<<"energy \t"<<"s-factor \t"<<"cross section\n";
     cout<<En<<"  \t"<<S_E<<" \t"<<crossSec<<endl;
     En += energyStep;
   } 
 
 Otro.close();
 
}

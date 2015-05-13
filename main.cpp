/*
 *                        main.cpp
 * FUSION ^Z1 A1 + ^Z2 A2 REACTIONS FOR ASTROPHYSICAL
 * CALCULATIONS.
 * BUILDING NOV 15,07.
 *
/********************************************************/

#include "a_fusion.h"
using namespace std;

int main()
{  
  int Z1, Z2, B1, B2;
  double Rnot1, Rnot2, massR;  
  cout<<"  ASYMMETRIC FUSION CALCULATION\n";
   
  cout<<" A1=? and Z1=?\n";
  cin>>B1;    cin>>Z1; 
  cout<<" A2=? and Z2=?\n";
  cin>>B2;    cin>>Z2;  
   
   Rnot1 = R_o(B1);
   Rnot2 = R_o(B2);
   massR = mReduce(B1, B2);
   
  cout<<"W-S parameter Ro1="<<Rnot1<<" fm\n";
  cout<<"W-S parameter Ro2="<<Rnot2<<" fm\n";
  cout<<"reduced mass mu="<<massR<<" MeV\n";
   
  fusionSemiClassical afusion(Rnot1, Z1, B1, Rnot2,  Z2, B2, massR);
  //Energy range     
  double minEnergy = 1.0; 
  double maxEnergy = 10.;
  int numberData = 10;
  afusion.saveAstrophysicalSfactor(minEnergy, maxEnergy, numberData);

  return 0;
 
}

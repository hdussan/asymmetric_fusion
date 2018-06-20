/*
 * Helber Dussan
 * 2007
 *                      
 * FUSION ^Z1 A1 + ^Z2 A2 REACTIONS FOR ASTROPHYSICAL CALCULATIONS.
 *
 *
 * Updated to run in parallel in the range of energies 
 * specified
 */

#include "a_fusion.h"
using namespace std;

int main()
{  
  int Z1, Z2, B1, B2;
  double Rnot1, Rnot2, massR;  
  cout<<"  ASYMMETRIC FUSION CALCULATION\n";    
  cin>>B1;    
  cin>>Z1; 

  cin>>B2;    
  cin>>Z2;  
   
  Rnot1 = R_o(B1);
  Rnot2 = R_o(B2);
  massR = mReduce(B1, B2);
   
  cout<<"W-S parameter Ro1="<<Rnot1<<" fm\n";
  cout<<"W-S parameter Ro2="<<Rnot2<<" fm\n";
  cout<<"reduced mass mu="<<massR<<" MeV\n";
   
  fusionSemiClassical afusion(Rnot1, Z1, B1, Rnot2,  Z2, B2, massR);
 
  //Energy range in MeV
  
  double minEnergy = 0.5; 
  double maxEnergy = 10.5;
  int numberData = 20;
  vector<double> energy,Sfactor,crossSection;
  solveFusion(afusion, minEnergy,maxEnergy, numberData, energy, Sfactor, crossSection);
  
  return 0;
}

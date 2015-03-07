/****************************************************************
                    Some standard constants  
 ****************************************************************/
#ifndef constants_h
#define constants_h
#include<stdlib.h>
#include<algorithm>
#include<cmath>

using namespace std;

double const pi = 3.1415926535;
double const e_ = 2.7182818284;
double const ln_2 = 0.6931471805;

/*   physical constants  */
double const M     = 939;     //baryon mass [MeV]
double const me    = 0.511;   //electron mass
double const hbarc = 197.33;  //conversion constant
double const alpha = 0.091725;//(e^2/(4pi*hbarc)) Structure fine constant
double const e2 = 1.4427713;     //[MeV*fm]  = 1.440292?
double const Uma = 931.5;//494013;//Atomic unit mass[MeV]
double const lp = 1.793;  //Proton's anomalous magnetic moment 
double const ln=-1.913;   //Neutron's anomalous magnetic moment 


#endif

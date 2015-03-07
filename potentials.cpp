/*
 *        Potentials.cpp
 *  Definitions of all the potentials
 *  for the semiclassical scattering
 *  of two nuclei
 */
 #include "Potentials.h"
 using namespace std;
/*
  Coulomb Potential between 2 nuclei
  inputs:
    r: radial distance [fm]
    Ro1: radius from charge dist. of nucleus 1 [fm]
    Z1:  Proton number of nucleus 1
    Ro2: radius from charge dist. of nucleus 2 [fm]
    Z2:  Proton number of nucleus 2
  output:
   Coulomb potential in MeV
 */
 double coulomb(double r, double Ro1, double Ro2, int Z1,int Z2)
{
  if(r > (Ro1 + Ro2))
  { 
    return e2 * Z1 * Z2 / r;
  }
  else
  {
    return e2 * Z1 * Z2 * (1.5 - 0.5 * r * r / (Ro1 * Ro1)) / Ro1;
  }
}

/*
   Centrifugal potential
   inputs:
     r: radial distance [fm]
     reducedMass: reduced mass [MeV]
     l: angular momentum quantum number
   output:
     centrifugal potential [MeV]
 */
double V_l(double r, double reducedMass, int l)
{
  return 0.5 * l * (l + 1) * hbarc * hbarc / (reducedMass * r * r);
}

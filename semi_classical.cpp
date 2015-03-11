/*
 *                       semi_classical.cpp
 *   FUSION ^Z1+^Z2 REACTIONS
 *
 ********************************************************/

#include "semi_classical.h"
using namespace std;

/*     
  Reduced mass
 */
double mReduce(int B1,  int B2)
{
  return B1 * B2 * Uma /(B1 + B2);
}

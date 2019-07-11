/** This is a CPP Header file **/
#include<cstdlib>
#include<cmath>
#include <limits>
typedef struct threadLocalDataBoxMuller {
  double z1;
  bool generate;
} tldBoxMuller;


/****
 * Box Muller method to generate random numbers 
 * in Gaussian (Normal) distribution. 
 * See https://en.wikipedia.org/wiki/Box-Muller_transform
 * for more details.
 ****/


double generateGaussiandistibution(double mu, double sigma, tldBoxMuller& d) {

  
  static const double epsilon = std::numeric_limits<double>::min();
  static const double two_pi = 2.0*3.14159265358979323846;

  d.generate = !d.generate;

  if (!d.generate)
    return d.z1 * sigma + mu;

  double u1, u2;
  do
    {
      u1 = rand() * (1.0 / RAND_MAX);
      u2 = rand() * (1.0 / RAND_MAX);
    }
  while ( u1 <= epsilon );

  double z0;
  z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
  d.z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
  return z0 * sigma + mu;
}

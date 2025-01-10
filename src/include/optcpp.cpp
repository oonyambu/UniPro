
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <stdlib.h>
#include <stdio.h>


double branin(double xx[2]) {

  double x1 = xx[0];
  double x2 = xx[1];

  // Define default values for parameters
  double a = 1.0;
  double b = 5.1 / (4 * M_PI * M_PI);
  double c = 5 / M_PI;
  double r = 6.0;
  double s = 10.0;
  double t = 1 / (8 * M_PI);

  // Compute intermediate terms
  double term1 = a * pow(x2 - b * pow(x1, 2) + c * x1 - r, 2);
  double term2 = s * (1 - t) * cos(x1);

  // Compute the final result
  double y = term1 + term2 + s;


  return y;
}

void braningr(int n, double xx[2], double y[2], void * args){
  double x1 = xx[0];
  double x2 = xx[1];

  // Define default values for parameters
  double a = 1.0;
  double b = 5.1 / (4 * M_PI * M_PI);
  double c = 5 / M_PI;
  double r = 6.0;
  double s = 10.0;
  double t = 1 / (8 * M_PI);
  y[0] = 2*a*(x2 - b*x1*x1 + c*x2 - r) *(c-2*b*x1) - s * (1-t)*sin(x1);
  y[1] = 2*a*(x2 - b*x1*x1 + c*x2 - r);

}

//typedef double (*optimfn)(int, double *, void *);
//typedef void (*optimgr)(int, double *, double *, void *);


// int main(){
//   double a[2]={1.};
//   printf("%f\n", branin(a));
//
// }

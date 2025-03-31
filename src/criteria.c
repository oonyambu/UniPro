#include "criteria.h"

double maxpro(int * data, paramsPtr p){
 double dsq = 0;
 int len=0;
  for(int i = 0; i < p->n - 1; ++i){
    for(int j = i + 1; j < p->n; ++j, ++len){
      double row = 1;
      for(int k = 0; k < p->m; ++k){
        row *= data[i + k * p->n] - data[j + k * p->n];
      	if (row == 0) return INFINITY;
      }
      dsq += 1/(row * row);
    }
  }
  return pow(dsq/len, 1./p->m);

}


double unipro(int * data, paramsPtr p){
  double dsq = 0, colsum = 0;
  double *col = (double*)calloc(p->n, sizeof(double));
  for(int i = 0; i < p->n; i++){
    for(int j = i + 1; j < p->n; j++){
      double row = 0;
      for(int k = 0; k < p->m; k++)
        row += abs(data[i + k * p->n] - data[j + k * p->n]);
      col[i] +=row;
      col[j] +=row;
      dsq += row * row;
    }
    colsum += col[i]*col[i];
  }
  free(col);
  return (2*(dsq - colsum/p->n)/p->denom_g + p->C)*1000;
}



double maximinLHD(int * data, paramsPtr p){
 double dsq = 0;
 int len=0;
  for(int i = 0; i < p->n - 1; ++i){
    for(int j = i + 1; j < p->n; ++j, ++len){
      double row = 0;
      for(int k = 0; k < p->m; ++k){
        double d =  data[i + k * p->n] - data[j + k * p->n];
      	row += d*d;
      }
      dsq += 1/pow(row, 0.5*p->r);
    }
  }
  return pow(dsq/len, 1./p->r);

}

// #include <stdio.h>


// #define n 5
// #define m 2
// #define s 5
// int main(int argc, char const *argv[])
// {
// 	int a[10] = {0, 1, 2, 3, 4, 1, 3, 0, 4, 2};
//     paramsPtr p = initializeParams(n, m, s, 0, 0, 0.,0., 0., 1, 1ll);
//     criteria funs[2] = {unipro, maxpro};

// 	printf("%f\n", funs[1](a, p));
// 	return 0;
// }

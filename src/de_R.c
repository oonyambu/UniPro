#include <Rinternals.h>
#include <R.h>
#include <R_ext/Random.h>
#include "de.h"



SEXP phi_Call(SEXP R_X, SEXP R_s){
  SEXP dim = getAttrib(R_X, R_DimSymbol);
  const int n = INTEGER(dim)[0], m = INTEGER(dim)[1];
  const int s = INTEGER(R_s)[0];
  int * X = INTEGER(R_X);
  SEXP R_out = PROTECT(allocVector(REALSXP, 1));
  double p = pow(s, 4);
  double numerator = 4*(5.0*m - 2)*p + 30.0*(3*m-5)*s*s + 15*m + 33;
  double denom = 720.0*(m - 1)*p;
  double C = numerator/ denom + (1.0 + pow(-1,s))/(64.0*p);
  double denom_g = 4.0*m * (m - 1) * n*n*s*s;

  double *out = REAL(R_out);
  *out = phi2D(X, n, m, s, denom_g, C);
  UNPROTECT(1);
  return R_out;
}


SEXP DE_Call(SEXP R_n, SEXP R_m, SEXP R_s, SEXP R_NP,
        SEXP R_itermax, SEXP R_pMut, SEXP R_pCR,
        SEXP R_pGBest, SEXP R_replicates){

  struct timespec start, end;
  clock_gettime(CLOCK_MONOTONIC, &start);
  int n = INTEGER(R_n)[0];
  int m = INTEGER(R_m)[0];
  int s = INTEGER(R_s)[0];
  int NP = INTEGER(R_NP)[0];
  int itermax = INTEGER(R_itermax)[0];
  int replicates = INTEGER(R_replicates)[0];
  double pMut = REAL(R_pMut)[0];
  double pCR = REAL(R_pCR)[0];
  double pGBest = REAL(R_pGBest)[0];

  SEXP R_X = PROTECT(allocVector(INTSXP, n*m));
  SEXP R_opt_val = PROTECT(allocVector(REALSXP, replicates));
  SEXP R_time_taken = PROTECT(allocVector(REALSXP, 1));

  double *opt_val = REAL(R_opt_val);
  int * X = INTEGER(R_X);
  double *time_taken = REAL(R_time_taken);

  int* array = (int*)malloc(n*m*NP * sizeof(int));
  int* potential = (int*)malloc(n*m*NP * sizeof(int));
  double *vals = (double*)malloc(NP *sizeof(double));
  int bestIndex;

  double p = pow(s, 4);
  double numerator = 4*(5.0*m - 2)*p + 30.0*(3*m-5)*s*s + 15*m + 33;
  double denom = 720.0*(m - 1)*p;
  double C = numerator/ denom + (1.0 + pow(-1,s))/(64.0*p);
  double denom_g = 4.0*m * (m - 1) * n*n*s*s;

  //#pragma omp parallel for private(it)
  for(int reps = 0; reps < replicates; reps++){
    gen_design(array, n, m, s, NP, vals, &bestIndex, denom_g, C);

    for(int it = 0; it< itermax; it++){
      for(int k = 0; k< NP; k++){
        int gl = rand_unif() < pGBest;
        for(int j = 0; j < m; j++){
          for (int i = 0; i < n; i++){
            potential[i + n*(j + m*k)] = gl? array[i + n*(j + m*bestIndex)] : array[i + n*(j + m*k)];
          }
          if (rand_unif() < pMut){
            int v[2];
            sample2(n, v);
            swap(&potential[v[0] + n*(j + m*k)], &potential[v[1] + n*(j + m*k)]);
          }
          int CR = rand_int(m);
          double pcr = rand_unif();
          for(int i = 0; i<n; i++){
            if( (pcr > pCR) & (j != CR)){
              potential[i + n*(j + m*k)] = array[i + n*(j + m*k)];
            }
          }
        }
        double phi_val = phi2D(potential + k*n*m, n, m, s, denom_g, C);
        if(phi_val < vals[k]) {
          vals[k] = phi_val;
          for (int i=0; i<n; i++){
            for(int j=0; j<m; j++) array[i + n*(j + m*k)] = potential[i + n*(j + m*k)];
          }
        }
        if(phi_val < vals[bestIndex]) bestIndex = k;
      }
    }
      opt_val[reps] = vals[bestIndex];
  }


  clock_gettime(CLOCK_MONOTONIC, &end);
  *time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

  int len = replicates == 1;
  SEXP R_out = PROTECT(allocVector(VECSXP, 2 + len));

  SET_VECTOR_ELT(R_out, 0, R_opt_val);
  SET_VECTOR_ELT(R_out, 1, R_time_taken);

  SEXP names = PROTECT(allocVector(STRSXP, 2 + len));
  SET_STRING_ELT(names, 0, mkChar("opt_val"));
  SET_STRING_ELT(names, 1, mkChar("time_taken"));

  if(len){

      for(int i=0;i<n;i++)
        for(int j=0; j<m; j++)
          X[i + n*j] = array[i + n*(j + m*bestIndex)];

    SET_STRING_ELT(names, 2, mkChar("X"));
    SEXP dims = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dims)[0] = n; // Number of rows
    INTEGER(dims)[1] = m; // Number of columns
    setAttrib(R_X, R_DimSymbol, dims); // Set dimension attribute
    SET_VECTOR_ELT(R_out, 2, R_X);
  }
  setAttrib(R_out, R_NamesSymbol, names);
  free(array);
  free(potential);
  free(vals);

  UNPROTECT(6);

  return R_out;
}



// int main(){
//   int n = 30, m = 3, s = 30, itermax = 1000, NP=100, replicates=10, X[90];
//   double pMut = 0.2, pCR=0.3, pGBest = 0.9, opt_val[10], time_taken;

//   DE_C(&n, &m, &s, &NP, &itermax,&pMut, &pCR, &pGBest, X,
//        opt_val, &time_taken, &replicates);
//   //DE_C(n, m, s, NP, itermax,pMut, pCR, pGBest, X,
//    //     opt_val, &time_taken, replicates);
//   for(int i=0; i<replicates;i++)printf("%f\n", opt_val[i]);
// }

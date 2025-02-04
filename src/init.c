#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "de.h"
#include "criteria.h"


static criteria PHI_FUNS[3] = {unipro, maxpro, maximinLHD};

SEXP phi(SEXP R_X, SEXP R_s, SEXP R_r, SEXP R_METHOD);

SEXP DE(SEXP R_n, SEXP R_m, SEXP R_s, SEXP R_NP,
        SEXP R_itermax, SEXP R_pMut, SEXP R_pCR,
        SEXP R_pGBest, SEXP R_replicates, SEXP R_seed,
        SEXP R_cores, SEXP R_METHOD, SEXP R_r);

SEXP detectCores(){
  SEXP R_out = PROTECT(allocVector(INTSXP, 1));
  INTEGER(R_out)[0] =  omp_get_max_threads();
  UNPROTECT(1);
  return R_out;
}

// Register the function
static const R_CallMethodDef CallMethods[] = {
  {"phi", (DL_FUNC) &phi, 3},
  {"DE", (DL_FUNC) &DE, 12},
  {"detectCores", (DL_FUNC) &detectCores, 0},
  {NULL, NULL, 0}
};

void R_init_RcppProgressExample(DllInfo *Info) {
  R_registerRoutines(Info, NULL, CallMethods, NULL, NULL);
  R_useDynamicSymbols(Info, FALSE);
}


SEXP phi(SEXP R_X, SEXP R_s,SEXP R_r,  SEXP R_method){
  if (!isMatrix(R_X)) {
    Rf_error("Input must be a matrix.");
  }
  R_len_t n = nrows(R_X);
  R_len_t m = ncols(R_X);
  int s = INTEGER(R_s)[0];
  int r = INTEGER(R_r)[0];
  int * data = INTEGER(R_X);
  SEXP R_out = PROTECT(allocVector(REALSXP, 1));

  criteria PHI = PHI_FUNS[INTEGER(R_method)[0] - 1];

  REAL(R_out)[0] = PHI(data,
       initializeParams(n, m, s, 0, 0, 0.,0., 0., 1, 1ll, r));
  UNPROTECT(1);
  return R_out;
}


SEXP DE(SEXP R_n, SEXP R_m, SEXP R_s, SEXP R_NP,
            SEXP R_itermax, SEXP R_pMut, SEXP R_pCR,
            SEXP R_pGBest, SEXP R_replicates, SEXP R_seed,
            SEXP R_cores, SEXP R_method, SEXP R_r){


  int n = INTEGER(R_n)[0];
  int m = INTEGER(R_m)[0];
  int s = INTEGER(R_s)[0];
  int r = INTEGER(R_r)[0];
  int NP = INTEGER(R_NP)[0];
  int itermax = INTEGER(R_itermax)[0];
  int replicates = INTEGER(R_replicates)[0];
  int cores = INTEGER(R_cores)[0];
  int seed = INTEGER(R_seed)[0];
  double pMut = REAL(R_pMut)[0];
  double pCR = REAL(R_pCR)[0];
  double pGBest = REAL(R_pGBest)[0];


  SEXP R_bestX = PROTECT(allocVector(INTSXP, n*m));
  SEXP R_phiValues = PROTECT(allocVector(REALSXP, replicates));
  SEXP R_timeTaken = PROTECT(allocVector(REALSXP, 1));
  SEXP R_out = PROTECT(allocVector(VECSXP, 3));
  SEXP dims = PROTECT(allocVector(INTSXP, 2));
  SEXP names = PROTECT(allocVector(STRSXP, 3));
  SEXP className = PROTECT(allocVector(STRSXP, 1));




  double *phiValues = REAL(R_phiValues);
  int * bestX = INTEGER(R_bestX);
  double *timeTaken = REAL(R_timeTaken);

  int method = INTEGER(R_method)[0] - 1;
  printf("%s\n", method == 0? "UniPro" : (method == 1? "MaxPro":"maximinLHD"));
  criteria phi = PHI_FUNS[method];
  DE_CC(n, m, s, NP, itermax, pMut, pCR, pGBest,
        replicates, seed, phiValues, timeTaken,
        bestX, cores, phi, r);

  INTEGER(dims)[0] = n;
  INTEGER(dims)[1] = m;
  setAttrib(R_bestX, R_DimSymbol, dims);

  SET_VECTOR_ELT(R_out, 0, R_phiValues);
  SET_VECTOR_ELT(R_out, 1, R_timeTaken);
  SET_VECTOR_ELT(R_out, 2, R_bestX);


  SET_STRING_ELT(names, 0, mkChar("measure"));
  SET_STRING_ELT(names, 1, mkChar("timeTaken"));
  SET_STRING_ELT(names, 2, mkChar("Design"));
  setAttrib(R_out, R_NamesSymbol, names);


  SET_STRING_ELT(className, 0, mkChar("DE"));
  classgets(R_out, className);


  UNPROTECT(7);

  return R_out;
}

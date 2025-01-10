#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Declare Fortran subroutine */
void F77_NAME(de_fortran)(const int *n, const int *m, int *s,
        const int *NP, const int *itermax,
        const double *pMut, const double *pCR,
        const double *pGBest, int * X, double* val,
        double* time_taken, const int * replicates);

void F77_NAME(phi_fortran)(const int *X, const int *n, const int *m,
              int *s, double * val);

/* Declare .C functions */
void phi_C(const int *X, const int *n, const int *m, const int *s, double*val);
void DE_C(int *n, int *m, int *s, int *NP, int *itermax,
          double *pMut, double *pCR, double *pGBest, int *X,
          double * opt_val, double *time_taken, int *replicates);


/*Declare .Call functions */
SEXP phi_Call(SEXP R_X, SEXP R_s);
SEXP DE_Call(SEXP R_n, SEXP R_m, SEXP R_s, SEXP R_NP,
             SEXP R_itermax, SEXP R_pMut, SEXP R_pCR,
             SEXP R_pGBest, SEXP R_replicates);

// SEXP optim2(SEXP par, SEXP fn, SEXP gr, SEXP method,
//             SEXP slower, SEXP supper,   SEXP options);
// //SEXP bran_min(SEXP par, SEXP options);
//
// /* Declare .External functions */
// SEXP optim(SEXP call, SEXP op, SEXP args, SEXP rho);



static const R_FortranMethodDef FortranEntries[] = {
    {"de_fortran", (DL_FUNC) &F77_NAME(de_fortran), 12},  // 2 = number of arguments
    {"phi_fortran", (DL_FUNC) &F77_NAME(phi_fortran), 5},
    {NULL, NULL, 0}
};

static const R_CMethodDef CEntries[] = {
    {"DE_C", (DL_FUNC) &DE_C, 12},  // Name, pointer, # of args
    {"phi_C", (DL_FUNC) &phi_C, 5},
    {NULL, NULL, 0} // Sentinel
};

static const R_CallMethodDef CallEntries[] = {
  {"DE_Call", (DL_FUNC) &DE_Call, 9},  // Name, pointer, # of args
  //{"optim_call", (DL_FUNC) &optim2, 7},
  {"phi_Call", (DL_FUNC) &phi_Call, 2},
  //{"bran_min_call", (DL_FUNC) &bran_min, 2},
  {NULL, NULL, 0} // Sentinel
};

/*
static const R_ExternalMethodDef ExtEntries[] = {
  {"optim_ex", (DL_FUNC) &optim, 7},
  //{"bran_min_ex", (DL_FUNC) &bran_min, 2},
  {NULL, NULL, 0} // Sentinel
};
*/
void R_init_UniPro(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, CallEntries, FortranEntries, NULL);
                       //FortranEntries, ExtEntries);
    R_useDynamicSymbols(dll, FALSE);
}

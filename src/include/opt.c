#include <R.h>
#include <Rinternals.h>
//#include <Rmath.h>
//#include <R_ext/Applic.h>
//#include <R_ext/Error.h>
//#include <R_ext/Memory.h>
//#include <R_ext/Print.h>
//#include <R_ext/Blas.h>
//#include <R_ext/Linpack.h>
//#include <R_ext/RS.h> /* for F77_CALL */
//#include <R_ext/Linpack.h>
//#include <R_ext/Applic.h>
//#include <R_ext/Print.h> /* Rprintf */







#include <stdlib.h>
#include <stdio.h>

#include "lbfgsb.c"

static double * vect(int n)
{
    return (double *)R_alloc(n, sizeof(double));
}


double branin(int n, double xx[2], void * args) {

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

/* include setulb() */

//#include "lbfgsb.c"
//#include <R_ext/Error.h>



typedef double (*optimfnC)(int, double *, void *);
typedef void (*optimgrC)(int, double *, double *, void *);

typedef struct opt_struct
{
    optimfn *fcall;    /* function */
    optimgr *gcall;    /* gradient */
    double* ndeps;   /* tolerances for numerical derivatives */
    double fnscale;  /* scaling for objective */
    double* parscale;/* scaling for parameters */
    int usebounds;
    double* lower, *upper;
} opt_struct, *OptStruct;



void lbfgsb(int n, int m, double *x, double *l, double *u, int *nbd,
        double *Fmin, optimfn fminfn, optimgr fmingr, int *fail,
        void *ex, double factr, double pgtol,
        int *fncount, int *grcount, int maxit, char *msg,
        int trace, int nREPORT)
{
    char task[60];
    double f, *g,  *wa;
    int tr = -1, iter = 0, *iwa, isave[21];
    isave[12] = 0; // -Wall

    if(n == 0) { /* not handled in setulb */
    *fncount = 1;
    *grcount = 0;
    *Fmin = fminfn(n, u, ex);
    strcpy(msg, "NOTHING TO DO");
    *fail = 0;
    return;
    }
    if (nREPORT <= 0)
    error(_("REPORT must be > 0 (method = \"L-BFGS-B\")"));
    switch(trace) {
    case 2: tr = 0; break;
    case 3: tr = nREPORT; break;
    case 4: tr = 99; break;
    case 5: tr = 100; break;
    case 6: tr = 101; break;
    default: tr = -1; break;
    }

    *fail = 0;
    g = vect(n);
    /* this needs to be zeroed for snd in mainlb to be zeroed */
    wa = (double *) S_alloc(2*m*n+4*n+11*m*m+8*m, sizeof(double));
    iwa = (int *) R_alloc(3*n, sizeof(int));
    strcpy(task, "START");
    while(1) {
    setulb(n, m, x, l, u, nbd, &f, g, factr, &pgtol, wa, iwa, task,
           tr, isave);
/*  Rprintf("in lbfgsb - %s\n", task);*/
    if (strncmp(task, "FG", 2) == 0) {
        f = fminfn(n, x, ex);
        if (!R_FINITE(f))
        error(_("L-BFGS-B needs finite values of 'fn'"));
        fmingr(n, x, g, ex);
    } else if (strncmp(task, "NEW_X", 5) == 0) {
        iter++;
        if(trace == 1 && (iter % nREPORT == 0)) {
        Rprintf("iter %4d value %f\n", iter, f);
        }
        if (iter > maxit) {
        *fail = 1;
        break;
        }
    } else if (strncmp(task, "WARN", 4) == 0) {
        *fail = 51;
        break;
    } else if (strncmp(task, "CONV", 4) == 0) {
        break;
    } else if (strncmp(task, "ERROR", 5) == 0) {
        *fail = 52;
        break;
    } else { /* some other condition that is not supposed to happen */
        *fail = 52;
        break;
    }
    }
    *Fmin = f;
    *fncount = *grcount = isave[12];
    if (trace) {
    Rprintf("final  value %f \n", *Fmin);
    if (iter < maxit && *fail == 0) Rprintf("converged\n");
    else Rprintf("stopped after %i iterations\n", iter);
    }
    strcpy(msg, task);
}














int main(){
    int npar = 2,lmm = 5;
    double dpar[2] = {1., 1};
    double lower[2] = {-5, 0};
    double upper[2] = {10, 15};
    int nbd[2]= {2, 2};
    double Fmin=1000;
    //  printf("%f", ML_NEGINF);

    int ifail = 0;

    int trace = 0, maxit=100, fncount = 0,
        grcount = 0, nREPORT=1;
    double factr= 1e+07, pgtol = 0;
    //const char *tn;
    char msg[60];

    OptStruct OS = (OptStruct)malloc(sizeof(opt_struct));;
    OS->fcall = branin;
    OS->gcall = braningr;
    OS->parscale = (double*)malloc(npar * sizeof(double));
    OS->ndeps = (double*)malloc(npar * sizeof(double));
    OS->lower = (double*)malloc(npar * sizeof(double));
    OS->upper = (double*)malloc(npar * sizeof(double));
    OS->usebounds = 2;

   for(int i = 0; i<npar; i++) {
      OS->ndeps[i] = dpar[i];
     OS->parscale[i] = 0.01;
     OS->lower[i] = lower[i];
     OS->upper[i] = upper[i];
   };
   // for(int i = 0; i<npar; i++) printf("%f ", OS->parscale[i]);

//   void lbfgsb(int n, int m, double *x, double *l, double *u, int *nbd,
//       double *Fmin, optimfn fminfn, optimgr fmingr, int *fail,
//       void *ex, double factr, double pgtol,
//       int *fncount, int *grcount, int maxit, char *msg,
//       int trace, int nREPORT)
// */


   //printf("%f\n", Fmin);

   //lbfgsb(npar, lmm, dpar, lower, upper, nbd, &Fmin, &branin, braningr,
    //     &ifail, (void *)OS, factr, pgtol, &fncount, &grcount,
     //    maxit, msg, trace, nREPORT);

    double g[2] = {0,0};
    double *wa = (double *) calloc(2*lmm*npar+4*npar+
            11*lmm*lmm+8*lmm, sizeof(double));

    int *iwa = (int *) calloc(3*npar, sizeof(int));
    int tr = -1, isave[21];
    isave[12] = 0;
    char task[60];

    setulb(npar, lmm, dpar, lower, upper,
        nbd, &Fmin, g,
        factr, &pgtol, wa, iwa, task, tr, isave);

   for(int i = 0; i<npar; i++)
        printf("%f ", dpar[i]);
    printf("\n");
   for(int i = 0; i<5; i++) printf("%c", task[i]);
        printf("\n");


     // Cleanup
//   free(OS->parscale);
//   free(OS->ndeps);
//   free(OS->lower);
//   free(OS->upper);
//   free(OS);
//  printf("%f", Fmin);

}

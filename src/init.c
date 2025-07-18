// Required R and C headers
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <time.h>
#include <omp.h>  // OpenMP for parallelism

// Structure to hold the minimum value and corresponding index
typedef struct {
  double value;
  int index;
} MinPair;


// Structure for Differential Evolution working memory
typedef struct {
  double *phiValues;   // Objective function values
  int *agent, *potential; // Current and proposed population
  MinPair globalMin;   // Best (minimum) result found
} *dePtr;

// Parameters structure for the DE algorithm
typedef struct {
  int n, m, s, size, NP, itermax;
  int seed, trace, cores, len, r, estimate;
  double pMut, pCR, pGBest, C, denom_g;
} params, *paramsPtr;


// Pointer to a function taking int* and paramsPtr, returning double
typedef double (*criteria)(int *, paramsPtr);


// Utility to allocate an int vector of given size
inline static int * allocVec(int size){
  int *vec = (int*)malloc(size * sizeof(int));
  return vec;
}


// MaxPro criterion function
// https://www.asc.ohio-state.edu/statistics/comp_exp/jour.club/max_projection_designs_roshan2015.pdf
// Equation 5
// \begin{equation}
//  \min _D \psi(D)=\left\{\frac{1}{\binom{n}{2}}
//  \sum_{i=1}^{n-1} \sum_{j=i+1}^n \frac{1}{\prod_{l=1}^m
//  \left(x_{i l}-x_{j l}\right)^2}\right\}^{1 / m}
//\end{equation}


double maxpro(int * data, paramsPtr p){
  // Penalize small distances between points
  double dsq = 0;
  int len=0;
  for(int i = 0; i < p->n - 1; ++i){
    for(int j = i + 1; j < p->n; ++j, ++len){
      double row = 1;
      for(int k = 0; k < p->m; ++k){
        row *= data[i + k * p->n] - data[j + k * p->n];
        if (row == 0) return INFINITY; // Degenerate design
      }
      dsq += 1/(row * row);
    }
  }
  return pow(dsq/len, 1./p->m);

}

// UniPro criterion function
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


// Maximin LHD criterion
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

// Available objective functions
static criteria PHI_FUNS[3] = {unipro, maximinLHD, maxpro};


// Uniform random integer and double utilities
#define intUniform(n) (n < 1 ? 0 : rand() % (n))
#define doubleUniform ((double)rand() / RAND_MAX)

// Select two distinct integers between 0 and n-1
int *sampleTwo(int n){
  int * vec =allocVec(n);
  vec[0] = intUniform(n);
  do vec[1] = intUniform(n); while(vec[0] == vec[1]);
  return vec;
}

// Swap elements of a vector
void swapVecAtPositions(int * vec, int i, int j){
  int temp = vec[i];
  vec[i] = vec[j];
  vec[j] = temp;
}
// Swap based on another vector (of 2 positions)
void swapVecUsingVec(int *vec, int *v){
  swapVecAtPositions(vec, v[0], v[1]);
}


// Generate random Latin Hypercube design
void randomlhs(dePtr de, paramsPtr p, criteria phi){
  int *index = allocVec(p->n);
#pragma omp parallel for simd simdlen(8)
  for(int i = 0; i < p->n; ++i) index[i] = i % p->s;
#pragma omp parallel for
  for (int j = 0; j < p->m*p->NP; ++j)
  {
    for (int i = p->n - 1; i > 0; --i)
    {
      swapVecAtPositions(index, i, intUniform(i));
      de->agent[i + j * p->n] = index[i];
    }
    de->agent[j * p->n] = index[0];
    if(j % p->m == p->m - 1){
      int k = j / p->m;
      double value = phi(de->agent + p->size * k, p);
      de->phiValues[k] = value;
      if(value < de->globalMin.value){
        de->globalMin.value = value;
        de->globalMin.index = k;
      }
    }
  }
  free(index);
}

// Initialize parameter struct
paramsPtr initParams(int n, int m, int s, int NP, int itermax,
                     double *pMut,double *pCR, double *pGBest,
                     int seed, int trace, int cores, int r){
  // Initialize and compute constants needed for evaluation
  paramsPtr p = (paramsPtr)malloc(sizeof(params));
  p->n = n;
  p->m = m;
  p->s = s;
  p->NP = NP;
  p->itermax = itermax;
  p->seed = seed;
  p->pMut = pMut == NULL? 0.2: *pMut;
  p->pCR = pCR == NULL? 0.3: *pCR;
  p->pGBest = pGBest == NULL? 0.95: *pGBest;
  p->estimate = pMut == NULL||pCR == NULL||pGBest == NULL;
  int numThreads = omp_get_max_threads();
  if(trace >1) printf("\n Total number of Cores: %d ---- Using: ", numThreads);
  p->cores = numThreads > cores? cores: numThreads - 1;
  if(trace>1) printf("%d\n", p->cores);
  p->trace = trace;
  p->size = n*m;
  p->len = p->size * NP;
  p->r = r;
  double p4 = pow(n, 4);
  double numerator = 4*(5.0*m - 2)*p4 + 30.0*(3*m-5)*s*s + 15*m + 33;
  double denom = 720.0*(m - 1)*p4;
  p->C =  numerator/ denom + (1.0 + pow(-1,s))/(64.0*p4);
  p->denom_g =   4.0*m * (m - 1) * n*n*s*s;
  return p;
}



// Allocate memory and initialize DE design
dePtr initialize(const paramsPtr p, criteria phi){
  dePtr de;
  de = (dePtr) malloc(sizeof(*de));
  de->agent = allocVec(p->len);
  de->potential = allocVec(p->len);
  de->phiValues = (double*)malloc(p->NP * sizeof(double));
  de->globalMin.value = 99999;
  randomlhs(de, p, phi);
  return de;
}


// Propose a new design from global best or current
void propose(dePtr de, int NP_INDEX, paramsPtr p){
  // with a probability pGBest, select between
  // the global best design and the current design
  // Other proposal methods could be implemented
  memcpy(de->potential + NP_INDEX * p->size,
         de->agent + (doubleUniform < p->pGBest?
                        de->globalMin.index : NP_INDEX) * p->size,
                        p->size * sizeof(int));

}

// Perform mutation on candidate design
void mutate(dePtr de, int offset, paramsPtr p){
  // With probability pMut, consider mutating a column.
  // mutation involves the swap operator
  // where two randomly selected values in the column are swaped
  if(doubleUniform < p->pMut)
    swapVecUsingVec(de->potential + offset, sampleTwo(p->n));
}

// Perform crossover from current to proposed
void crossOver(dePtr de, int offset, int changeJ, paramsPtr p){
  // Using the proposed design,
  // With a probability pCR choose between the mutated column
  // and the original column. Note that one column in the proposed
  // must be maintained.

  if(doubleUniform > p->pCR && changeJ)
    memcpy(de->potential+offset,
           de->agent + offset, p->n*sizeof(int));
}


// Full trial generation and selection step
void getTrial(dePtr de, int NP_INDEX, paramsPtr p, criteria phi){
  // Take an agent design, propose, mutate and do crossover.
  // This gives us the trial design. We then select between
  // the trial design and the original design based n the phi value

  propose(de, NP_INDEX, p);
  int start =  NP_INDEX * p->size;
  int maintain = intUniform(p->m);
  for(int j = 0; j < p->m; ++j){
    int offset = start + j * p->n;
    mutate(de, offset, p);
    crossOver(de, offset, maintain!=j, p);
  }
  double value = phi(de->potential + start, p);
  if(value < de->phiValues[NP_INDEX]){
    de->phiValues[NP_INDEX] = value;
    memcpy(de->agent + start, de->potential + start,
           p->size * sizeof(int));
  }
  if(value < de->globalMin.value){
    de->globalMin.value = value;
    de->globalMin.index = NP_INDEX;
  }
}

// Print and update a text-based progress bar
void print_progress(int completed, int total) {
  int bar_width = 50;
  float progress = (float)completed / total;

  Rprintf("\r[");
  int pos = bar_width * progress;
  for (int i = 0; i < bar_width; ++i) {
    if (i < pos) Rprintf("=");
    else if (i == pos) Rprintf(">");
    else printf(" ");
  }
  Rprintf("] %d%%", (int)(progress * 100));
  fflush(stdout);
}

void progressbar(int * completed, int replications){
  //Update progress atomically
#pragma omp atomic
  (*completed)++;
  //Print progress when enough progress has been made
  if (*completed % (replications / 100 + 1) == 0 ||
      *completed == replications)
  {
#pragma omp critical
    print_progress(*completed, replications);
  }
}

// Run DE for a single replicate
double DE1(paramsPtr p, int *bestX, criteria phi) {
  double value = 0;

#pragma omp parallel for num_threads(p->cores)
  for (int reps = 0; reps < 1; reps++) {
    dePtr de = initialize(p, phi);
#pragma omp parallel for
    for (int it = 0; it < p->itermax; it++) {
      for (int i = 0; i < p->NP; ++i) {
        getTrial(de, i, p, phi);
      }
    }
    memcpy(bestX,
           de->agent + p->size * de->globalMin.index,
           p->size * sizeof(int));

    value = de->globalMin.value;
  }

  return value;
}

// Estimate the best (pMut, pCR) for this problem using grid search
void estimate_params(paramsPtr shared_p, criteria phi) {
  double min_val = 999;
  double best_pmut = 0, best_pcr = 0;
#pragma omp parallel for reduction(min:min_val)
  for (int i = 0; i < 5; i++) {
    params p = *shared_p;
    p.itermax = 100;
    p.trace = 0;
    int *localX = malloc(p.n * p.m * sizeof(int));
    p.pMut = 0.1 + 0.2 * i;
    p.pCR  = 0.9 - 0.2 * i;
    double val = DE1(&p, localX, phi);

#pragma omp critical
{
  if (val < min_val) {
    min_val = val;
    best_pmut = p.pMut;
    best_pcr = p.pCR;
  }
}
free(localX);
  }

  shared_p->pMut = best_pmut;
  shared_p->pCR = best_pcr;
}

// Core DE loop with OpenMP parallelization across replicates
void DE_C(const paramsPtr p, criteria phi, int replications, double * vals, int * bestX){
  double min_val;
  int completed = 0;
#pragma omp parallel for num_threads(p->cores)
  for(int reps = 0; reps < replications; reps++){
    int *localX = allocVec(p->n*p->m);
    params local_p = *p;
    srand(p->seed + reps);
    if(p->estimate) estimate_params(&local_p, phi);
    double value = DE1(&local_p,  localX, phi);
    vals[reps] = value;

    if(value < min_val){
      min_val = value;
      memcpy(bestX, localX, p->size*sizeof(int));
    }
    if(p->trace > 1) progressbar(&completed, replications);
    free(localX);
  }
}

// R interface: detect number of available CPU cores
SEXP detectCores(){
  SEXP R_out = PROTECT(allocVector(INTSXP, 1));
  INTEGER(R_out)[0] =  omp_get_max_threads();
  UNPROTECT(1);
  return R_out;
}

// R interface for evaluating a criterion function
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
       initParams(n, m, s, 0, 0, 0,
                  0, 0, 1, 0, 2, r));
  UNPROTECT(1);
  return R_out;
}



// R interface to DE() optimization
SEXP DE(SEXP R_n, SEXP R_m, SEXP R_s, SEXP R_NP,
        SEXP R_itermax, SEXP R_pMut, SEXP R_pCR,
        SEXP R_pGBest, SEXP R_replicates,
        SEXP R_cores, SEXP R_method, SEXP R_r,
        SEXP R_trace, SEXP R_seed, SEXP R_shared){

  struct timespec start, end;
  clock_gettime(CLOCK_MONOTONIC, &start);

  // Unpack all inputs, initialize, run DE_C, wrap output
  int n = INTEGER(R_n)[0];
  int m = INTEGER(R_m)[0];
  int s = INTEGER(R_s)[0];
  int r = INTEGER(R_r)[0];
  int NP = INTEGER(R_NP)[0];
  int itermax = INTEGER(R_itermax)[0];
  int replicates = INTEGER(R_replicates)[0];
  int cores = INTEGER(R_cores)[0];
  int trace = INTEGER(R_trace)[0];
  int seed = INTEGER(R_seed)[0];

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
  if(trace)
    switch (method) {
    case 0:
      Rprintf("UniPro\n");
      break;
    case 2:
      Rprintf("MaxPro\n");
      break;
    case 1:
      Rprintf("maximinLHD\n");
      break;
    default:
      Rprintf("UniPro\n");
    break;
    }
  criteria phi = PHI_FUNS[method];


  double mut, cr, gbest;
  double *pMut = NULL, *pCR = NULL, *pGBest = NULL;

  if (!Rf_isNull(R_pMut)) {
    mut = REAL(R_pMut)[0];
    pMut = &mut;
  }

  if (!Rf_isNull(R_pCR)) {
    cr = REAL(R_pCR)[0];
    pCR = &cr;
  }

  if (!Rf_isNull(R_pGBest)) {
    gbest = REAL(R_pGBest)[0];
    pGBest = &gbest;
  }


  paramsPtr p = initParams(n, m, s, NP, itermax, pMut,
                           pCR, pGBest, seed, trace, cores, r);

  DE_C(p, phi, replicates, phiValues, bestX);

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
  clock_gettime(CLOCK_MONOTONIC, &end);
  *timeTaken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
  if(trace > 1) Rprintf("  %8.3f Secs\n", *timeTaken);
  UNPROTECT(7);
  return R_out;
}



// Register the function
static const R_CallMethodDef CallMethods[] = {
  {"phi", (DL_FUNC) &phi, 3},
  {"DE", (DL_FUNC) &DE, 15},
  {"detectCores", (DL_FUNC) &detectCores, 0},
  {NULL, NULL, 0}
};

// Called automatically when R loads the package
void R_init_RcppProgressExample(DllInfo *Info) {
  R_registerRoutines(Info, NULL, CallMethods, NULL, NULL);
  R_useDynamicSymbols(Info, FALSE);
}




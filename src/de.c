#include "criteria.h"
#include "de.h"



#define intUniform(n) (n < 1 ? 0: rand() % (n))
#define doubleUniform ((double)rand() / RAND_MAX)

inline static int * allocVec(int size){
  int *vec = (int*)malloc(size * sizeof(int));
  return vec;
}

int *sampleTwo(int n){
  int * vec = (int*)malloc(n * sizeof(int));
  vec[0] = intUniform(n);
  do vec[1] = intUniform(n); while(vec[0] == vec[1]);
  return vec;
}


void swapVecAtPositions(int * vec, int i, int j){
  int temp = vec[i];
  vec[i] = vec[j];
  vec[j] = temp;
}

void swapVecUsingVec(int *vec, int *v){
  swapVecAtPositions(vec, v[0], v[1]);
}



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


paramsPtr initializeParams(int n, int m, int s, int NP, int itermax,
                           double pMut,double pCR, double pGbest, int replications,
                           long int seed, int r){
  paramsPtr p = (paramsPtr)malloc(sizeof(params));
  p->n = n;
  p->m = m;
  p->s = s;
  p->NP = NP;
  p->itermax = itermax;
  p->pMut = pMut;
  p->pCR = pCR;
  p->pGbest = pGbest;
  p->replications = replications;
  p->seed = seed;
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



void propose(dePtr de, int NP_INDEX, paramsPtr p){
  // with a probability pGbest, select between
  // the global best design and the current design
  // Other proposal methods could be implemented
  memcpy(de->potential + NP_INDEX * p->size,
         de->agent + (doubleUniform < p->pGbest?
                        de->globalMin.index : NP_INDEX) * p->size,
                        p->size * sizeof(int));

}

void mutate(dePtr de, int offset, paramsPtr p){
  // With probability pMut, consider mutating a column.
  // mutatin involves the swap operator
  // where two randomly selected values in the column are swaped
  if(doubleUniform < p->pMut)
    swapVecUsingVec(de->potential + offset, sampleTwo(p->n));
}

void crossOver(dePtr de, int offset, int changeJ, paramsPtr p){
  // Using the proposed design,
  // With a probability pCR choose between the mutated column
  // and the original column. Note that ne column in the proposed
  // must be maintained.

  if(doubleUniform < p->pCR && changeJ)
    memcpy(de->potential+offset,
           de->agent + offset, p->n*sizeof(int));
}

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

void print_progress(int completed, int total) {
    int bar_width = 50;
    float progress = (float)completed / total;

    printf("\r[");
    int pos = bar_width * progress;
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) printf("=");
        else if (i == pos) printf(">");
        else printf(" ");
    }
    printf("] %d%%", (int)(progress * 100));
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



void DE_CC(int n, int m, int s, int NP, int itermax,
          double pMut, double pCR, double pGbest,
          int replications, unsigned int seed, double * vals,
          double *timeTaken, int * bestX, int numCores, criteria phi, int r)
{
  struct timespec start, end;
  clock_gettime(CLOCK_MONOTONIC, &start);
  srand(seed);
  double Min = 999.;
  int numThreads = omp_get_max_threads();
  printf("\n Total number of Cores: %d ---- Using: ", numThreads);
  paramsPtr p = initializeParams(n, m, s, NP, itermax,
                                 pMut, pCR, pGbest, replications, seed, r);
  int completed = 0;
  numThreads = numThreads > numCores? numCores: numThreads - 1;
  printf("%d\n", numThreads);
  #pragma omp parallel for num_threads(numThreads)
  for(int reps = 0; reps < replications; reps++){
    dePtr de = initialize(p, phi);
    #pragma omp parallel for
    for(int it = 0; it < itermax; it++)
      for(int i = 0; i<NP; ++i) getTrial(de, i, p, phi);
    double value =  de->globalMin.value;
    vals[reps] = value;
    if(value < Min){
      Min = value;
      memcpy(bestX, de->agent + p->size * de->globalMin.index,
             p->size*sizeof(int));
    }
    progressbar(&completed, replications);
  }
  clock_gettime(CLOCK_MONOTONIC, &end);
  *timeTaken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
  printf("  %8.3f Secs\n", *timeTaken);
}




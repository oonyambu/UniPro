#include "de.h"

void sample2(int n, int vec[2]){
  vec[0] = rand_int(n);
  do vec[1] = rand_int(n); while(vec[1] == vec[0]);
}

void print_design(int *X, int n, int m, int NP){
  for(int k=0; k<NP; k++){
    for(int i = 0; i<n; i++){
      for(int j = 0; j<m; j++)
        printf("%d ", X[i + n*(j + k*m)]);
      printf("\n");
    }
    printf("\n\n");
  }
}



double phi2D(const int *X, const int n, const int m, const int s, double denom_g, double C){
  double dist = 0;
  double dist_row = 0;
  double v[1000] = {0};
  //#pragma omp parallel for reduction(+:dist_row) private(j)
  for (int i = 0; i < n; i++){
    for (int j = i+1; j < n; j++){
      double d = 0;
      for(int k = 0; k<m; k++)
        d += abs(X[i + k*n] - X[j + k*n]);
      dist += 2*d*d;
      v[j] += d;
      v[i] += d;
    }
    dist_row += v[i]*v[i];
  }
  return ((dist - 2.0*dist_row/n)/denom_g + C)*1000;
}




void gen_design(int *array, const int n, const int m,
                const int s, const int NP, double *phi_vals,
                int *bestIndex, double denom_g, double C){

  double val = 9999999.;
  int k;
   MinInfo global_min = {val, -1};
  for(k = 0; k<NP; k++){
    for(int i = 0; i<n; i++){
      for(int j=0; j<m; j++)
        array[i + n*(j + k*m)] = i % n + 1;
    }
    for(int i = 0; i < n-1; i++)
      for (int j = 0; j < m; j++)
        swap(&array[i + n*(j + m*k)], &array[rand_int(n) + n*(j + m*k)]);
    phi_vals[k] = phi2D(array + k*n*m, n, m, s, denom_g, C);

    if (phi_vals[k] < val) {
      global_min.value = phi_vals[k];
      global_min.index = k;
    }
  }
  *bestIndex = global_min.index;
}

void phi_C(const int *X, const int *n, const int *m, const int *s, double*val){
  double p = pow(*s, 4);
  double numerator = 4*(5.0**m - 2)*p + 30.0*(3**m-5)**s**s + 15**m + 33;
  double denom = 720.0*(*m - 1)*p;
  double C = numerator/ denom + (1.0 + pow(-1,*s))/(64.0*p);
  double denom_g = 4.0**m * (*m - 1) * *n**n**s**s;
  *val = phi2D(X, *n, *m, *s, denom_g, C);
}


void DE_C(int *n, int *m, int *s, int *NP, int *itermax,
        double *pMut, double *pCR, double *pGBest, int *X,
        double * opt_val, double *time_taken, int *replicates){

  struct timespec start, end;
  clock_gettime(CLOCK_MONOTONIC, &start);
  double p = pow(*s, 4);
  double numerator = 4*(5.0**m - 2)*p + 30.0*(3**m-5)**s**s + 15**m + 33;
  double denom = 720.0*(*m - 1)*p;
  double C = numerator/ denom + (1.0 + pow(-1,*s))/(64.0*p);
  double denom_g = 4.0**m * (*m - 1) * *n**n**s**s;

  int* array = (int*)malloc(*n**m**NP * sizeof(int));
  int* potential = (int*)malloc(*n**m**NP * sizeof(int));
  double *vals = (double*)malloc(*NP *sizeof(double));
  int reps, final;
  //#pragma omp parallel for private(reps)
  for(reps = 0; reps < *replicates; reps++){
    int bestIndex;
    gen_design(array, *n, *m, *s, *NP, vals, &bestIndex, denom_g, C);
    for(int it = 0; it< *itermax; it++){
      for(int k = 0; k< *NP; k++){
        int gl = rand_unif() < *pGBest;
        for(int j = 0; j < *m; j++){
          for (int i = 0; i < *n; i++)
            potential[i + *n*(j + *m*k)] = gl? array[i + *n*(j + *m*bestIndex)] : array[i + *n*(j + *m*k)];
          if (rand_unif() < *pMut){
            int v[2];
            sample2(*n, v);
            swap(&potential[v[0] + *n*(j + *m*k)], &potential[v[1] + *n*(j + *m*k)]);
          }
          int CR = rand_int(*m);
          double pcr = rand_unif();
          for(int i = 0; i<*n; i++)
            if( (pcr > *pCR) & (j != CR))
              potential[i + *n*(j + *m*k)] = array[i + *n*(j + *m*k)];
        }
        double phi_val = phi2D(potential + k**n**m, *n, *m, *s, denom_g, C);
        if(phi_val < vals[k]) {
          vals[k] = phi_val;
          for (int i=0; i<*n; i++)
            for(int j=0; j<*m; j++) array[i + *n*(j + *m*k)] = potential[i + *n*(j + *m*k)];
        }
        if(phi_val < vals[bestIndex]) bestIndex = k;
      }
    }
    opt_val[reps] = vals[bestIndex];
  }

  // if(*replicates == 1){
  //   for(int i=0;i<*n;i++)
  //     for(int j=0; j<*m; j++)
  //       X[i + *n*j] = array[i + *n*(j + *m*bestIndex)];
  // }

  free(array);
  free(potential);
  free(vals);

  clock_gettime(CLOCK_MONOTONIC, &end);
  *time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
}


int main(){
  int n = 30, m = 3, s = 30, itermax = 1000, NP=100, replicates=10, X[90];
  double pMut = 0.2, pCR=0.3, pGBest = 0.9, opt_val[10], time_taken;

  DE_C(&n, &m, &s, &NP, &itermax,&pMut, &pCR, &pGBest, X,
       opt_val, &time_taken, &replicates);
  //DE_C(n, m, s, NP, itermax,pMut, pCR, pGBest, X,
   //     opt_val, &time_taken, replicates);
  for(int i=0; i<replicates;i++)printf("%f\n", opt_val[i]);
}

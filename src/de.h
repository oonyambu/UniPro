#ifndef __UNIPRO__
#define __UNIPRO__

#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <string.h>


typedef struct {
  double value;
  int index;
} MinPair;


typedef struct {
  double *phiValues;
  int *agent, *potential;
  MinPair globalMin;
} *dePtr;


typedef struct {
  int n, m, s, size, NP,  itermax;
  int replications, seed, len, r;
  double pMut, pCR, pGbest, C, denom_g;
} params,*paramsPtr;


paramsPtr initializeParams(int n, int m, int s, int NP, int itermax,
                           double pMut,double pCR, double pGbest, int replications,
                           long int seed, int r);

typedef double (*criteria)(int *, paramsPtr);

void DE_CC(int n, int m, int s, int NP, int itermax,
           double pMut, double pCR, double pGbest,
           int replications, unsigned int seed, double * vals,
           double *timeTaken, int * bestX,  int numCores, criteria phi, int r);
#endif

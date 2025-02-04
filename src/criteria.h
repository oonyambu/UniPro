#ifndef __CRITERIA__
#define __CRITERIA__
#pragma once
#include "de.h"


double maxpro(int * data, paramsPtr p);
double unipro(int * data, paramsPtr p);
double maximinLHD(int * data, paramsPtr p); //  mean(1/dist(D)^15)^(1/15)

#endif //__CRITERIA__

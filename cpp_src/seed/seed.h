#ifndef SIM_2D_CPP_SEED_H
#define SIM_2D_CPP_SEED_H

#include "../config.h"

long double unif_rand();

void set_seed(unsigned int seed);

void runif_seq(long double *seq, int len, long double min, long double max);

int unif_index(int dn);

#endif //SIM_2D_CPP_SEED_H

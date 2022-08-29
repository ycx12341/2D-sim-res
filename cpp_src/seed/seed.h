#ifndef SIM_2D_CPP_SEED_H
#define SIM_2D_CPP_SEED_H

#include "../config.h"

#include <array>

DBL_T unif_rand();

void set_seed(unsigned int seed);

void runif_seq(DBL_T *seq, int len, DBL_T min, DBL_T max);

int unif_index(int dn);

std::array<int, 2> unif_index2(int dn);

int sample_prob1(int dn, const DBL_T *prob);

#endif //SIM_2D_CPP_SEED_H

#ifndef SIM_2D_CPP_SEED_H
#define SIM_2D_CPP_SEED_H

#include "../config.h"

#include <vector>

DBL_T unif_rand();

void set_seed(unsigned int seed);

void runif_seq(DBL_T *seq, int len, DBL_T min, DBL_T max);

int unif_index(int dn);

std::vector<int> unif_index(int index_num, int dn);

int sample_int_index(int dn, const DBL_T *prob);

std::vector<int> sample_indices(int sample_num, const std::vector<DBL_T> &prob, bool replace);

/**
 * Default N01_kind: INVERSION
 * @return
 */
DBL_T norm_rand();

#endif //SIM_2D_CPP_SEED_H

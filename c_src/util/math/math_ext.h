/**
 * Math extensions.
 * Derived from https://svn.r-project.org/R/
 */

#ifndef C_SRC_MATH_EXT_H
#define C_SRC_MATH_EXT_H

#include <stdbool.h>

/**
 * Set the initial seed for variate generators.
 */
void set_seed(unsigned int seed);

/**
 * Provides information about the uniform distribution on the interval from min to max.
 * @param min lower and upper limits of the distribution. Must be finite.
 * @param max upper limits of the distribution. Must be finite.
 * @return runif(1, min, max)[1]
 */
double runif(double a, double b);

void runif_seq(double *seq, int len, double min, double max);

double unif_rand();

int unif_index(int dn);

typedef struct {
    int arr[2];
} int_arr_2_t;

int_arr_2_t unif_index2(int dn);

int sample_prob1(int dn, const double prob[dn]);

#endif //C_SRC_MATH_EXT_H

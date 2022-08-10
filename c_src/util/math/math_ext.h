/**
 * Math extensions.
 * Derived from https://svn.r-project.org/R/
 */

#ifndef C_SRC_MATH_EXT_H
#define C_SRC_MATH_EXT_H

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

double unif_rand();

/**
 * Generate regular sequences.
 * @param buf   A buffer to store the generated sequence.
 * @param from  Lower bound of the sequence.
 * @param to    Upper bound of the sequence.
 * @param num   Length of sequence.
 * @return The number of sequence.
 */
int seq(double *buf, double from, double to, int num);

#endif //C_SRC_MATH_EXT_H

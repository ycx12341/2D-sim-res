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

#endif //C_SRC_MATH_EXT_H

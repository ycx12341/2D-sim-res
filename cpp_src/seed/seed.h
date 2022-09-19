/**
 * Random Samples
 * Permutations
 * Uniform Distribution
 * Randomness
 * Seeding
 *
 * Derived From https://svn.r-project.org/R/ branch R-4-2-branch.
 */

#ifndef SIM_2D_CPP_SEED_H
#define SIM_2D_CPP_SEED_H

#include "../config.h"

#include <vector>

/**
 * Continuous uniform random numbers.
 * @return DBL_T random number.
 */
DBL_T unif_rand();

/**
 * Seeding Random Variate Generators.
 * Using "Mersenne-Twister" approach.
 * @param seed Seed.
 *             When set to 0, use system time for seed.
 */
void set_seed(unsigned int seed);

/**
 * Random Uniform Distribution.
 * Equivalent to:
 *      seq = runif(n = len, min = min, max = max)
 * @param len Number of Observations.
 * @param seq DBL_T array to hold the random deviates.
 * @param min Lower limit of the distribution.
 * @param max Upper limit of the distribution.
 */
void runif_seq(unsigned len, DBL_T *seq, DBL_T min, DBL_T max);

/**
 * Takes a sample from [0, dn-1].
 * Equivalent to:
 *      sample(x = [0:dn-1], size = 1, replace = FALSE, prob = NULL)
 * @param dn The upper limit of sample area.
 * @return An integer in range of [0, dn-1].
 */
unsigned int unif_index(unsigned int dn);

/**
 * Takes multiple samples from [0, dn-1]
 * Equivalent to:
 *      sample(x = [0:dn-1], size = index_num, replace = FALSE, prob = NULL)
 * @param index_num The number of samples.
 * @param dn        The upper limit of sample area.
 * @return A vector containing samples with size `index_num`.
 */
std::vector<unsigned int> unif_index(unsigned int index_num, unsigned int dn);

/**
 * Takes a sample from [1, dn] based on their probability weights.
 * Equivalent to:
 *      sample.int(n = dn, size = 1, replace = FALSE, prob = prob)
 * @param dn   The upper limit of sample area.
 * @param prob Array indicating probability weights for obtaining the elements of the vector being sampled.
 *             The size of the array must be equal to dn.
 * @return An integer in range [1, dn]
 */
unsigned int sample_int_index(unsigned int dn, const DBL_T *prob);

/**
 * Takes samples from [1, prob.size()-1] based on their probability weights.
 * Equivalent to:
 *      sample(x = [0, prob.size()-1], size = sample_num, replace = replace, prob = prob)
 * @param sample_num The number of samples.
 * @param prob       Array indicating probability weights for obtaining the elements of the vector being sampled.
 *                   (The size of the array - 1) will be the upperbound of the sample area.
 * @param replace    Should sampling be with replacement?
 *                   If false, sample_num must be equal to or less than prob.size().
 * @return A vector containing samples with size `sample_num`.
 */
std::vector<unsigned int>
sample_indices(unsigned int sample_num, const std::vector<DBL_T> &prob, bool replace);

/**
 * Random Normal Distribution.
 * Equivalent to:
 *      rnorm(1, mean = mu, sd = sigma)
 * @param mu    Mean.
 * @param sigma Standard Deviation.
 * @return observation.
 */
DBL_T rnorm(DBL_T mu, DBL_T sigma);

#endif //SIM_2D_CPP_SEED_H

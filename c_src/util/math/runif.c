/**
 * runif implementation.
 * Derived from https://svn.r-project.org/R/
 */

#include <math.h>

#include "math_ext.h"

double runif(const double min, const double max) {
    if (!isfinite(min) || !isfinite(max) || max < min) {
        return NAN;
    }
    if (min == max) { return min; }

    double u;
    do {
        u = unif_rand();
    } while (u <= 0 || u >= 1);
    return min + (max - min) * u;
}

void runif_seq(double* seq, const int len, const double min, const double max) {
    for (int i = 0; i < len; i++) {
        seq[i] = runif(min, max);
    }
}

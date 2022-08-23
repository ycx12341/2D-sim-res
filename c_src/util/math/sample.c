/**
 * runif implementation.
 * Derived from https://svn.r-project.org/R/
 */

#include <math.h>
#include <stdint.h>

#include "math_ext.h"

static double rbits(const int bits) {
    uint64_t v = 0;

    for (int n = 0; n <= bits; n += 16) {
        int v1 = (int) floor(unif_rand() * 65536);
        v = 65536 * v + v1;
    }

    const int_least64_t one64 = 1L;
    return (double) (v & ((one64 << bits) - 1));
}

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

void runif_seq(double *seq, const int len, const double min, const double max) {
    for (int i = 0; i < len; i++) {
        seq[i] = runif(min, max);
    }
}

int unif_index(const int dn) {
    if (dn <= 0) { return 0; }
    int    bits = (int) ceil(log2(dn));
    double dv;
    do { dv = rbits(bits); } while (dn <= dv);
    return (int) dv;
}

int_arr_2_t unif_index2(const int dn) {
    int_arr_2_t ry;
    ry.arr[0] = -1;
    ry.arr[1] = -1;
    if (dn < 2) { return ry; }

    int x[dn], n = dn;

    for (int i = 0; i < dn; i++) { x[i] = i; }
    for (int i = 0; i < 2; i++) {
        int j     = unif_index(n);
        ry.arr[i] = x[j];
        x[j]      = x[--n];
    }
    return ry;
}

int sample_prob1(const int dn, const double prob[dn]) {
    int    x[dn];
    double prob_sort[dn];
    double rT, mass, totalmass;

//    for (int i = 0; i < dn; ++i) { x[i] = i + 1; }
    for (int i = 0; i < dn; ++i) { prob_sort[i] = prob[i]; }
    quickSort(dn, prob_sort, false);


    totalmass = 1;
    rT        = totalmass * unif_rand();
    mass      = 0;

    int i;
    for (i = 0; i < dn - 1; i++) {
        mass += prob_sort[i];
        if (rT <= mass)
            break;
    }
    return dn - i;
}
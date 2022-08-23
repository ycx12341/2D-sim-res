/**
 * runif implementation.
 * Derived from https://svn.r-project.org/R/
 */

#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <float.h>

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

void fixupProb(double *p, int n, int require_k) {
    double sum  = 0.0, pi;
    int    npos = 0;

    for (int i = 0; i < n; i++, npos++) {
        pi = p[i];
        assert(0 <= pi && pi < DBL_MAX);
        sum += p[i];
    }

    assert(npos != 0 && require_k <= npos);
    for (int i = 0; i < n; i++) { p[i] /= sum; }
}

/**
 * Sort a[] into descending order by "heapsort";
 * Sort ib[] alongside;
 * If initially, ib[] = 1...n, it will contain the permutation finally
 */
void revsort(double *a, int *ib, int n) {
    int    l, j, ir, i;
    double ra;
    int    ii;

    if (n <= 1) { return; }

    a--;
    ib--;

    l  = (n >> 1) + 1;
    ir = n;

    for (;;) {
        if (l > 1) {
            l  = l - 1;
            ra = a[l];
            ii = ib[l];
        } else {
            ra = a[ir];
            ii = ib[ir];
            a[ir]  = a[1];
            ib[ir] = ib[1];
            if (--ir == 1) {
                a[1]  = ra;
                ib[1] = ii;
                return;
            }
        }
        i = l;
        j = l << 1;
        while (j <= ir) {
            if (j < ir && a[j] > a[j + 1]) { ++j; }
            if (ra > a[j]) {
                a[i]  = a[j];
                ib[i] = ib[j];
                j += (i = j);
            } else {
                j = ir + 1;
            }
        }
        a[i]  = ra;
        ib[i] = ii;
    }
}

int sample_prob1(const int dn, const double prob[dn]) {
    int    perm[dn];
    double prob_cpy[dn];

    for (int i = 0; i < dn; i++) {
        prob_cpy[i] = prob[i];
        perm[i]     = i + 1;
    }

    fixupProb(prob_cpy, dn, false);
    revsort(prob_cpy, perm, dn);

    for (int i = 1; i < dn; i++) {
        prob_cpy[i] += prob_cpy[i - 1];
    }

    double rU = unif_rand();
    int    i;
    for (i = 0; i < dn - 1; i++) {
        if (rU <= prob_cpy[i]) { break; }
    }
    return perm[i];
}
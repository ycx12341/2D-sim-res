/**
 * Random Samples, Permutations and Uniform Distribution.
 *
 * Derived From https://svn.r-project.org/R/ branch R-4-2-branch.
 */

#include <cmath>
#include <cstdint>
#include <cassert>
#include <climits>

#include "seed.h"

static double rbits(const int bits) {
    int_least64_t v = 0;
    for (int      n = 0; n <= bits; n += 16) {
        int v1 = (int) floor(unif_rand() * 65536);
        v = 65536 * v + v1;
    }
    // mask out the bits in the result that are not needed
    return (double) (v & ((1L << bits) - 1));
}

DBL_T runif(const DBL_T min, const DBL_T max) {
    assert(std::isfinite(min) && std::isfinite(max) && min <= max);
    if (min == max) { return min; }

    DBL_T u;
    do {
        u = unif_rand();
    } while (u <= 0 || u >= 1);
    return min + (max - min) * u;
}

void runif_seq(const unsigned len, DBL_T *seq, const DBL_T min, const DBL_T max) {
    for (int i = 0; i < len; i++) {
        seq[i] = runif(min, max);
    }
}

unsigned int unif_index(const unsigned int dn) {
    if (dn <= 0) { return 0; }
    int    bits = (int) ceil(log2(dn));
    double dv;
    do { dv = rbits(bits); } while (dn <= dv);
    return (int) dv;
}

std::vector<int> unif_index(const int index_num, const int dn) {
    std::vector<int> ry;
    if (dn <= 0) { return ry; }

    int x[dn], n = dn;

    for (int i = 0; i < dn; i++) { x[i] = i; }
    for (int i = 0; i < index_num; i++) {
        int j = unif_index(n);
        ry.push_back(x[j]);
        x[j]  = x[--n];
    }
    return ry;
}

void fixupProb(DBL_T *p, int n, int require_k) {
    DBL_T sum  = 0.0, pi;
    int   npos = 0;

    for (int i = 0; i < n; i++, npos++) {
        pi = p[i];
        assert(0 <= pi && pi < INFINITY);
        sum += pi;
    }

    assert(npos != 0 && require_k <= npos);
    for (int i = 0; i < n; i++) { p[i] /= sum; }
}

/**
 * Sort a[] into descending order by "heapsort";
 * Sort ib[] alongside;
 * If initially, ib[] = 1...n, it will contain the permutation finally
 */
void revsort(DBL_T *a, int *ib, int n) {
    int   l, j, ir, i;
    DBL_T ra;
    int   ii;

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

int sample_int_index(const int dn, const DBL_T *prob) {
    int   perm[dn], nm1 = dn - 1;
    DBL_T prob_cpy[dn];

    for (int i = 0; i < dn; i++) {
        prob_cpy[i] = prob[i];
        perm[i]     = i + 1;
    }

    fixupProb(prob_cpy, dn, false);
    revsort(prob_cpy, perm, dn);

    DBL_T rT = unif_rand(), mass = 0;
    int      i;
    for (i = 0; i < nm1; i++) {
        mass += prob_cpy[i];
        if (rT <= mass) { break; }
    }
    return perm[i];
}

std::vector<int> sample_indices(
        const int sample_num,
        const std::vector<DBL_T> &prob,
        const bool replace
) {
    int              dn = (int) prob.size(), nm1 = dn - 1;
    assert(dn <= INT_MAX);

    if (!replace) { assert(sample_num <= dn); }

    DBL_T            prob_cpy[dn];
    int              i  = 0, j = 0;
    int              perm[dn];
    std::vector<int> res;

    for (const DBL_T p: prob) {
        prob_cpy[i] = p;
        perm[i]     = i;
        i++;
    }

    fixupProb(prob_cpy, dn, replace);
    revsort(prob_cpy, perm, dn);

    if (replace) {
        double rU;

        for (i = 1; i < dn; i++) { prob_cpy[i] += prob_cpy[i - 1]; }


        /* compute the sample */
        for (i = 0; i < sample_num; i++) {
            rU     = unif_rand();
            for (j = 0; j < nm1; j++) {
                if (rU <= prob_cpy[j]) { break; }
            }
            res.push_back(perm[j]);
        }
    } else {
        double rT, mass, totalmass = 1;
        int    k, n1;

        for (i = 0, n1 = dn - 1; i < sample_num; i++, n1--) {
            rT     = totalmass * unif_rand();
            mass   = 0;
            for (j = 0; j < n1; j++) {
                mass += prob_cpy[j];
                if (rT <= mass) { break; }
            }
            res.push_back(perm[j]);
            totalmass -= prob_cpy[j];
            for (k = j; k < n1; k++) {
                prob_cpy[k] = prob_cpy[k + 1];
                perm[k]     = perm[k + 1];
            }
        }
    }

    assert(res.size() == sample_num);
    return res;
}

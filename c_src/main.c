#include <stdio.h>

#include "util/math/math_ext.h"
#include "util/matrix/matrix.h"
#include "algorithm/scc.h"
#include "util/collection/collection.h"

void revsort(double *a, int *ib, int n)
{
/* Sort a[] into descending order by "heapsort";
 * sort ib[] alongside;
 * if initially, ib[] = 1...n, it will contain the permutation finally
 */

    int l, j, ir, i;
    double ra;
    int ii;

    if (n <= 1) return;

    a--; ib--;

    l = (n >> 1) + 1;
    ir = n;

    for (;;) {
        if (l > 1) {
            l = l - 1;
            ra = a[l];
            ii = ib[l];
        }
        else {
            ra = a[ir];
            ii = ib[ir];
            a[ir] = a[1];
            ib[ir] = ib[1];
            if (--ir == 1) {
                a[1] = ra;
                ib[1] = ii;
                return;
            }
        }
        i = l;
        j = l << 1;
        while (j <= ir) {
            if (j < ir && a[j] > a[j + 1]) ++j;
            if (ra > a[j]) {
                a[i] = a[j];
                ib[i] = ib[j];
                j += (i = j);
            }
            else
                j = ir + 1;
        }
        a[i] = ra;
        ib[i] = ii;
    }
}

void ProbSampleNoReplace(int n, double *p, int *perm,
                         int nans, int *ans) {
    double rU;
    int i, j;
    int nm1 = n - 1;

    /* record element identities */
    for (i = 0; i < n; i++)
        perm[i] = i + 1;

    /* sort the probabilities into descending order */
    revsort(p, perm, n);

    /* compute cumulative probabilities */
    for (i = 1 ; i < n; i++)
        p[i] += p[i - 1];

    /* compute the sample */
    for (i = 0; i < nans; i++) {
        rU = unif_rand();
        for (j = 0; j < nm1; j++) {
            if (rU <= p[j])
                break;
        }
        ans[i] = perm[j];
    }
}

int main() {
    set_seed(SEED);

    double prob1[6] = {0.1, 0.2, 0.1, 0.2, 0.1, 0.3};
    double prob2[6] = {0.1, 0.2, 0.1, 0.2, 0.1, 0.3};
    double prob3[6] = {0.1, 0.2, 0.1, 0.2, 0.1, 0.3};
    double prob4[6] = {0.1, 0.2, 0.1, 0.2, 0.1, 0.3};
    double prob5[6] = {0.1, 0.2, 0.1, 0.2, 0.1, 0.3};
//
//    printf("%d \n", sample_prob1(5, prob));
//    printf("%d \n", sample_prob1(5, prob));
//    printf("%d \n", sample_prob1(5, prob));
//    printf("%d \n", sample_prob1(5, prob));
//    printf("%d \n", sample_prob1(5, prob));

    int ans[1];
    int perm[6];
    ProbSampleNoReplace(6, prob1, perm, 1, ans);
    printf("[%d]\n", ans[0]);
    ProbSampleNoReplace(6, prob2, perm, 1, ans);
    printf("[%d]\n", ans[0]);
    ProbSampleNoReplace(6, prob3, perm, 1, ans);
    printf("[%d]\n", ans[0]);
    ProbSampleNoReplace(6, prob4, perm, 1, ans);
    printf("[%d]\n", ans[0]);
    ProbSampleNoReplace(6, prob5, perm, 1, ans);
    printf("[%d]\n", ans[0]);

//    sse_pars_t pars;
//    init_pars(&pars);
//    for (int i = 0; i < 1; ++i) {
//        calculate_scc(&pars, i);
//    }

    return 0;
}

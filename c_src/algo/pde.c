#include "algo.h"
#include "../matrix/matrix.h"

#include <stdio.h>

static const int PDE_THRESHOLD = 594;

void pde(
        const int t,
        const int idx,
        const sse_pars_t *pars,
        double f[Y_LEN][X_LEN],
        double m[Y_LEN][X_LEN]
) {
    if ((t + 1) > PDE_THRESHOLD) {
        for (int i = 1; i < Y_LEN - 1; ++i) {
            for (int j = 1; j < X_LEN - 1; ++j) {
                f[i][j] *= 1.0 - DT * pars->eta[idx] * m[i][j];
            }
        }
    }
}
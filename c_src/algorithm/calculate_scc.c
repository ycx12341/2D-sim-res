#include "scc.h"
#include "../util/math/math_ext.h"
#include "../util/matrix/matrix.h"

#include <stdio.h>
#include <math.h>

/* Space discretization */
#define H               (1.0 / 59.0)
#define SPACE_LENGTH_Y  ((1.0 / H) + 1.0)
#define SPACE_LENGTH_X  round(SPACE_LENGTH_Y * (280.0 / 480.0))

/* Time discretization */
#define T               4.52
#define DT              0.0025      // dt chosen as 0.0025 so the stability condition (dn < ((h^2)/4dt)) can be maintained.
#define TIME_STEPS      (T / DT)    // Total dimensionless timesteps and the dimensionless timesteps for one single day.
#define INT_TIME_STEPS  (1 / DT)
#define DAY_TIME_STEPS  600

#define BETA            0

double n0_to_m0(double x) {
    return 0.5 * x;
}

double n0_to_f0(double x) {
    return 1 - 0.5 * x;
}

/**
 * Purpose: Generate a SCC invasion pattern with the given parameters.
 * O(n^2)
 * @param pars Parameters
 */
void generate_pattern(sse_pars_t *pars) {
    /* Space discretization: Create a 60*35 domain. */
    int    x_len = (int) SPACE_LENGTH_X, y_len = (int) SPACE_LENGTH_Y;
    double x[x_len], y[y_len];
    if (seq(x, 0.0, H * (SPACE_LENGTH_X - 1), x_len) <= 0) { return; }
    if (seq(y, 0.0, 1, y_len) <= 0) { return; }

    /*
     * Initial condition
     */
    double n0[y_len][x_len], n[y_len][x_len], n0_sort[y_len][x_len];
    double f0[y_len][x_len], f[y_len][x_len];
    double m0[y_len][x_len], m[y_len][x_len];

    for (int j = 0; j < x_len; ++j) {
        n0[0][j] = x[j] <= 0.1 ? cos(M_PI * x[j] * 5) : 0;
    }

    MATRIX_ITR(y_len, x_len, n0[i][j] = n0[0][j];)

    MATRIX_MAP(n0, m0, y_len, x_len, n0_to_m0)
    MATRIX_MAP(n0, f0, y_len, x_len, n0_to_f0)

    MATRIX_COPY(n0, n, y_len, x_len)
    MATRIX_COPY(f0, f, y_len, x_len)
    MATRIX_COPY(m0, m, y_len, x_len)

//    MATRIX_PRINT(n, y_len, x_len, "[%.7f]")

    /* Sort the initial cells */
    MATRIX_COPY(n0, n0_sort, y_len, x_len)
}

/**
 * Calculates the least square difference between the simulated
 * invasion pattern and the reference invasion pattern.
 * @param pars Parameters
 */
void calculate_scc(sse_pars_t *pars) {
    generate_pattern(pars);
}
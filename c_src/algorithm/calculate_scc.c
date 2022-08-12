#include "scc.h"
#include "../util/matrix/matrix.h"
#include "../util/math/math_ext.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

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

/**
 * Purpose: Generate a SCC invasion pattern with the given parameters.
 * O(n^2)
 * @param pars Parameters
 */
void generate_pattern(sse_pars_t *pars, const int idx) {
    int i, j;

    /* Space discretization: Create a 60*35 domain. */
    int    x_len = (int) SPACE_LENGTH_X, y_len = (int) SPACE_LENGTH_Y;
    double x[x_len], y[y_len];
    assert(seq(x, 0.0, H * (SPACE_LENGTH_X - 1), x_len) == x_len);
    assert(seq(y, 0.0, 1, y_len) == y_len);

    /*
     * Initial condition
     */
    // double n0[y_len][x_len], n[y_len][x_len], n0_sort[y_len][x_len];
    // double f0[y_len][x_len], f[y_len][x_len];
    // double m0[y_len][x_len], m[y_len][x_len];

    double n0[y_len][x_len], n0_sort[y_len][x_len];
    double f0[y_len][x_len];
    double m0[y_len][x_len];

    for (j = 0; j < x_len; ++j) {
        n0[0][j] = x[j] <= 0.1 ? cos(M_PI * x[j] * 5) : 0;
    }

    MATRIX_ITR(i, y_len, j, x_len, {
        n0[i][j] = n0[0][j];
        m0[i][j] = 0.5 * n0[i][j];
        f0[i][j] = 1 - m0[i][j];
    })

    // MATRIX_MAP(n0, m0, y_len, x_len, n0_to_m0)
    // MATRIX_MAP(n0, f0, y_len, x_len, n0_to_f0)

    // MATRIX_COPY(n0, n, y_len, x_len)
    // MATRIX_COPY(f0, f, y_len, x_len)
    // MATRIX_COPY(m0, m, y_len, x_len)

    /* Sort the initial cells */
    MATRIX_COPY(n0, n0_sort, i, y_len, j, x_len)
//    MATRIX_PRINT(n0_sort, i, y_len, j, x_len, "[%.7f]")

    /* Initial glioma cells will be allocated to the locations with the highest densities in the domain (left boundary). */
    double      n_cells = round(SPACE_LENGTH_Y * pars->init_cells_cols[idx]);
    arraylist_t *coord  = new_arraylist(false);
    while (coord->size < n_cells) {
        int    sample_idx = 0;
        pair_t res        = matrix_max(y_len, x_len, n0_sort);
        if (res.y._int > 1) {
            sample_idx = unif_index(res.y._int);
        }
        pair_t sample = matrix_find(y_len, x_len, n0_sort, res.x._double, sample_idx + 1);

        node_t node;
        node._intPair[0] = sample.x._int;
        node._intPair[1] = sample.y._int;
        arraylist_append(coord, node);
        n0_sort[sample.x._int][sample.y._int] = -DBL_MAX;
    }

    double ind_position[y_len][x_len];
    for (i = 0; i < coord->size; ++i) {
        node_t pos  = arraylist_get(coord, i);
        ind_position[pos._intPair[0]][pos._intPair[1]] = 1;
    }
    // TODO: need test

    double ind_position_init[y_len][x_len];
    MATRIX_COPY(ind_position, ind_position_init, i, y_len, j, x_len)


    /* Set the cut points for the domain (to be used later for density matching & discrepancy calculation.) */
    int mat_size = SPACE_LENGTH_Y / 12;
    double x_cut[mat_size], y_cut[mat_size];
    assert(seq(&y_cut, 1, SPACE_LENGTH_Y, mat_size) == mat_size);
    assert(seq(&x_cut, 1, SPACE_LENGTH_X, mat_size) == mat_size);
    //  todo Test seq

    /* Numerical scheme which solves the PDE system */
    for (int k = 0; k < TIME_STEPS; k++) {

        /* At the end of everyday, some of the current cells in the domain will undergo extinction or mitosis. */
        if (k % DAY_TIME_STEPS == 0) {
            double cell_den[coord->size];
            for (int l = 0; l < coord->size; l++) {
                /* code */
            }
        }

        /* Solving the PDE model numerically */
        if (k > round(DAY_TIME_STEPS * (95.0 / 96.0))) {

        } else {

        }


    }
    


    // test print
    for (int k = 0; k < coord->size; ++k) {
        node_t node  = arraylist_get(coord, k);
        printf("[%d,%d]", node._intPair[0], node._intPair[1]);
    }

    arraylist_free(coord);
}

/**
 * Calculates the least square difference between the simulated
 * invasion pattern and the reference invasion pattern.
 * @param pars Parameters
 */
void calculate_scc(sse_pars_t *pars, const int idx) {
    generate_pattern(pars, idx);
}
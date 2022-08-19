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
#define Y_LEN           (int) SPACE_LENGTH_Y
#define X_LEN           (int) SPACE_LENGTH_X

/* Time discretization */
#define T               4.52
#define DT              0.0025      // dt chosen as 0.0025 so the stability condition (dn < ((h^2)/4dt)) can be maintained.
#define TIME_STEPS      (T / DT)    // Total dimensionless timesteps and the dimensionless timesteps for one single day.
#define INT_TIME_STEPS  (1 / DT)
#define DAY_TIME_STEPS  600

#define BETA            0

#define PDE_THRESHOLD   round(DAY_TIME_STEPS * (95.0 / 96.0))

/**
 * If the cell has more than two neighbouring positions which are
 * not occupied, it will proliferate. The original cell will vanish
 * and split into two daughter cells, which will be randomly
 * distributed into two unoccupied neighbouring locations.
 * @param nbr_num number of neighbouring positions
 */
void cell_proliferate(const int nbr_num, const int nbr_temp[nbr_num], const int *nbr_coord[2],
                      double ind_position[(int) SPACE_LENGTH_Y][(int) SPACE_LENGTH_X],
                      const int cell_position[2]) {
    int zeros = int_array_count(nbr_num, nbr_temp, 0);
    if (zeros >= 2) {
        int_arr_2_t sample_idx = unif_index2(zeros);

        const int *sample_a = nbr_coord[int_array_find(nbr_num, nbr_temp, 0, sample_idx.arr[0] + 1)];
        const int *sample_b = nbr_coord[int_array_find(nbr_num, nbr_temp, 0, sample_idx.arr[1] + 1)];

        ind_position[cell_position[0]][cell_position[1]] = 0;   // Original cell vanishes.
        ind_position[sample_a[0]][sample_a[1]]           = 1;   // Daughter cells being allocated.
        ind_position[sample_b[0]][sample_b[1]]           = 1;
    }
}

/**
 * Solving the PDE model numerically.
 */
void solve_PDE(const int idx, const int t, const sse_pars_t *pars,
               double n[(int) SPACE_LENGTH_Y][(int) SPACE_LENGTH_X],
               double f[(int) SPACE_LENGTH_Y][(int) SPACE_LENGTH_X],
               const double m[(int) SPACE_LENGTH_Y][(int) SPACE_LENGTH_X]) {
    if ((t + 1) > PDE_THRESHOLD) {
        /* Diffusion starts having an impact after a certain amount of time in day 1. */

        /*
         * f[i][j] = f[i][j] * ( 1 - (F_F) * m[i][j] )
         *      in range [1:Y_LEN-1, 1:X_LEN-1]
         *
         * F_F = DT * eta
         */
        const double F_F = DT * pars->eta[idx];
        for (int     i   = 1, is = Y_LEN - 1, js = X_LEN - 1; i < is; ++i) {
            for (int j = 1; j < js; ++j) {
                f[i][j] = f[i][j] * (1 - F_F * m[i][j]);
            }
        }

        /*
         * n[i][j] = n[i][j] * (A) + n[i][j+1] * (B) + n[i][j-1] * (C) + n[i-1][j] * (D) + n[i+1][j] * (E)
         *      in range [1:Y_LEN-1, 1:X_LEN-1]
         *
         * N_F1    = dt / (h^2)
         * N_F2    = dn * (N_F1)
         * N_F3    = gamma * (N_F1)
         *
         * A       = (A_F1) - ( (N_F3) * (A2) ) + (A_F3) * (A3)
         * A_F1    = 1 - ( 4 * (N_F2) )
         * A_F3    = r * dt
         * A2      = f[i][j+1] + f[i][j-1] + 4 * f[i][j] + f[i-1][j] + f[i+1][j]
         * A3      = 1 - n[i][j] - f[i][j]
         *
         * B       = (N_F2) - (N_F3 / 4) * ( f[i][j+1] - f[i][j-1] )
         * C       = (N_F2) + (N_F3 / 4) * ( f[i][j+1] - f[i][j-1] )
         * D       = (N_F2) - (N_F3 / 4) * ( f[i-1][j] - f[i+1][j] )
         * E       = (N_F2) + (N_F3 / 4) * ( f[i-1][j] - f[i+1][j] )
         */
        const double N_F1 = DT / pow(H, 2.0);
        const double N_F2 = pars->dn[idx] * N_F1;
        const double N_F3 = pars->gamma[idx] * N_F1;
        const double A_F1 = 1.0 - 4.0 * N_F2;
        const double A_F3 = pars->rn[idx] * DT;

        double n_cpy[Y_LEN][X_LEN];
        MATRIX_COPY(n, n_cpy, Y_LEN, X_LEN)

        double   A, A2, A3, B, C, D, E;
        for (int i = 1, is = Y_LEN - 1, js = X_LEN - 1; i < is; ++i) {
            for (int j = 1; j < js; ++j) {
                A2 = f[i][j + 1] + f[i][j - 1] + f[i - 1][j] + f[i + 1][j] - 4.0 * f[i][j];
                A3 = 1 - n[i][j] - f[i][j];
                A  = A_F1 - N_F3 * A2 + A_F3 * A3;
                B  = N_F2 - N_F3 / 4.0 * (f[i][j + 1] - f[i][j - 1]);
                C  = N_F2 + N_F3 / 4.0 * (f[i][j + 1] - f[i][j - 1]);
                D  = N_F2 - N_F3 / 4.0 * (f[i - 1][j] - f[i + 1][j]);
                E  = N_F2 + N_F3 / 4.0 * (f[i - 1][j] - f[i + 1][j]);

                n[i][j] = n_cpy[i][j] * A +
                          n_cpy[i][j + 1] * B + n_cpy[i][j - 1] * C +
                          n_cpy[i - 1][j] * D + n_cpy[i + 1][j] * E;
            }
        }
    } else {

    }
}

/**
 * Proliferation mechanism
 */
void proliferation(const int prof_cells_num, double prof_cells[prof_cells_num], arraylist_t *coord,
                   double ind_position[Y_LEN][X_LEN]) {
    for (int q = 0; q < prof_cells_num; q++) {
        node_t cell_pos         = arraylist_get(coord, (int) prof_cells[q]);
        int    cell_position[2] = {cell_pos._intPair[0], cell_pos._intPair[1]};

        /* Possible locations for daughter cells (8 surrounding points.) */
        const int right_position[2]      = {cell_position[0], cell_position[1] + 1};
        const int right_down_position[2] = {cell_position[0] + 1, cell_position[1] + 1};
        const int down_position[2]       = {cell_position[0] + 1, cell_position[1]};
        const int up_right_position[2]   = {cell_position[0] - 1, cell_position[1] + 1};
        const int up_position[2]         = {cell_position[0] - 1, cell_position[1]};
        const int up_left_position[2]    = {cell_position[0] - 1, cell_position[1] - 1};
        const int left_position[2]       = {cell_position[0], cell_position[1] - 1};
        const int left_down_position[2]  = {cell_position[0] + 1, cell_position[1] - 1};

        const int right      = (int) ind_position[right_position[0]][right_position[1]];
        const int right_down = (int) ind_position[right_down_position[0]][right_down_position[1]];
        const int down       = (int) ind_position[down_position[0]][down_position[1]];
        const int up_right   = (int) ind_position[up_right_position[0]][up_right_position[1]];
        const int up         = (int) ind_position[up_position[0]][up_position[1]];
        const int up_left    = (int) ind_position[up_left_position[0]][up_left_position[1]];
        const int left       = (int) ind_position[left_position[0]][left_position[1]];
        const int left_down  = (int) ind_position[left_down_position[0]][left_down_position[1]];

        if (cell_position[0] == 0) {

            if (cell_position[1] == 0) {                    /* Special case: top left corner */
                /* Possible directions to move, check if these points are occupied. */
                const int neighbouring_temp[3]   = {right, right_down, down};
                const int *neighbouring_coord[3] = {right_position, right_down_position, down_position};
                cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, cell_position);
            } else if (cell_position[1] == X_LEN - 1) {     /* Special case: top right corner */
                const int neighbouring_temp[3]   = {left, left_down, down};
                const int *neighbouring_coord[3] = {left_position, left_down_position, down_position};
                cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, cell_position);
            } else {
                const int neighbouring_temp[5]   = {left, right, down, left_down, right_down};
                const int *neighbouring_coord[5] = {left_position, right_position, down_position,
                                                    left_down_position, right_down_position};
                cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, cell_position);
            }

        } else if (cell_position[0] == Y_LEN - 1) {         /* Lower boundary */

            if (cell_position[1] == 0) {                    /* Special case: bottom left corner */
                const int neighbouring_temp[3]   = {up, up_right, right};
                const int *neighbouring_coord[3] = {up_position, up_right_position, right_position};
                cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, cell_position);
            } else if (cell_position[1] == X_LEN - 1) {     /* Special case: bottom right corner */
                const int neighbouring_temp[3]   = {left, up, up_left};
                const int *neighbouring_coord[3] = {left_position, up_position, up_left_position};
                cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, cell_position);
            } else {                                        // WARNING: BE CAREFUL with the order (keep for consistency)
                const int neighbouring_temp[5]   = {left, up_left, up, up_right, right};
                const int *neighbouring_coord[5] = {left_position, up_left_position, up_position,
                                                    up_right_position, right_position};
                cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, cell_position);
            }

        } else if (cell_position[1] == 0) {                 /* Left boundary */

            if (cell_position[0] == 0) {                    // do nothing so far
            } else if (cell_position[0] == Y_LEN - 1) {     // do nothing so far
            } else {
                const int neighbouring_temp[5]   = {up, up_right, right, right_down, down};
                const int *neighbouring_coord[5] = {
                        up_position, up_right_position, right_position, right_down_position, down_position
                };
                cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, cell_position);
            }

        } else if (cell_position[1] == X_LEN - 1) {         /* Right boundary */

            if (cell_position[0] == 0) {
            } else if (cell_position[0] == Y_LEN - 1) {
            } else {
                const int neighbouring_temp[5]   = {up, up_left, left, left_down, down};
                const int *neighbouring_coord[5] = {
                        up_position, up_left_position, left_position, left_down_position,
                        down_position
                };
                cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, cell_position);
            }

        } else {                                            /* Other cells that are not at the boundary */
            const int neighbouring_temp[8]   = {left, right, up, down, up_left, up_right, left_down,
                                                right_down};
            const int *neighbouring_coord[8] = {
                    left_position, right_position, up_position, down_position,
                    up_left_position, up_right_position, left_down_position, right_down_position
            };
            cell_proliferate(8, neighbouring_temp, neighbouring_coord, ind_position, cell_position);
        }
    }
}

/**
 * Purpose: Generate a SCC invasion pattern with the given parameters.
 * @param pars Parameters
 */
double generate_pattern(sse_pars_t *pars, const int idx) {
    /* Space discretization: Create a 60*35 domain. */
    double x[X_LEN], y[Y_LEN];
    assert(seq_length_out(x, 0.0, H * (X_LEN - 1), X_LEN) == X_LEN);
    assert(seq_length_out(y, 0.0, 1, Y_LEN) == Y_LEN);

    /*
     * Initial condition
     */
    double n0[Y_LEN][X_LEN], n0_sort[Y_LEN][X_LEN], n[Y_LEN][X_LEN];
    double f0[Y_LEN][X_LEN], f[Y_LEN][X_LEN];
    double m0[Y_LEN][X_LEN], m[Y_LEN][X_LEN];

    for (int j = 0; j < X_LEN; ++j) {
        n0[0][j] = x[j] <= 0.1 ? cos(M_PI * x[j] * 5) : 0;
    }

    MATRIX_ITR(Y_LEN, X_LEN, {
        n0[_i_][_j_] = n0[0][_j_];
        m0[_i_][_j_] = 0.5 * n0[_i_][_j_];
        f0[_i_][_j_] = 1 - m0[_i_][_j_];
    })

    MATRIX_COPY(n0, n, Y_LEN, X_LEN)
    MATRIX_COPY(f0, f, Y_LEN, X_LEN)
    MATRIX_COPY(m0, m, Y_LEN, X_LEN)

    /* Sort the initial cells */
    MATRIX_COPY(n0, n0_sort, Y_LEN, X_LEN)

    /* Initial glioma cells will be allocated to the locations with the highest densities in the domain (left boundary). */
    double      n_cells = round(SPACE_LENGTH_Y * pars->init_cells_cols[idx]);
    arraylist_t *coord  = new_arraylist(false);
    while (coord->size < n_cells) {
        int    sample_idx = 0;
        pair_t res        = matrix_max(Y_LEN, X_LEN, n0_sort);

        if (res.y._int > 1) { sample_idx = unif_index(res.y._int); }

        pair_t sample = matrix_find(Y_LEN, X_LEN, n0_sort, res.x._double, sample_idx + 1);
        assert(sample.x._int >= 0 && sample.y._int >= 0);

        node_t node;
        node._intPair[0] = sample.x._int;
        node._intPair[1] = sample.y._int;
        arraylist_append(coord, node);
        n0_sort[sample.x._int][sample.y._int] = -DBL_MAX;
    }

    double   ind_position[Y_LEN][X_LEN];
    for (int i = 0; i < coord->size; ++i) {
        node_t pos = arraylist_get(coord, i);
        ind_position[pos._intPair[0]][pos._intPair[1]] = 1;
    }

    double ind_position_init[Y_LEN][X_LEN];
    MATRIX_COPY(ind_position, ind_position_init, Y_LEN, X_LEN)

    /* Set the cut points for the domain (to be used later for density matching & discrepancy calculation.) */
    int    mat_size  = SPACE_LENGTH_Y / 12;
    int    y_cut_len = ceil((SPACE_LENGTH_Y - 1.0) / mat_size);
    int    x_cut_len = ceil((SPACE_LENGTH_X - 1.0) / mat_size);
    double x_cut[x_cut_len], y_cut[y_cut_len];
    assert(seq_by(y_cut, 1, SPACE_LENGTH_Y, mat_size) == y_cut_len);
    assert(seq_by(x_cut, 1, SPACE_LENGTH_X, mat_size) == x_cut_len);

    /* Numerical scheme which solves the PDE system */
    for (int t = 0; t < TIME_STEPS; t++) {

        /* At the end of every day, some current cells in the domain will undergo extinction or mitosis. */
        if ((t + 1) % DAY_TIME_STEPS == 0) {

            /* Extract the density values at the positions of current cells */
            int    cell_den_len = coord->size;
            double cell_den[cell_den_len];

            for (int i = 0; i < cell_den_len; i++) {
                node_t pos = arraylist_get(coord, i);
                cell_den[i] = n[pos._intPair[0]][pos._intPair[1]];
            }

            /* Cell extinction: some the current cells at the locations with the
             * lowest cell densities will undergo extinction. */
            int dead_cells_num = (int) round((double) coord->size * pars->prob_death[idx]);
            int dead_cells[dead_cells_num];

            for (int uu = 0; uu < dead_cells_num; ++uu) {
                /* Find the cells at locations with the lowest densities. */
                int    sample_idx = 0;
                pair_t ind_dead   = array_min(cell_den_len, cell_den);
                if (ind_dead.y._int > 1) {
                    sample_idx = unif_index(ind_dead.y._int);
                }
                int sample = double_array_find(cell_den_len, cell_den, ind_dead.x._double, sample_idx + 1);
                assert(sample >= 0);
                dead_cells[uu]   = sample;
                cell_den[sample] = DBL_MAX;
            }

            /* Update the cell coordinates and the density vector. */
            coord        = arraylist_remove_many(coord, dead_cells_num, dead_cells);
            cell_den_len = double_array_delete_many(cell_den_len, cell_den, dead_cells_num, dead_cells);

            /* If all cells are dead, terminate the algorithm. */
            if (cell_den_len == 0) {
                return NAN;
            }

            /* Cell mitosis: some current cells at the locations
             * with the highest densities will undergo mitosis. */
            int    prof_cells_num = (int) round((double) coord->size * pars->prob_prof[idx]);
            double prof_cells[prof_cells_num];

            for (int uu = 0; uu < prof_cells_num; ++uu) {
                /* Find the cells at the locations with the highest densities. */
                int    sample_idx = 0;
                pair_t ind_prof   = array_max(cell_den_len, cell_den);
                if (ind_prof.y._int > 1) {
                    sample_idx = unif_index(ind_prof.y._int);
                }
                int sample = double_array_find(cell_den_len, cell_den, ind_prof.x._double, sample_idx + 1);
                assert(sample >= 0);
                prof_cells[uu]   = sample;
                cell_den[sample] = -DBL_MAX;
            }

            /* Record the current positions of the remaining cells. */
            MATRIX_INIT(ind_position, Y_LEN, X_LEN, 0)
            for (int i  = 0; i < coord->size; ++i) {
                node_t pos = arraylist_get(coord, i);
                ind_position[pos._intPair[0]][pos._intPair[1]] = 1;
            }

            /* Proliferation mechanism */
            proliferation(prof_cells_num, prof_cells, coord, ind_position);

            /* Update the cell coordinates and positions.*/
            arraylist_clear(coord);
            MATRIX_ITR2(Y_LEN, X_LEN, {
                if (ind_position[_i_][_j_] == 1) {
                    node_t pos;
                    pos._intPair[0] = _i_;
                    pos._intPair[1] = _j_;
                    arraylist_append(coord, pos);
                } else {
                    ind_position[_i_][_j_] = 0;
                }
            })
        }

        /* Solving the PDE model numerically */
        solve_PDE(idx, t, pars, n, f, m);

        // Testing
        if ((t + 1) > PDE_THRESHOLD) {
            MATRIX_PRINT(n, Y_LEN, X_LEN, " %.7f ");
            return NAN;
        }

    }

    arraylist_free(coord);
}
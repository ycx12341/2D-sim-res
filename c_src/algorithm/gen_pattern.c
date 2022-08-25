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
#define MAT_SIZE        (int) (SPACE_LENGTH_Y / 12)
#define Y_CUT_LEN       (int) ceil((SPACE_LENGTH_Y - 1.0) / MAT_SIZE)
#define X_CUT_LEN       (int) ceil((SPACE_LENGTH_X - 1.0) / MAT_SIZE)

/* Time discretization */
#define T               4.52
#define DT              0.0025      // dt chosen as 0.0025 so the stability condition (dn < ((h^2)/4dt)) can be maintained.
#define TIME_STEPS      (T / DT)    // Total dimensionless timesteps and the dimensionless timesteps for one single day.
#define INT_TIME_STEPS  (1 / DT)
#define DAY_TIME_STEPS  600
#define N_F1            DT / pow(H, 2.0)

#define BETA            0

#define PDE_THRESHOLD   round(DAY_TIME_STEPS * (95.0 / 96.0))

/**
 * If the cell has more than two neighbouring positions which are
 * not occupied, it will proliferate. The original cell will vanish
 * and split into two daughter cells, which will be randomly
 * distributed into two unoccupied neighbouring locations.
 * @param nbr_num number of neighbouring positions
 */
void cell_proliferate(
        const int nbr_num,
        const int nbr_temp[nbr_num],
        const int *nbr_coord[2],
        double ind_position[Y_LEN][X_LEN],
        const int x_coord,
        const int y_coord
) {
    const int zeros = int_array_count(nbr_num, nbr_temp, 0);
    if (zeros >= 2) {
        int_arr_2_t sample_idx = unif_index2(zeros);

        const int *sample_a = nbr_coord[int_array_find(nbr_num, nbr_temp, 0, sample_idx.arr[0] + 1)];
        const int *sample_b = nbr_coord[int_array_find(nbr_num, nbr_temp, 0, sample_idx.arr[1] + 1)];

        ind_position[x_coord][y_coord]         = 0;     // Original cell vanishes.
        ind_position[sample_a[0]][sample_a[1]] = 1;     // Daughter cells being allocated.
        ind_position[sample_b[0]][sample_b[1]] = 1;
    }
}

/**
 * Solving the PDE model numerically.
 */
void solve_PDE(
        const int idx,
        const int t,
        const sse_pars_t *pars,
        double n[Y_LEN][X_LEN],
        double f[Y_LEN][X_LEN],
        double m[Y_LEN][X_LEN]
) {
    /*
     * f[i][j] = f[i][j] * ( 1 - (F_F) * m[i][j] )
     *      in range [1:Y_LEN-1, 1:X_LEN-1]
     *
     * F_F = DT * eta
     */
    const double F_F1 = DT * pars->eta[idx];
    for (int     i    = 1, is = Y_LEN - 1, js = X_LEN - 1; i < is; ++i) {
        for (int j = 1; j < js; ++j) {
            f[i][j] = f[i][j] * (1 - F_F1 * m[i][j]);
        }
    }

    /*
     * G_F1    = r * dt
     * G_F2    = dt * alpha
     */
    const double G_F1 = pars->rn[idx] * DT;
    const double G_F2 = DT * pars->alpha[idx];

    if ((t + 1) > PDE_THRESHOLD) {
        /* Diffusion starts having an impact after a certain amount of time in day 1. */

        /*
         * n[i][j] = n[i][j] * (A) + n[i][j+1] * (B) + n[i][j-1] * (C) + n[i-1][j] * (D) + n[i+1][j] * (E)
         *      in range [1:Y_LEN-1, 1:X_LEN-1]
         *
         * N_F1    = dt / (h^2)
         * N_F2    = dn * (N_F1)
         * N_F3    = gamma * (N_F1)
         * N_F4    = N_F3 / 4.0;
         *
         * A       = (A_F1) - ( (N_F3) * (A2) ) + (G_F1) * (A3)
         * A_F1    = 1 - ( 4 * (N_F2) )
         * A2      = f[i][j+1] + f[i][j-1] + 4 * f[i][j] + f[i-1][j] + f[i+1][j]
         * A3      = 1 - n[i][j] - f[i][j]
         *
         * B       = (N_F2) - (N_F3 / 4) * ( f[i][j+1] - f[i][j-1] )
         * C       = (N_F2) + (N_F3 / 4) * ( f[i][j+1] - f[i][j-1] )
         * D       = (N_F2) - (N_F3 / 4) * ( f[i-1][j] - f[i+1][j] )
         * E       = (N_F2) + (N_F3 / 4) * ( f[i-1][j] - f[i+1][j] )
         */
        const double N_F2 = pars->dn[idx] * N_F1;
        const double N_F3 = pars->gamma[idx] * N_F1;
        const double N_F4 = DT / pow(H, 2.0) * pars->gamma[idx] / 4.0;
        const double A_F1 = 1.0 - 4.0 * N_F2;

        double n_cpy[Y_LEN][X_LEN];
        MATRIX_COPY(n, n_cpy, Y_LEN, X_LEN)

        double   A, A2, A3, B, C, D, E;
        for (int i = 1, is = Y_LEN - 1, js = X_LEN - 1; i < is; ++i) {
            for (int j = 1; j < js; ++j) {
                A2 = f[i][j + 1] + f[i][j - 1] - 4 * f[i][j] + f[i - 1][j] + f[i + 1][j];
                A3 = 1 - n_cpy[i][j] - f[i][j];
                A  = A_F1 - N_F3 * A2 + G_F1 * A3;
                B  = N_F2 - N_F4 * (f[i][j + 1] - f[i][j - 1]);
                C  = N_F2 + N_F4 * (f[i][j + 1] - f[i][j - 1]);
                D  = N_F2 - N_F4 * (f[i - 1][j] - f[i + 1][j]);
                E  = N_F2 + N_F4 * (f[i - 1][j] - f[i + 1][j]);

                n[i][j] = n_cpy[i][j] * A +
                          n_cpy[i][j + 1] * B + n_cpy[i][j - 1] * C +
                          n_cpy[i - 1][j] * D + n_cpy[i + 1][j] * E;
            }
        }

        /*
         * m[i][j] = m[i][j] * (M_F2) + (G_F2) * n[i][j] + (M_F1) * (F)
         *      in range [1:Y_LEN-1, 1:X_LEN-1]
         *
         * M_F1    = (N_F1) * dm
         * M_F2    = 1 - 4 * M_F1
         * F       = m[i][j+1] + m[i][j-1] + m[i-1][j] + m[i+1][j]
         */
        const double M_F1 = N_F1 * pars->dm[idx];
        const double M_F2 = 1.0 - 4.0 * M_F1;

        double m_cpy[Y_LEN][X_LEN];
        MATRIX_COPY(m, m_cpy, Y_LEN, X_LEN)

        double   F;
        for (int i = 1, is = Y_LEN - 1, js = X_LEN - 1; i < is; ++i) {
            for (int j = 1; j < js; ++j) {
                F = m_cpy[i][j + 1] + m_cpy[i][j - 1] + m_cpy[i - 1][j] + m_cpy[i + 1][j];
                m[i][j] = m[i][j] * M_F2 + G_F2 * n[i][j] + M_F1 * F;
            }
        }
    } else {
        /*
         * n[i][j] = n[i][j] * (1 + G_F1 * (1 - n[i][j] - f[i][j]))
         * m[i][j] = m[i][j] + G_F2 * n[i][j]
         */
        for (int i = 1, is = Y_LEN - 1, js = X_LEN - 1; i < is; ++i) {
            for (int j = 1; j < js; ++j) {
                n[i][j] = n[i][j] * (1 + G_F1 * (1 - n[i][j] - f[i][j]));
                m[i][j] = m[i][j] + G_F2 * n[i][j];
            }
        }
    }
}

/**
 * Boundary condition
 */
void boundary_condition(
        double n[Y_LEN][X_LEN],
        double f[Y_LEN][X_LEN],
        double m[Y_LEN][X_LEN]
) {
    for (int j = 0, x_1 = Y_LEN - 1, x_2 = x_1 - 1; j < X_LEN; ++j) {
        n[0][j]   = n[1][j];        // n[1, ] <- n[2, ]
        f[0][j]   = f[1][j];        // f[1, ] <- f[2, ]
        m[0][j]   = m[1][j];        // m[1, ] <- m[2, ]
        n[x_1][j] = n[x_2][j];      // n[length(y), ] <- n[(length(y) - 1), ]
        f[x_1][j] = f[x_2][j];      // f[length(y), ] <- f[(length(y) - 1), ]
        m[x_1][j] = m[x_2][j];      // m[length(y), ] <- m[(length(y) - 1), ]
    }
    for (int i = 0, y_1 = X_LEN - 1, y_2 = y_1 - 1; i < Y_LEN; ++i) {
        n[i][0]   = n[i][1];        // n[, 1] <- n[, 2]
        f[i][0]   = f[i][1];        // f[, 1] <- f[, 2]
        m[i][0]   = m[i][1];        // m[, 1] <- m[, 2]
        n[i][y_1] = n[i][y_2];      // n[, length(x)] <- n[, (length(x) - 1)]
        f[i][y_1] = f[i][y_2];      // f[, length(x)] <- f[, (length(x) - 1)]
        m[i][y_1] = m[i][y_2];      // m[, length(x)] <- m[, (length(x) - 1)]
    }
}

/**
 * Proliferation mechanism
 */
void proliferation(
        const int prof_cells_num,
        double prof_cells[prof_cells_num],
        arraylist_t *coord,
        double ind_position[Y_LEN][X_LEN]
) {
    for (int q = 0; q < prof_cells_num; q++) {
        node_t    cell_pos = arraylist_get(coord, (int) prof_cells[q]);
        const int x_coord  = cell_pos._intPair[0];
        const int y_coord  = cell_pos._intPair[1];

        /* Possible locations for daughter cells (8 surrounding points.) */
        const int right_position[2]      = {x_coord, y_coord + 1};
        const int right_down_position[2] = {x_coord + 1, y_coord + 1};
        const int down_position[2]       = {x_coord + 1, y_coord};
        const int up_right_position[2]   = {x_coord - 1, y_coord + 1};
        const int up_position[2]         = {x_coord - 1, y_coord};
        const int up_left_position[2]    = {x_coord - 1, y_coord - 1};
        const int left_position[2]       = {x_coord, y_coord - 1};
        const int left_down_position[2]  = {x_coord + 1, y_coord - 1};

        const int right      = (int) ind_position[right_position[0]][right_position[1]];
        const int right_down = (int) ind_position[right_down_position[0]][right_down_position[1]];
        const int down       = (int) ind_position[down_position[0]][down_position[1]];
        const int up_right   = (int) ind_position[up_right_position[0]][up_right_position[1]];
        const int up         = (int) ind_position[up_position[0]][up_position[1]];
        const int up_left    = (int) ind_position[up_left_position[0]][up_left_position[1]];
        const int left       = (int) ind_position[left_position[0]][left_position[1]];
        const int left_down  = (int) ind_position[left_down_position[0]][left_down_position[1]];

        if (x_coord == 0) {
            if (y_coord == 0) {                    /* Special case: top left corner */
                /* Possible directions to move, check if these points are occupied. */
                const int neighbouring_temp[3]   = {right, right_down, down};
                const int *neighbouring_coord[3] = {right_position, right_down_position, down_position};
                cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
            } else if (y_coord == X_LEN - 1) {     /* Special case: top right corner */
                const int neighbouring_temp[3]   = {left, left_down, down};
                const int *neighbouring_coord[3] = {left_position, left_down_position, down_position};
                cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
            } else {
                const int neighbouring_temp[5]   = {left, right, down, left_down, right_down};
                const int *neighbouring_coord[5] = {left_position, right_position, down_position,
                                                    left_down_position, right_down_position};
                cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
            }
        } else if (x_coord == Y_LEN - 1) {         /* Lower boundary */
            if (y_coord == 0) {                    /* Special case: bottom left corner */
                const int neighbouring_temp[3]   = {up, up_right, right};
                const int *neighbouring_coord[3] = {up_position, up_right_position, right_position};
                cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
            } else if (y_coord == X_LEN - 1) {     /* Special case: bottom right corner */
                const int neighbouring_temp[3]   = {left, up, up_left};
                const int *neighbouring_coord[3] = {left_position, up_position, up_left_position};
                cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
            } else {                                        // WARNING: BE CAREFUL with the order (keep for consistency)
                const int neighbouring_temp[5]   = {left, up_left, up, up_right, right};
                const int *neighbouring_coord[5] = {left_position, up_left_position, up_position,
                                                    up_right_position, right_position};
                cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
            }
        } else if (y_coord == 0 && x_coord != Y_LEN - 1) {              /* Left boundary */
            const int neighbouring_temp[5]   = {up, up_right, right, right_down, down};
            const int *neighbouring_coord[5] = {
                    up_position, up_right_position, right_position, right_down_position, down_position
            };
            cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
        } else if (y_coord == X_LEN - 1 && x_coord != Y_LEN - 1) {      /* Right boundary */
            const int neighbouring_temp[5]   = {up, up_left, left, left_down, down};
            const int *neighbouring_coord[5] = {
                    up_position, up_left_position, left_position, left_down_position,
                    down_position
            };
            cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
        } else {                                                        /* Other cells that are not at the boundary */
            const int neighbouring_temp[8]   = {
                    left, right, up, down, up_left, up_right, left_down, right_down
            };
            const int *neighbouring_coord[8] = {
                    left_position, right_position, up_position, down_position,
                    up_left_position, up_right_position, left_down_position, right_down_position
            };
            cell_proliferate(8, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
        }
    }
}

bool cell_density(
        const sse_pars_t *pars,
        const int idx,
        arraylist_t *coord,
        double n[Y_LEN][X_LEN],
        double ind_position[Y_LEN][X_LEN]
) {
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
    arraylist_remove_many(coord, dead_cells_num, dead_cells);
    cell_den_len = double_array_delete_many(cell_den_len, cell_den, dead_cells_num, dead_cells);

    /* If all cells are dead, terminate the algorithm. */
    if (cell_den_len == 0 || coord->size == 0) { return false; }

    /* Cell mitosis: some current cells at the locations
             * with the highest densities will undergo mitosis. */
    const int prof_cells_num = (int) round((double) coord->size * pars->prob_prof[idx]);
    double    prof_cells[prof_cells_num];

    for (int uu = 0; uu < prof_cells_num; ++uu) {
        /* Find the cells at the locations with the highest densities. */
        int    sample_idx = 0;
        pair_t ind_prof   = array_max(cell_den_len, cell_den);
        if (ind_prof.y._int > 1) {
            sample_idx = unif_index(ind_prof.y._int);
        }
        const int sample = double_array_find(cell_den_len, cell_den, ind_prof.x._double, sample_idx + 1);
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

    return true;
}


#define MOVEMENT_F_SUM(f_ip1j, f_im1j, f_ijp1, f_ijm1, f_xy) \
    ((f_ip1j) + (f_im1j) + (f_ijp1) + (f_ijm1) - 4.0 * (f_xy))
#define MAX_MOVEMENT 5

void movement(
        const int t,
        const int idx,
        const sse_pars_t *pars,
        arraylist_t *coord,
        double n[Y_LEN][X_LEN],
        double f[Y_LEN][X_LEN],
        double ind_position[Y_LEN][X_LEN]
) {
    if ((t + 1) <= PDE_THRESHOLD) { return; }

    /* Movements only occur after diffusion starts having an impact on cells. */
    const double G_F1 = pars->rn[idx] * DT;
    const double N_F2 = pars->dn[idx] * N_F1;
    const double N_F3 = pars->gamma[idx] * N_F1;
    const double N_F4 = N_F3 / 4.0;
    const double A_F1 = 1.0 - 4.0 * N_F2;

    for (int i = 0, len = coord->size; i < len; ++i) {
        node_t    cell_pos = arraylist_get(coord, i);
        const int x_coord  = cell_pos._intPair[0];
        const int y_coord  = cell_pos._intPair[1];

        const double A = 1 - n[x_coord][y_coord] - f[x_coord][y_coord];

        double p0, p1 = NAN, p2 = NAN, p3 = NAN, p4 = NAN;
        double F;
        double f_ip1j = f[x_coord][y_coord + 1];
        double f_im1j = f[x_coord][y_coord - 1];
        double f_ijp1 = f[x_coord - 1][y_coord];
        double f_ijm1 = f[x_coord + 1][y_coord];

        if (y_coord == 0) {
            f_im1j = 0;
            p1     = 0;
            if (x_coord == 0) {
                f_ijp1 = 0;
                p4     = 0;
            } else if (x_coord == Y_LEN - 1) {
                f_ijm1 = 0;
                p3     = 0;
            }
        } else if (y_coord == X_LEN - 1) {
            f_ip1j = 0;
            p2     = 0;
            if (x_coord == 0) {
                f_ijp1 = 0;
                p4     = 0;
            } else if (x_coord == Y_LEN - 1) {
                f_ijm1 = 0;
                p3     = 0;
            }
        } else if (x_coord == 0) {
            f_ijp1 = 0;
            p4     = 0;
        } else if (x_coord == Y_LEN - 1) {
            f_ijm1 = 0;
            p3     = 0;
        }

        F  = MOVEMENT_F_SUM(f_ip1j, f_im1j, f_ijp1, f_ijm1, f[x_coord][y_coord]);
        p0 = A_F1 - N_F3 * F + G_F1 * A;
        p1 = isnan(p1) ? N_F2 - N_F4 * (f_ip1j - f_im1j) : 0;
        p2 = isnan(p2) ? N_F2 + N_F4 * (f_ip1j - f_im1j) : 0;
        p3 = isnan(p3) ? N_F2 - N_F4 * (f_ijp1 - f_ijm1) : 0;
        p4 = isnan(p4) ? N_F2 + N_F4 * (f_ijp1 - f_ijm1) : 0;

        double p[MAX_MOVEMENT] = {p0, p1, p2, p3, p4};
        int    zeros           = 0;

        for (int j = 0; j < MAX_MOVEMENT; ++j) {
            if (p[j] < 0) {
                p[j] = 0;
            } else if (p[j] > 1) {
                p[j] = 1;
            }
            if (p[j] == 0) { zeros++; }
        }

        int mvment = zeros < MAX_MOVEMENT ? sample_prob1(MAX_MOVEMENT, p) : 0;

        assert(0 <= mvment && mvment <= 5);
        if (mvment == 1) {
            ind_position[x_coord][y_coord] = 1;
        } else if (mvment == 2) {
            if (ind_position[x_coord][y_coord - 1] == 0) {
                ind_position[x_coord][y_coord] = 0;
                node_t val;
                val._intPair[0] = x_coord;
                val._intPair[1] = y_coord - 1;
                arraylist_set(coord, i, val);
                ind_position[x_coord][y_coord - 1] = 1;
            } else {
                ind_position[x_coord][y_coord] = 1;
            }
        } else if (mvment == 3) {
            if (ind_position[x_coord][y_coord + 1] == 0) {
                ind_position[x_coord][y_coord] = 0;
                node_t val;
                val._intPair[0] = x_coord;
                val._intPair[1] = y_coord + 1;
                arraylist_set(coord, i, val);
                ind_position[x_coord][y_coord + 1] = 1;
            } else {
                ind_position[x_coord][y_coord] = 1;
            }
        } else if (mvment == 4) {
            if (ind_position[x_coord + 1][y_coord] == 0) {
                ind_position[x_coord][y_coord] = 0;
                node_t val;
                val._intPair[0] = x_coord + 1;
                val._intPair[1] = y_coord;
                arraylist_set(coord, i, val);
                ind_position[x_coord + 1][y_coord] = 1;
            } else {
                ind_position[x_coord][y_coord] = 1;
            }
        } else if (mvment == 5) {
            if (ind_position[x_coord - 1][y_coord] == 0) {
                ind_position[x_coord][y_coord] = 0;
                node_t val;
                val._intPair[0] = x_coord - 1;
                val._intPair[1] = y_coord;
                arraylist_set(coord, i, val);
                ind_position[x_coord - 1][y_coord] = 1;
            } else {
                ind_position[x_coord][y_coord] = 1;
            }
        }
    }
}

#undef MOVEMENT_F_SUM

/**
 * Purpose: Generate a SCC invasion pattern with the given parameters.
 * @param pars Parameters
 */
double generate_pattern(const sse_pars_t *pars, const int idx) {
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
    arraylist_t *coord  = new_arraylist();
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

    double ind_position[Y_LEN][X_LEN];
    MATRIX_INIT(ind_position, Y_LEN, X_LEN, 0)

    for (int i = 0; i < coord->size; ++i) {
        node_t pos = arraylist_get(coord, i);
        ind_position[pos._intPair[0]][pos._intPair[1]] = 1;
    }

    double ind_position_init[Y_LEN][X_LEN];
    MATRIX_COPY(ind_position, ind_position_init, Y_LEN, X_LEN)

    /* Set the cut points for the domain (to be used later for density matching & discrepancy calculation.) */
    double x_cut[X_CUT_LEN], y_cut[Y_CUT_LEN];
    assert(seq_by(y_cut, 1, SPACE_LENGTH_Y, MAT_SIZE) == Y_CUT_LEN);
    assert(seq_by(x_cut, 1, SPACE_LENGTH_X, MAT_SIZE) == X_CUT_LEN);

    /* Numerical scheme which solves the PDE system */
    for (int t = 0; t < TIME_STEPS; t++) {

        /* At the end of every day, some current cells in the domain will undergo extinction or mitosis. */
        if ((t + 1) % DAY_TIME_STEPS == 0
            && !cell_density(pars, idx, coord, n, ind_position)) {
            arraylist_free(coord);
            return NAN;
        }

        solve_PDE(idx, t, pars, n, f, m);   // Solving the PDE model numerically

        boundary_condition(n, f, m);        // Boundary condition

        /* If singularity occurs while solving the PDE model, terminate the simulation. */
        double n_itr, f_itr, m_itr;
        MATRIX_ITR(Y_LEN, X_LEN, {
            n_itr = n[_i_][_j_];
            f_itr = f[_i_][_j_];
            m_itr = m[_i_][_j_];
            if (n_itr < 0 || isnan(n_itr) ||
                f_itr < 0 || isnan(f_itr) ||
                m_itr < 0 || isnan(m_itr)) {
                arraylist_free(coord);
                return NAN;
            }
        })

        movement(t, idx, pars, coord, n, f, ind_position);

        if ((t + 1) == (DAY_TIME_STEPS * 3)) {  // TODO replace 3

            for (int i = 0; i < coord->size; ++i) {
                node_t node = arraylist_get(coord, i);
                printf("%d ( %d, %d )\n", i, node._intPair[0], node._intPair[1]);
            }

//            double   density_mat[Y_CUT_LEN][X_CUT_LEN];
//            for (int i = 0; i < Y_CUT_LEN; ++i) {
//                for (int j = 0; j < X_CUT_LEN; ++j) {
//
//                    int ones = 0;
//                    int yc   = (int) y_cut[i], xc = (int) x_cut[j];

//                    for (int k = yc - 1; k < yc + MAT_SIZE - 1; ++k) {
//                        for (int l = xc - 1; l < xc + MAT_SIZE - 1; ++l) {
//                            if (ind_position[k][l] == 1) ones++;
//                            printf("%f ", ind_position[k][l]);
//                        }
//                        printf("\n");
//                    }
//                    printf("ones: %d\n", ones);
//                    printf("%d %d %d %d\n", yc, yc + MAT_SIZE - 1, xc, xc + MAT_SIZE - 1);
//                }
//            }
            return t;
        }
    }

    arraylist_free(coord);
    return NAN;
}
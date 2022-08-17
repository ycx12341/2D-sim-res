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
 * If the cell has more than two neighbouring positions which are
 * not occupied, it will proliferate. The original cell will vanish
 * and split into two daughter cells, which will be randomly
 * distributed into two unoccupied neighbouring locations.
 * @param nbr_num number of neighbouring positions
 */
void cell_proliferate(int nbr_num, int nbr_temp[nbr_num], int nbr_coord[nbr_num][2],
                      double ind_position[(int) SPACE_LENGTH_Y][(int) SPACE_LENGTH_X],
                      const int cell_position[2]) {
    int zeros = int_array_count(nbr_num, nbr_temp, 0);
    if (zeros >= 2) {
        int_arr_2_t sample_idx = unif_index2(zeros);

        int *sample_a = nbr_coord[int_array_find(nbr_num, nbr_temp, 0, sample_idx.arr[0])];
        int *sample_b = nbr_coord[int_array_find(nbr_num, nbr_temp, 0, sample_idx.arr[1])];

        ind_position[cell_position[0]][cell_position[1]] = 0;   // Original cell vanishes.
        ind_position[sample_a[0]][sample_a[1]]           = 1;   // Daughter cells being allocated.
        ind_position[sample_b[0]][sample_b[1]]           = 1;

        printf("a: [%d,%d] b: [%d,%d]\n", sample_a[0], sample_a[1], sample_b[0], sample_b[1]);
    }
}

/**
 * Purpose: Generate a SCC invasion pattern with the given parameters.
 * O(n^2)
 * @param pars Parameters
 */
double generate_pattern(sse_pars_t *pars, const int idx) {
    /* Space discretization: Create a 60*35 domain. */
    const int x_len = (int) SPACE_LENGTH_X, y_len = (int) SPACE_LENGTH_Y;
    double    x[x_len], y[y_len];
    assert(seq_length_out(x, 0.0, H * (SPACE_LENGTH_X - 1), x_len) == x_len);
    assert(seq_length_out(y, 0.0, 1, y_len) == y_len);

    /*
     * Initial condition
     */
    // double n0[y_len][x_len], n[y_len][x_len], n0_sort[y_len][x_len];
    // double f0[y_len][x_len], f[y_len][x_len];
    // double m0[y_len][x_len], m[y_len][x_len];

    double n0[y_len][x_len], n0_sort[y_len][x_len], n[y_len][x_len];
    double f0[y_len][x_len];
    double m0[y_len][x_len];

    for (int j = 0; j < x_len; ++j) {
        n0[0][j] = x[j] <= 0.1 ? cos(M_PI * x[j] * 5) : 0;
    }

    MATRIX_ITR(y_len, x_len, {
        n0[_i_][_j_] = n0[0][_j_];
        m0[_i_][_j_] = 0.5 * n0[_i_][_j_];
        f0[_i_][_j_] = 1 - m0[_i_][_j_];
    })

    // MATRIX_MAP(n0, m0, y_len, x_len, n0_to_m0)
    // MATRIX_MAP(n0, f0, y_len, x_len, n0_to_f0)

    MATRIX_COPY(n0, n, y_len, x_len)
    // MATRIX_COPY(f0, f, y_len, x_len)
    // MATRIX_COPY(m0, m, y_len, x_len)

    /* Sort the initial cells */
    MATRIX_COPY(n0, n0_sort, y_len, x_len)
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
        assert(sample.x._int >= 0 && sample.y._int >= 0);

        node_t node;
        node._intPair[0] = sample.x._int;
        node._intPair[1] = sample.y._int;
        arraylist_append(coord, node);
        n0_sort[sample.x._int][sample.y._int] = -DBL_MAX;
    }

    double   ind_position[y_len][x_len];
    for (int i = 0; i < coord->size; ++i) {
        node_t pos = arraylist_get(coord, i);
        ind_position[pos._intPair[0]][pos._intPair[1]] = 1;
    }

    double ind_position_init[y_len][x_len];
    MATRIX_COPY(ind_position, ind_position_init, y_len, x_len)

    /* Set the cut points for the domain (to be used later for density matching & discrepancy calculation.) */
    int    mat_size  = SPACE_LENGTH_Y / 12;
    int    y_cut_len = ceil((SPACE_LENGTH_Y - 1.0) / mat_size);
    int    x_cut_len = ceil((SPACE_LENGTH_X - 1.0) / mat_size);
    double x_cut[x_cut_len], y_cut[y_cut_len];
    assert(seq_by(y_cut, 1, SPACE_LENGTH_Y, mat_size) == y_cut_len);
    assert(seq_by(x_cut, 1, SPACE_LENGTH_X, mat_size) == x_cut_len);

    /* Numerical scheme which solves the PDE system */
    for (int i = 0; i < TIME_STEPS; i++) {

        /* At the end of every day, some current cells in the domain will undergo extinction or mitosis. */
        if ((i + 1) % DAY_TIME_STEPS == 0) {
            /* Extract the density values at the positions of current cells */
            int    cell_den_len = coord->size;
            double cell_den[cell_den_len];

            for (int j = 0; j < coord->size; j++) {
                node_t pos = arraylist_get(coord, j);
                cell_den[j] = n[pos._intPair[0]][pos._intPair[1]];
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
            int    prof_cells_num = (int) round(coord->size * pars->prob_prof[idx]);
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
            MATRIX_INIT(ind_position, y_len, x_len, 0)
            for (int j  = 0; j < coord->size; ++j) {
                node_t pos = arraylist_get(coord, j);
                ind_position[pos._intPair[0]][pos._intPair[1]] = 1;
            }

            /* Proliferation mechanism */
            for (int q = 0; q < prof_cells_num; q++) {
                node_t cell_pos         = arraylist_get(coord, (int) prof_cells[q]);
                int    cell_position[2] = {cell_pos._intPair[0], cell_pos._intPair[1]};

                /* Possible locations for daughter cells (8 surrounding points.) */
                int right_position[2]      = {cell_position[0], cell_position[1] + 1};
                int right_down_position[2] = {cell_position[0] + 1, cell_position[1] + 1};
                int down_position[2]       = {cell_position[0] + 1, cell_position[1]};
                int up_right_position[2]   = {cell_position[0] - 1, cell_position[1] + 1};
                int up_position[2]         = {cell_position[0] - 1, cell_position[1]};
                int up_left_position[2]    = {cell_position[0] - 1, cell_position[1] - 1};
                int left_position[2]       = {cell_position[0], cell_position[1] - 1};
                int left_down_position[2]  = {cell_position[0] + 1, cell_position[1] - 1};

                if (cell_position[0] == 0) {

                    if (cell_position[1] == 0) {
                        /* Special case: top left corner */

                        /* Possible directions to move, check if these points are occupied. */
                        int right      = (int) ind_position[right_position[0]][right_position[1]];
                        int right_down = (int) ind_position[right_down_position[0]][right_down_position[1]];
                        int down       = (int) ind_position[down_position[0]][down_position[1]];

                        /* Neighbouring status*/
                        int neighbouring_temp[3]     = {right, right_down, down};
                        int neighbouring_coord[3][2] = {
                                {right_position[0],      right_position[1]},
                                {right_down_position[0], right_down_position[1]},
                                {down_position[0],       down_position[1]}
                        };

                        cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, cell_position);

                    } else if (cell_position[1] == x_len - 1) {
                        /* Special case: top right corner (same reasoning is followed...) */

                        int left      = (int) ind_position[left_position[0]][left_position[1]];
                        int left_down = (int) ind_position[left_down_position[0]][left_down_position[1]];
                        int down      = (int) ind_position[down_position[0]][down_position[1]];

                        int neighbouring_temp[3]     = {left, left_down, down};
                        int neighbouring_coord[3][2] = {
                                {left_position[0],      left_position[1]},
                                {left_down_position[0], left_down_position[1]},
                                {down_position[0],      down_position[1]}
                        };

                        cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, cell_position);

                    } else {

                        int left       = (int) ind_position[left_position[0]][left_position[1]];
                        int right      = (int) ind_position[right_position[0]][right_position[1]];
                        int down       = (int) ind_position[down_position[0]][down_position[1]];
                        int left_down  = (int) ind_position[left_down_position[0]][left_down_position[1]];
                        int right_down = (int) ind_position[right_down_position[0]][right_down_position[1]];

                        int neighbouring_temp[5]     = {left, right, down, left_down, right_down};
                        int neighbouring_coord[5][2] = {
                                {left_position[0],       left_position[1]},
                                {right_position[0],      right_position[1]},
                                {down_position[0],       down_position[1]},
                                {left_down_position[0],  left_down_position[1]},
                                {right_down_position[0], right_down_position[1]}
                        };

                        cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, cell_position);

                    }
                    /* Lower boundary */
                } else if (cell_position[0] == y_len - 1) {
                    /* Special case: bottom left corner */
                    if (cell_position[1] == 0) {

                    } else if (cell_position[1] == x_len - 1) {

                    } else {

                    }
                } else if (cell_position[1] == 0) {

                } else if (cell_position[1] == x_len - 1) {

                } else {

                }
            }

        }

        /* Solving the PDE model numerically */
        if (i > round(DAY_TIME_STEPS * (95.0 / 96.0))) {

        } else {

        }


    }

// test print
//    for (int k = 0; k < coord->size; ++k) {
//        node_t node  = arraylist_get(coord, k);
//        printf("[%d,%d]", node._intPair[0], node._intPair[1]);
//    }

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
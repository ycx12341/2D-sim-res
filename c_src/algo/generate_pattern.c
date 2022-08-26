#include "algo.h"
#include "../matrix/matrix.h"
#include "../seed/seed.h"

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>

/**
 * Space discretization (60*35) and time discretization.
 */
void discretization(double x[X_LEN], double y[Y_LEN]) {
    assert(seq_length_out(x, 0, H * (X_LEN - 1), X_LEN) == X_LEN);
    assert(seq_length_out(y, 0, 1, Y_LEN) == Y_LEN);
}

/**
 * Initial Condition.
 */
void initial_condition(
        double n[Y_LEN][X_LEN],
        double f[Y_LEN][X_LEN],
        double m[Y_LEN][X_LEN],
        arraylist_t *coord,
        const sse_pars_t *pars,
        const int idx
) {
    assert(coord);
    double x[X_LEN], y[Y_LEN];
    discretization(x, y);

    MATRIX_INIT(n, Y_LEN, X_LEN, 0)
    MATRIX_INIT(f, Y_LEN, X_LEN, 0)
    MATRIX_INIT(m, Y_LEN, X_LEN, 0)

    for (int j = 0; j < X_LEN; ++j) {
        for (int i = 0; i < Y_LEN; ++i) {
            n[i][j] = x[j] <= 0.1 ? cos(M_PI * x[j] * 5) : 0;
            f[i][j] = 1 - 0.5 * n[i][j];
            m[i][j] = 0.5 * n[i][j];
        }
    }

    double n_sort[Y_LEN][X_LEN];
    MATRIX_COPY(n, n_sort, Y_LEN, X_LEN)

    const int   N_CELLS = (int) round(Y_LEN * pars->init_cells_cols[idx]);
    double      max;
    int         count, ind_unique_idx;
    node_t      crd;
    for (pair_t p, ind; coord->size < N_CELLS;) {
        p     = matrix_max(Y_LEN, X_LEN, n_sort);
        max   = p.x._double;
        count = p.y._int;

        ind_unique_idx = count > 1 ? unif_index(count) : 0;
        assert(ind_unique_idx >= 0);
        ind = matrix_find(Y_LEN, X_LEN, n_sort, max, ind_unique_idx + 1);

        crd._intPair[0] = ind.x._int;
        crd._intPair[1] = ind.y._int;
        arraylist_append(coord, crd);

        n_sort[ind.x._int][ind.y._int] = -DBL_MAX;
    }
}

void ind_position_init(
        double ind_position[Y_LEN][X_LEN],
        arraylist_t *coord) {
    MATRIX_INIT(ind_position, Y_LEN, X_LEN, 0)

    node_t   crd;
    for (int i = 0, len = coord->size; i < len; ++i) {
        crd = arraylist_get(coord, i);
        ind_position[crd._intPair[0]][crd._intPair[1]] = 1;
    }
}

bool end_of_day(
        const int t,
        const int idx,
        arraylist_t *coord,
        const sse_pars_t *pars,
        double n[Y_LEN][X_LEN]
) {
    if ((t + 1) % DAY_TIME_STEPS == 0) {
        int    cell_den_len = coord->size;
        double cell_den[cell_den_len];
        ARRAY_INIT(coord->size, cell_den, 0)

        node_t   node;
        for (int i = 0; i < coord->size; ++i) {
            node = arraylist_get(coord, i);
            cell_den[i] = n[node._intPair[0]][node._intPair[1]];
        }

        const int DEAD_CELLS_NUM = (int) round(coord->size * pars->prob_death[idx]);
        int       dead_cells[DEAD_CELLS_NUM];
        ARRAY_INIT(DEAD_CELLS_NUM, dead_cells, 0)

        pair_t   ind_dead;
        double   min;
        for (int i = 0, count, nth, index; i < DEAD_CELLS_NUM; ++i) {
            ind_dead = array_min(coord->size, cell_den);
            min      = ind_dead.x._double;
            count    = ind_dead.y._int;
            nth      = count > 1 ? unif_index(count) : 0;
            index    = double_array_find(coord->size, cell_den, min, nth + 1);
            assert(index >= 0);

            dead_cells[i]   = index;
            cell_den[index] = INFINITY;
        }

        cell_den_len = double_array_delete_many(cell_den_len, cell_den, DEAD_CELLS_NUM, dead_cells);
        arraylist_remove_many(coord, DEAD_CELLS_NUM, dead_cells);
        assert(cell_den_len == coord->size);

        if (cell_den_len == 0) {
            return false;
        }

        const int PROF_CELLS_NUM = (int) round(coord->size * pars->prob_prof[idx]);
        assert(PROF_CELLS_NUM > 0);
        double prof_cells[PROF_CELLS_NUM];
        ARRAY_INIT(PROF_CELLS_NUM, prof_cells, 0)

        pair_t   ind_prof;
        double   max;
        for (int i = 0, count, nth, index; i < PROF_CELLS_NUM; ++i) {
            ind_prof = array_max(coord->size, cell_den);
            max      = ind_prof.x._double;
            count    = ind_prof.y._int;
            nth      = count > 1 ? unif_index(count) : 0;
            index    = double_array_find(coord->size, cell_den, max, nth + 1);
            assert(index >= 0);

            prof_cells[i]   = index;
            cell_den[index] = -INFINITY;
        }

        double ind_position[Y_LEN][X_LEN];
        MATRIX_INIT(ind_position, Y_LEN, X_LEN, 0)

        for (int i = 0; i < coord->size; ++i) {
            node = arraylist_get(coord, i);
            ind_position[node._intPair[0]][node._intPair[1]] = 1;
        }

        /* Proliferation mechanism */
        proliferation(PROF_CELLS_NUM, prof_cells, ind_position, coord);
    }
    return true;
}

bool solve_PDE(
        const int idx,
        arraylist_t *coord,
        const sse_pars_t *pars,
        double n[Y_LEN][X_LEN],
        double f[Y_LEN][X_LEN],
        double m[Y_LEN][X_LEN]
) {
    const int Y_CUT_LEN = (int) ceil((SPACE_LENGTH_Y - 1.0) / MAT_SIZE);
    const int X_CUT_LEN = (int) ceil((SPACE_LENGTH_X - 1.0) / MAT_SIZE);

    double x_cut[X_CUT_LEN], y_cut[Y_CUT_LEN];
    ARRAY_INIT(MAT_SIZE, x_cut, 0)
    ARRAY_INIT(MAT_SIZE, y_cut, 0)

    assert(seq_by(y_cut, 1, SPACE_LENGTH_Y, MAT_SIZE) == Y_CUT_LEN);
    assert(seq_by(x_cut, 1, SPACE_LENGTH_X, MAT_SIZE) == X_CUT_LEN);

    for (int t = 0; t < TIME_STEPS; ++t) {
        if (!end_of_day(t, idx, coord, pars, n)) { return false; }
        pde(t, idx, pars, f, m);
    }
    return true;
}

double generate_pattern(const sse_pars_t *pars, const int idx) {
    double      n[Y_LEN][X_LEN], f[Y_LEN][X_LEN], m[Y_LEN][X_LEN];
    arraylist_t *coord = new_arraylist();
    initial_condition(n, f, m, coord, pars, idx);

    double ind_position[Y_LEN][X_LEN];
    ind_position_init(ind_position, coord);

    if (!solve_PDE(idx, coord, pars, n,f,m)) { return NAN; };

    MATRIX_PRINT(f, Y_LEN, X_LEN, "%.11f ")
    printf("--------------\n");

    arraylist_free(coord);
    return 0;
}
#include "scc.h"

#include <iostream>
#include <cassert>

#include "scc.h"

#define H               parent->h
#define SPACE_LENGTH_X  parent->space_length_x
#define SPACE_LENGTH_Y  parent->space_length_y
#define PARS            parent->pars
#define X_CUT_LEN       parent->x_cut_len
#define Y_CUT_LEN       parent->y_cut_len
#define MAT_SIZE        parent->mat_size

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::Dimension::initial_condition() {
    DBL_T x[X_LEN], y[Y_LEN];
    assert(seq_length_out<DBL_T>(x, 0, H * (SPACE_LENGTH_X - 1), X_LEN) == X_LEN);
    assert(seq_length_out<DBL_T>(y, 0, 1, Y_LEN) == Y_LEN);

    n = new MatrixS<DBL_T, Y_LEN, X_LEN>(0);
    f = new MatrixS<DBL_T, Y_LEN, X_LEN>(0);
    m = new MatrixS<DBL_T, Y_LEN, X_LEN>(0);

    for (int j = 0; j < X_LEN; ++j) {
        (*n)(0, j) = x[j] <= 0.1 ? cos(M_PI * x[j] * 5) : 0;
    }

    for (int i = 1; i < Y_LEN; ++i) {
        for (int j = 0; j < X_LEN; ++j) {
            (*n)(i, j) = (*n)(0, j);
        }
    }

    f->iter_index([&](int i, int j) { (*f)(i, j) = 1 - 0.5 * (*n)(i, j); });
    m->iter_index([&](int i, int j) { (*m)(i, j) = 0.5 * (*n)(i, j); });

    MatrixS<DBL_T, Y_LEN, X_LEN> n_sort(*n);
    const int                    N_CELLS = (int) round((double) SPACE_LENGTH_Y * (double) PARS->INIT_CELLS_COLS[IDX]);
    std::vector<COORD_T >        maxes;

    for (int ind_idx, maxes_size; coord.size() < N_CELLS;) {
        maxes      = n_sort.matrix_which_max();
        maxes_size = (int) maxes.size();
        ind_idx    = maxes_size > 1 ? unif_index(maxes_size) : 0;
        assert(0 <= ind_idx && ind_idx < maxes_size);

        COORD_T ind = maxes.at(ind_idx);
        coord.push_back(ind);
        n_sort(ind[0], ind[1]) = (DBL_T) -INFINITY;
    }

    ind_pos = new MatrixS<DBL_T, Y_LEN, X_LEN>(0);
    for (COORD_T &c: coord) {
        (*ind_pos)(c[0], c[1]) = 1;
    }

    ind_pos_init = new MatrixS<DBL_T, Y_LEN, X_LEN>(*ind_pos);

    x_cut = new DBL_T[X_CUT_LEN];
    y_cut = new DBL_T[Y_CUT_LEN];
    assert(seq_by<DBL_T>(y_cut, 1, SPACE_LENGTH_Y, MAT_SIZE) == Y_CUT_LEN);
    assert(seq_by<DBL_T>(x_cut, 1, SPACE_LENGTH_X, MAT_SIZE) == X_CUT_LEN);
}

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::Dimension::generate_pattern() {
    initial_condition();
    pde();
}

#undef H
#undef SPACE_LENGTH_X
#undef SPACE_LENGTH_Y
#undef PARS
#undef X_CUT_LEN
#undef Y_CUT_LEN
#undef MAT_SIZE

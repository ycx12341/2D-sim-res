#include "scc.h"

#include <iostream>
#include <cassert>

#include "scc.h"

#define H           parent->h
#define SLX         parent->space_length_x
#define SLY         parent->space_length_y
#define PARS        parent->pars
#define X_CUT_LEN   parent->x_cut_len
#define Y_CUT_LEN   parent->y_cut_len
#define MAT_SIZE    parent->mat_size

#define F           (*f)
#define M           (*m)
#define N           (*n)
#define IND_POS     (*ind_pos)

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
void Sim_2D<N_DIMS, Y_LEN, X_LEN>::Dimension::initial_condition() {
    DBL_T x[X_LEN], y[Y_LEN];
    assert(seq_length_out<DBL_T>(x, 0, H * (SLX - 1), X_LEN) == X_LEN);
    assert(seq_length_out<DBL_T>(y, 0, 1, Y_LEN) == Y_LEN);

    n = new MatrixS<DBL_T, Y_LEN, X_LEN>(0);
    f = new MatrixS<DBL_T, Y_LEN, X_LEN>(0);
    m = new MatrixS<DBL_T, Y_LEN, X_LEN>(0);

    for (unsigned j = 0; j < X_LEN; ++j) {
        N(0, j) = x[j] <= 0.1 ? cos(M_PI * x[j] * 5) : 0;
    }

    for (unsigned i = 1; i < Y_LEN; ++i) {
        for (unsigned j = 0; j < X_LEN; ++j) {
            N(i, j) = N(0, j);
        }
    }

    f->iter_index([&](unsigned i, unsigned j) { F(i, j) = 1 - 0.5 * N(i, j); });
    m->iter_index([&](unsigned i, unsigned j) { M(i, j) = 0.5 * N(i, j); });

    MatrixS<DBL_T, Y_LEN, X_LEN> n_sort(N);
    const unsigned               N_CELLS = round((double) SLY * (double) PARS->INIT_CELLS_COLS[IDX]);
    std::vector<COORD_T >        maxes;

    for (unsigned ind_idx, maxes_size; coord.size() < N_CELLS;) {
        maxes      = n_sort.matrix_which_max();
        maxes_size = (unsigned) maxes.size();
        ind_idx    = maxes_size > 1 ? unif_index(maxes_size) : 0;
        assert(0 <= ind_idx && ind_idx < maxes_size);

        COORD_T ind = maxes.at(ind_idx);
        coord.push_back(ind);
        n_sort(ind[0], ind[1]) = (DBL_T) -INFINITY;
    }

    ind_pos = new MatrixS<DBL_T, Y_LEN, X_LEN>(0);
    for (COORD_T &c: coord) {
        IND_POS(c[0], c[1]) = 1;
    }

    ind_pos_init = new MatrixS<DBL_T, Y_LEN, X_LEN>IND_POS;

    x_cut = new DBL_T[X_CUT_LEN];
    y_cut = new DBL_T[Y_CUT_LEN];
    assert(seq_by<DBL_T>(y_cut, 1, SLY, MAT_SIZE) == Y_CUT_LEN);
    assert(seq_by<DBL_T>(x_cut, 1, SLX, MAT_SIZE) == X_CUT_LEN);
}

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
void Sim_2D<N_DIMS, Y_LEN, X_LEN>::Dimension::generate_pattern() {
    initial_condition();
    pde();
}

#undef H
#undef SLX
#undef SLY
#undef PARS
#undef X_CUT_LEN
#undef Y_CUT_LEN
#undef MAT_SIZE
#undef F
#undef M
#undef N
#undef IND_POS

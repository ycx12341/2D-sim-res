#include "calculate_sse.h"

#include <iostream>
#include <cassert>

void Sim_2D::initial_condition() {
    DBL_T x[X_LEN], y[Y_LEN];
    assert(seq_length_out<DBL_T>(x, 0, H * (SPACE_LENGTH_X - 1), X_LEN) == X_LEN);
    assert(seq_length_out<DBL_T>(y, 0, 1, Y_LEN) == Y_LEN);

    n = *new Matrix<DBL_T>(Y_LEN, X_LEN, 0);
    f = *new Matrix<DBL_T>(Y_LEN, X_LEN, 0);
    m = *new Matrix<DBL_T>(Y_LEN, X_LEN, 0);

    for (int j = 0; j < X_LEN; ++j) {
        n(0, j) = x[j] <= 0.1 ? cos(M_PI * x[j] * 5) : 0;
    }

    for (int i = 1; i < Y_LEN; ++i) {
        for (int j = 0; j < X_LEN; ++j) {
            n(i, j) = n(0, j);
        }
    }

    f.iterate_by_index([&](int i, int j) { return 1 - 0.5 * n(i, j); });
    m.iterate_by_index([&](int i, int j) { return 0.5 * n(i, j); });

    Matrix<DBL_T>         n_sort(n);
    const int             N_CELLS = (int) round((double) SPACE_LENGTH_Y * (double) pars->INIT_CELLS_COLS[IDX]);
    std::vector<COORD_T > maxes;

    for (int ind_idx, maxes_size; coord.size() < N_CELLS;) {
        maxes      = n_sort.matrix_which_max();
        maxes_size = (int) maxes.size();
        ind_idx    = maxes_size > 1 ? unif_index(maxes_size) : 0;
        assert(0 <= ind_idx && ind_idx < maxes_size);

        COORD_T ind = maxes.at(ind_idx);
        coord.push_back(ind);
        n_sort(ind[0], ind[1]) = (DBL_T) -INFINITY;
    }

    ind_pos = *new Matrix<DBL_T>(Y_LEN, X_LEN, 0);
    for (COORD_T &c: coord) {
        ind_pos(c[0], c[1]) = 1;
    }

    DBL_T MAT_SIZE         = SPACE_LENGTH_Y / 12.0;
    const int    Y_CUT_LEN = (int) ceil((double) ((SPACE_LENGTH_Y - 1.0) / MAT_SIZE));
    const int    X_CUT_LEN = (int) ceil((double) ((SPACE_LENGTH_X - 1.0) / MAT_SIZE));
    DBL_T x_cut[X_CUT_LEN], y_cut[Y_CUT_LEN];
    assert(seq_by<DBL_T>(y_cut, 1, SPACE_LENGTH_Y, MAT_SIZE) == Y_CUT_LEN);
    assert(seq_by<DBL_T>(x_cut, 1, SPACE_LENGTH_X, MAT_SIZE) == X_CUT_LEN);
}

void *Sim_2D::generate_pattern() {
    initial_condition();
    if (!pde()) { return nullptr; }
}

//for (auto &i: coord) {
//    std::cout << i[0] << " " << i[1] << std::endl;
//}

// std::cout << m << std::endl;

//for (int i = 0; i < Y_CUT_LEN; ++i) {
//std::cout << x_cut[0] << " ";
//}
//std::cout << std::endl;
#include "calculate_sse.h"

#include <iostream>

void Sim_2D::initial_condition() {
    LDBL x[X_LEN], y[Y_LEN];
    assert(seq_length_out<LDBL>(x, 0, H * (SPACE_LENGTH_X - 1), X_LEN) == X_LEN);
    assert(seq_length_out<LDBL>(y, 0, 1, Y_LEN) == Y_LEN);

    n = MATRIX_ZERO(LDBL, Y_LEN, X_LEN);
    f = MATRIX_ONES(LDBL, Y_LEN, X_LEN);
    m = MATRIX_ZERO(LDBL, Y_LEN, X_LEN);

    for (int j = 0; j < X_LEN; ++j) {
        n(0, j) = x[j] <= 0.1 ? cos(M_PI * x[j] * 5) : 0;
    }

    for (int              i = 1; i < Y_LEN; ++i) {
        n.block(i, 0, 1, X_LEN) = n.block(0, 0, 1, X_LEN);
    }

    f = f - 0.5 * n;
    m = 0.5 * n;

    MATRIX_T(LDBL) n_sort         = n;
    const int             N_CELLS = (int) round((double) SPACE_LENGTH_Y * (double) pars->INIT_CELLS_COLS[IDX]);
    std::vector<COORD_T > maxes;

    for (int ind_idx, maxes_size; coord.size() < N_CELLS;) {
        maxes      = matrix_which_max<LDBL>(&n_sort);
        maxes_size = (int) maxes.size();
        ind_idx    = maxes_size > 1 ? unif_index(maxes_size) : 0;

        COORD_T ind = maxes.at(ind_idx);
        coord.push_back(ind);
        n_sort(ind[0], ind[1]) = -LDBL_MAX;
    }

    ind_pos = MATRIX_ZERO(LDBL, Y_LEN, X_LEN);
    for (COORD_T &c: coord) {
        ind_pos(c[0], c[1]) = 1;
    }

    LDBL MAT_SIZE          = SPACE_LENGTH_Y / 12.0;
    const int    Y_CUT_LEN = (int) ceil((double) ((SPACE_LENGTH_Y - 1.0) / MAT_SIZE));
    const int    X_CUT_LEN = (int) ceil((double) ((SPACE_LENGTH_X - 1.0) / MAT_SIZE));
    LDBL x_cut[X_CUT_LEN], y_cut[Y_CUT_LEN];
    assert(seq_by<LDBL>(y_cut, 1, SPACE_LENGTH_Y, MAT_SIZE) == Y_CUT_LEN);
    assert(seq_by<LDBL>(x_cut, 1, SPACE_LENGTH_X, MAT_SIZE) == X_CUT_LEN);
}

void *Sim_2D::generate_pattern() {
    initial_condition();
    if (!solve_pde()) { return nullptr; }
}

//for (auto &i: coord) {
//    std::cout << i[0] << " " << i[1] << std::endl;
//}

// std::cout << m << std::endl;

//for (int i = 0; i < Y_CUT_LEN; ++i) {
//std::cout << x_cut[0] << " ";
//}
//std::cout << std::endl;
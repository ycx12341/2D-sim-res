#include "calculate_sse.h"

#include <iostream>

bool Sim_2D::end_of_day(const int t) {
    if ((t + 1) % DAY_TIME_STEPS == 0) {
        std::vector<LDBL> cell_den;
        for (auto         c: coord) {
            cell_den.push_back(n(c[0], c[1]));
        }

        const int        DEAD_CELLS_NUM = (int) round((double) (coord.size() * pars->PROB_DEATH[IDX]));
        int              dead_cells[DEAD_CELLS_NUM];
        std::vector<int> mins;

        for (int i = 0, mins_size, ind_idx, ind; i < DEAD_CELLS_NUM; ++i) {
            mins      = vector_which_min<LDBL>(cell_den);
            mins_size = (int) mins.size();
            ind_idx   = mins_size > 1 ? unif_index(mins_size) : 0;

            ind = mins.at(ind_idx);
            dead_cells[i] = ind;
            cell_den.at(ind) = (LDBL) INFINITY;
        }

        if ((int) coord.size() == 0) {
            return false;
        }

        coord    = vector_remove_many_by_index<COORD_T >(coord, DEAD_CELLS_NUM, dead_cells);
        cell_den = vector_remove_many_by_index<LDBL>(cell_den, DEAD_CELLS_NUM, dead_cells);
        assert(coord.size() == cell_den.size());

        if ((int) cell_den.size() == 0) { return false; }

        const int        PROF_CELLS_NUM = (int) round((double) (coord.size() * pars->PROB_PROF[IDX]));
        int              prof_cells[PROF_CELLS_NUM];
        std::vector<int> maxes;

        for (int i = 0, maxes_size, ind_idx, ind; i < PROF_CELLS_NUM; ++i) {
            maxes      = vector_which_max<LDBL>(cell_den);
            maxes_size = (int) maxes.size();
            ind_idx    = maxes_size > 1 ? unif_index(maxes_size) : 0;

            ind = maxes.at(ind_idx);
            prof_cells[i] = ind;
            cell_den.at(ind) = (LDBL) -INFINITY;
        }

        ind_pos = MATRIX_ZERO(LDBL, Y_LEN, X_LEN);
        for (COORD_T &i: coord) {
            ind_pos(i[0], i[1]) = 1;
        }
    }

    return true;
}

bool Sim_2D::solve_pde() {
    for (int t = 0; t < TIME_STEPS; ++t) {
        if (!end_of_day(t)) { return false; }
    }
    return true;
}

//for (auto &i: coord) {
//    std::cout << i[0] << " " << i[1] << std::endl;
//}

// std::cout << m << std::endl;

//for (int i = 0; i < Y_CUT_LEN; ++i) {
//std::cout << x_cut[0] << " ";
//}
//std::cout << std::endl;
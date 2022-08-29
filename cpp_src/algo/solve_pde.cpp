#include "calculate_sse.h"

#include <iostream>

template<int Nbr_Num>
void Sim_2D::cell_proliferate(
        std::array<int, Nbr_Num> nbr_temp,
        std::array<COORD_T, Nbr_Num> nghr_cord,
        COORD_T cell_pos
) {
    std::vector<int> zeros = std_array_which_equals<int, Nbr_Num>(nbr_temp, 0);
    if ((int) zeros.size() >= 2) {
        std::array<int, 2> sample = unif_index2((int) zeros.size());
        assert(sample[0] != -1 && sample[1] != -1);

        COORD_T            sample_a = nghr_cord[sample[0]];
        COORD_T            sample_b = nghr_cord[sample[1]];

        ind_pos(cell_pos[0], cell_pos[1]) = 0;
        ind_pos(sample_a[0], sample_a[1]) = 1;
        ind_pos(sample_b[0], sample_b[1]) = 1;
    }
}

void Sim_2D::proliferation(const int PROF_CELLS_NUM, int *prof_cells) {
    for (int i = 0; i < PROF_CELLS_NUM; ++i) {
        COORD_T cell_pos = coord.at(prof_cells[i]);
        const int x      = cell_pos[0];
        const int y      = cell_pos[1];

        COORD_T   _r_pos = {x, y + 1};
        COORD_T   rd_pos = {x + 1, y + 1};
        COORD_T   _d_pos = {x + 1, y};
        COORD_T   ru_pos = {x - 1, y + 1};
        COORD_T   _u_pos = {x - 1, y};
        COORD_T   lu_pos = {x - 1, y - 1};
        COORD_T   _l_pos = {x, y - 1};
        COORD_T   ld_pos = {x + 1, y - 1};

#define R_ (int) ind_pos(_r_pos[0], _r_pos[1])
#define RD (int) ind_pos(rd_pos[0], rd_pos[1])
#define D_ (int) ind_pos(_d_pos[0], _d_pos[1])
#define RU (int) ind_pos(ru_pos[0], ru_pos[1])
#define U_ (int) ind_pos(_u_pos[0], _u_pos[1])
#define LU (int) ind_pos(lu_pos[0], lu_pos[1])
#define L_ (int) ind_pos(_l_pos[0], _l_pos[1])
#define LD (int) ind_pos(ld_pos[0], ld_pos[1])

        if (x == 0) {
            if (y == 0) {
                cell_proliferate<3>({R_, RD, D_}, {_r_pos, rd_pos, _d_pos}, cell_pos);
            } else if (y == X_LEN - 1) {
                cell_proliferate<3>({L_, LD, D_}, {_l_pos, ld_pos, _d_pos}, cell_pos);
            } else {
                cell_proliferate<5>({L_, R_, D_, LD, RD}, {_l_pos, _r_pos, _d_pos, ld_pos, rd_pos}, cell_pos);
            }
        } else if (x == Y_LEN - 1) {
            if (y == 0) {
                cell_proliferate<3>({U_, RU, R_}, {_u_pos, ru_pos, _r_pos}, cell_pos);
            } else if (y == X_LEN - 1) {
                cell_proliferate<3>({L_, U_, LU}, {_l_pos, _u_pos, lu_pos}, cell_pos);
            } else {
                cell_proliferate<5>({L_, LU, U_, RU, R_}, {_l_pos, lu_pos, _u_pos, ru_pos, _r_pos}, cell_pos);
            }
        } else if (y == 0 && x != Y_LEN - 1) {
            cell_proliferate<5>({U_, RU, R_, RD, D_}, {_u_pos, ru_pos, _r_pos, rd_pos, _d_pos}, cell_pos);
        } else if (y == X_LEN - 1 && x != Y_LEN - 1) {
            cell_proliferate<5>({U_, LU, L_, LD, D_}, {_u_pos, lu_pos, _l_pos, ld_pos, _d_pos}, cell_pos);
        } else {
            cell_proliferate<8>({L_, R_, U_, D_, LU, RU, LD, RD},
                                {_l_pos, _r_pos, _u_pos, _d_pos, lu_pos, ru_pos, ld_pos, rd_pos}, cell_pos);
        }

#undef R_
#undef RD
#undef D_
#undef RU
#undef U_
#undef LU
#undef L_
#undef LD
    }
}

bool Sim_2D::end_of_day(const int t) {
    if ((t + 1) % DAY_TIME_STEPS == 0) {
        std::vector<DBL_T> cell_den;
        for (auto          c: coord) {
            cell_den.push_back(n(c[0], c[1]));
        }

        const int        DEAD_CELLS_NUM = (int) round((double) ((DBL_T) coord.size() * pars->PROB_DEATH[IDX]));
        int              dead_cells[DEAD_CELLS_NUM];
        std::vector<int> mins;

        for (int i = 0, mins_size, ind_idx, ind; i < DEAD_CELLS_NUM; ++i) {
            mins      = vector_which_min<DBL_T>(cell_den);
            mins_size = (int) mins.size();
            ind_idx   = mins_size > 1 ? unif_index(mins_size) : 0;

            ind = mins.at(ind_idx);
            dead_cells[i] = ind;
            cell_den.at(ind) = (DBL_T) INFINITY;
        }

        if ((int) coord.size() == 0) {
            return false;
        }

        coord    = vector_remove_many_by_index<COORD_T >(coord, DEAD_CELLS_NUM, dead_cells);
        cell_den = vector_remove_many_by_index<DBL_T>(cell_den, DEAD_CELLS_NUM, dead_cells);
        assert(coord.size() == cell_den.size());

        if ((int) cell_den.size() == 0) { return false; }

        const int        PROF_CELLS_NUM = (int) round((double) ((DBL_T) coord.size() * pars->PROB_PROF[IDX]));
        int              prof_cells[PROF_CELLS_NUM];
        std::vector<int> maxes;

        for (int i = 0, maxes_size, ind_idx, ind; i < PROF_CELLS_NUM; ++i) {
            maxes      = vector_which_max<DBL_T>(cell_den);
            maxes_size = (int) maxes.size();
            ind_idx    = maxes_size > 1 ? unif_index(maxes_size) : 0;

            ind = maxes.at(ind_idx);
            prof_cells[i] = ind;
            cell_den.at(ind) = (DBL_T) -INFINITY;
        }

        ind_pos.setAll(0);
        for (COORD_T &i: coord) {
            ind_pos(i[0], i[1]) = 1;
        }

        proliferation(PROF_CELLS_NUM, prof_cells);

        coord.clear();
        for (COORD_T c: ind_pos.matrix_which_equals(1.0L)) {
            coord.push_back(c);
        }

        ind_pos.setAll(0);
        for (COORD_T c: coord) {
            ind_pos(c[0], c[1]) = 1;
        }
    }

    return true;
}

void Sim_2D::solve_pde(const int t) {
    if ((t + 1) > PDE_TIME_STEPS) {
        f.iterate_range_index(1, 1, Y_LEN - 2, X_LEN - 2, [&](int i, int j) {
            return f(i, j) * (1.0L - DT * pars->ETA[IDX] * m(i, j));
        });

        Matrix<DBL_T> n_cpy(n);
        n.iterate_range_index(1, 1, Y_LEN - 2, X_LEN - 2, [&](int i, int j) {
            return
                    n_cpy(i, j) * (
                            1.0 - (4.0 * DT * pars->DN[IDX] / (H * H)) - (
                                    DT * pars->GAMMA[IDX] / (H * H) * (
                                            f(i, j + 1) + f(i, j - 1) - 4 * f(i, j) + f(i - 1, j) + f(i + 1, j)
                                    )
                            ) + pars->RN[IDX] * (1.0 - n_cpy(i, j) - f(i, j)) * DT
                    )
                    +
                    n_cpy(i, j + 1) * (
                            DT * pars->DN[IDX] / (H * H) - (
                                    DT * pars->GAMMA[IDX] / (4.0 * (H * H)) * (
                                            f(i, j + 1) - f(i, j - 1)
                                    )
                            )
                    )
                    + n_cpy(i, j - 1) * (
                            DT * pars->DN[IDX] / (H * H) + (
                                    DT * pars->GAMMA[IDX] / (4.0 * (H * H)) * (
                                            f(i, j + 1) - f(i, j - 1)
                                    )
                            )
                    )
                    + n_cpy(i - 1, j) * (
                            DT * pars->DN[IDX] / (H * H) - (
                                    DT * pars->GAMMA[IDX] / (4.0 * (H * H)) * (
                                            f(i - 1, j) - f(i + 1, j)
                                    )
                            )
                    )
                    + n_cpy(i + 1, j) * (
                            DT * pars->DN[IDX] / (H * H) + (
                                    DT * pars->GAMMA[IDX] / (4.0 * (H * H)) * (
                                            f(i - 1, j) - f(i + 1, j)
                                    )
                            )
                    )
                    ;
        });
        Matrix<DBL_T> m_cpy(m);
        m.iterate_range_index(1, 1, Y_LEN - 2, X_LEN - 2, [&](int i, int j) {
            return m_cpy(i, j) * (1.0 - (4.0 * DT * pars->DM[IDX] / (H * H))) +
                   DT * pars->ALPHA[IDX] * n(i, j) +
                   DT * pars->DM[IDX] / (H * H) * (
                           m_cpy(i, j + 1) + m_cpy(i, j - 1) + m_cpy(i - 1, j) + m_cpy(i + 1, j)
                   );
        });
    } else {
        f.iterate_range_index(1, 1, Y_LEN - 2, X_LEN - 2, [&](int i, int j) {
            return f(i, j) * (1.0L - DT * pars->ETA[IDX] * m(i, j));
        });
        n.iterate_range_index(1, 1, Y_LEN - 2, X_LEN - 2, [=](int i, int j) {
            return n(i, j) * (1 + pars->RN[IDX] * (1 - n(i, j) - f(i, j)) * DT);
        });
        m.iterate_range_index(1, 1, Y_LEN - 2, X_LEN - 2, [=](int i, int j) {
            return m(i, j) + pars->ALPHA[IDX] * DT * n(i, j);
        });
    }
}

bool Sim_2D::pde() {
    for (int t = 0; t < TIME_STEPS; ++t) {
        if (!end_of_day(t)) { return false; }
        solve_pde(t);
    }
    m.print("%.7f ");
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
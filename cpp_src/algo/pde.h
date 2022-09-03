#include <iostream>

#include "scc.h"

#define H               parent->h
#define DAY_TIME_STEPS  parent->day_time_steps
#define PARS            parent->pars
#define DT              parent->dt
#define PDE_TIME_STEPS  parent->pde_time_steps
#define X_CUT_LEN       parent->x_cut_len
#define Y_CUT_LEN       parent->y_cut_len
#define MAT_SIZE        parent->mat_size
#define TIME_STEPS      parent->time_steps

template<int Y_LEN, int X_LEN>
template<int Nbr_Num>
void Sim_2D<Y_LEN, X_LEN>::Dimension::cell_proliferate(
        const std::array<int, Nbr_Num> &nbr_temp,
        const std::array<COORD_T, Nbr_Num> &nghr_cord,
        COORD_T cell_pos
) {
    std::vector<int> zeros = std_array_which_equals<int, Nbr_Num>(nbr_temp, 0);
    if ((int) zeros.size() >= 2) {
        std::vector<int> sample = unif_index(2, (int) zeros.size());
        assert(sample.size() == 2);

        const int sample_a_idx = zeros.at(sample.at(0));
        const int sample_b_idx = zeros.at(sample.at(1));
        assert(sample_a_idx != sample_b_idx);

        COORD_T          sample_a = nghr_cord[sample_a_idx];
        COORD_T          sample_b = nghr_cord[sample_b_idx];

        (*ind_pos)(cell_pos[0], cell_pos[1]) = 0;
        (*ind_pos)(sample_a[0], sample_a[1]) = 1;
        (*ind_pos)(sample_b[0], sample_b[1]) = 1;
    }
}

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::Dimension::proliferation(const int PROF_CELLS_NUM, int *prof_cells) {
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

#define R_ (int) (*ind_pos)(_r_pos[0], _r_pos[1])
#define RD (int) (*ind_pos)(rd_pos[0], rd_pos[1])
#define D_ (int) (*ind_pos)(_d_pos[0], _d_pos[1])
#define RU (int) (*ind_pos)(ru_pos[0], ru_pos[1])
#define U_ (int) (*ind_pos)(_u_pos[0], _u_pos[1])
#define LU (int) (*ind_pos)(lu_pos[0], lu_pos[1])
#define L_ (int) (*ind_pos)(_l_pos[0], _l_pos[1])
#define LD (int) (*ind_pos)(ld_pos[0], ld_pos[1])

        if (x == 0) {
            if (y == 0) {
                cell_proliferate<3>({R_, RD, D_}, {_r_pos, rd_pos, _d_pos}, cell_pos);
            } else if (y == X_LEN - 1) {
                cell_proliferate<3>({L_, LD, D_}, {_l_pos, ld_pos, _d_pos}, cell_pos);
            } else {
                cell_proliferate<5>({L_, R_, D_, LD, RD}, {_l_pos, _r_pos, _d_pos, ld_pos, rd_pos}, cell_pos);
            }
        } else if (x == Y_LEN - 1) {
            if (y == 0) {   // fixme ORDER DOESN"T MATCH
//                cell_proliferate<3>({U_, RU, R_}, {_u_pos, ru_pos, _r_pos}, cell_pos);
                cell_proliferate<3>({U_, R_, RU}, {_u_pos, _r_pos, ru_pos}, cell_pos);
            } else if (y == X_LEN - 1) {
                cell_proliferate<3>({L_, U_, LU}, {_l_pos, _u_pos, lu_pos}, cell_pos);
            } else {        // fixme ORDER DOESN"T MATCH
                cell_proliferate<5>({L_, LU, U_, RU, R_}, {_l_pos, lu_pos, _u_pos, ru_pos, _r_pos}, cell_pos);
            }
        } else if (y == 0) {
            cell_proliferate<5>({U_, RU, R_, RD, D_}, {_u_pos, ru_pos, _r_pos, rd_pos, _d_pos}, cell_pos);
        } else if (y == X_LEN - 1) {
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

template<int Y_LEN, int X_LEN>
bool Sim_2D<Y_LEN, X_LEN>::Dimension::end_of_day(const int time) {
    if ((time + 1) % DAY_TIME_STEPS == 0) {
        std::vector<DBL_T> cell_den;
        for (COORD_T       c: coord) {
            cell_den.push_back((*n)(c[0], c[1]));
        }

        const int        DEAD_CELLS_NUM = (int) round((double) ((DBL_T) coord.size() * PARS->PROB_DEATH[IDX]));
        int              dead_cells[DEAD_CELLS_NUM];
        std::vector<int> mins;

        for (int i = 0, mins_size, ind_idx, ind; i < DEAD_CELLS_NUM; ++i) {
            mins = vector_which_min<DBL_T>(cell_den);

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

        const int        PROF_CELLS_NUM = (int) round((double) ((DBL_T) coord.size() * PARS->PROB_PROF[IDX]));
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

        (*ind_pos).setAll(0);
        for (COORD_T &i: coord) {
            (*ind_pos)(i[0], i[1]) = 1;
        }

        proliferation(PROF_CELLS_NUM, prof_cells);

        coord.clear();
        for (COORD_T c: (*ind_pos).matrix_which_equals(1)) {
            coord.push_back(c);
        }

        (*ind_pos).setAll(0);
        for (COORD_T c: coord) {
            (*ind_pos)(c[0], c[1]) = 1;
        }
    }
    return true;
}

template<int Y_LEN, int X_LEN>
bool Sim_2D<Y_LEN, X_LEN>::Dimension::solve_pde(const int time) {
    if ((time + 1) > PDE_TIME_STEPS) {
        (*f).iter_range_index(1, 1, Y_LEN - 2, X_LEN - 2, [&](int i, int j) {
            (*f)(i, j) = (*f)(i, j) * (1.0 - DT * PARS->ETA[IDX] * (*m)(i, j));
        });

        MatrixS<DBL_T, Y_LEN, X_LEN> n_cpy(*n);
        (*n).iter_range_index(1, 1, Y_LEN - 2, X_LEN - 2, [&](int i, int j) {
            (*n)(i, j) =
                    n_cpy(i, j) * (
                            1.0 - (4.0 * DT * PARS->DN[IDX] / (H * H)) - (
                                    DT * PARS->GAMMA[IDX] / (H * H) * (
                                            (*f)(i, j + 1) +
                                            (*f)(i, j - 1) -
                                            4 * (*f)(i, j) +
                                            (*f)(i - 1, j) +
                                            (*f)(i + 1, j)
                                    )
                            ) + PARS->RN[IDX] * (1.0 - n_cpy(i, j) - (*f)(i, j)) * DT
                    )
                    +
                    n_cpy(i, j + 1) * (
                            DT * PARS->DN[IDX] / (H * H) - (
                                    DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (
                                            (*f)(i, j + 1) - (*f)(i, j - 1)
                                    )
                            )
                    )
                    + n_cpy(i, j - 1) * (
                            DT * PARS->DN[IDX] / (H * H) + (
                                    DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (
                                            (*f)(i, j + 1) - (*f)(i, j - 1)
                                    )
                            )
                    )
                    + n_cpy(i - 1, j) * (
                            DT * PARS->DN[IDX] / (H * H) - (
                                    DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (
                                            (*f)(i - 1, j) - (*f)(i + 1, j)
                                    )
                            )
                    )
                    + n_cpy(i + 1, j) * (
                            DT * PARS->DN[IDX] / (H * H) + (
                                    DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (
                                            (*f)(i - 1, j) - (*f)(i + 1, j)
                                    )
                            )
                    );
        });
        MatrixS<DBL_T, Y_LEN, X_LEN> m_cpy(*m);
        (*m).iter_range_index(1, 1, Y_LEN - 2, X_LEN - 2, [&](int i, int j) {
            (*m)(i, j) = m_cpy(i, j) * (1.0 - (4.0 * DT * PARS->DM[IDX] / (H * H))) +
                         DT * PARS->ALPHA[IDX] * (*n)(i, j) +
                         DT * PARS->DM[IDX] / (H * H) * (
                                 m_cpy(i, j + 1) + m_cpy(i, j - 1) + m_cpy(i - 1, j) + m_cpy(i + 1, j)
                         );
        });
    } else {
        (*f).iter_range_index(1, 1, Y_LEN - 2, X_LEN - 2, [&](int i, int j) {
            (*f)(i, j) = (*f)(i, j) * (1.0L - DT * PARS->ETA[IDX] * (*m)(i, j));
        });
        (*n).iter_range_index(1, 1, Y_LEN - 2, X_LEN - 2, [&](int i, int j) {
            (*n)(i, j) = (*n)(i, j) * (1 + PARS->RN[IDX] * (1 - (*n)(i, j) - (*f)(i, j)) * DT);
        });
        (*m).iter_range_index(1, 1, Y_LEN - 2, X_LEN - 2, [&](int i, int j) {
            (*m)(i, j) = (*m)(i, j) + PARS->ALPHA[IDX] * DT * (*n)(i, j);
        });
    }

    (*n).iter_cols([&](int j) {
        (*n)(0, j)         = (*n)(1, j);
        (*n)(Y_LEN - 1, j) = (*n)(Y_LEN - 2, j);
    });
    (*n).iter_rows([&](int i) {
        (*n)(i, 0)         = (*n)(i, 1);
        (*n)(i, X_LEN - 1) = (*n)(i, X_LEN - 2);
    });

    (*f).iter_cols([&](int j) {
        (*f)(0, j)         = (*f)(1, j);
        (*f)(Y_LEN - 1, j) = (*f)(Y_LEN - 2, j);
    });
    (*f).iter_rows([&](int i) {
        (*f)(i, 0)         = (*f)(i, 1);
        (*f)(i, X_LEN - 1) = (*f)(i, X_LEN - 2);
    });

    (*m).iter_cols([&](int j) {
        (*m)(0, j)         = (*m)(1, j);
        (*m)(Y_LEN - 1, j) = (*m)(Y_LEN - 2, j);
    });
    (*m).iter_rows([&](int i) {
        (*m)(i, 0)         = (*m)(i, 1);
        (*m)(i, X_LEN - 1) = (*m)(i, X_LEN - 2);
    });

    DBL_T    nv, fv, mv;
    for (int i = 0; i < Y_LEN; ++i) {
        for (int j = 0; j < X_LEN; ++j) {
            nv = (*n)(i, j);
            fv = (*f)(i, j);
            mv = (*m)(i, j);
            if (std::isnan(nv) || nv < 0 ||
                std::isnan(fv) || fv < 0 ||
                std::isnan(mv) || mv < 0) {
                return false;
            }
        }
    }
    return true;
}

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::Dimension::movement(const int time) {
    if ((time + 1) > round(DAY_TIME_STEPS * (95.0 / 96.0))) {
        int x, y;
        DBL_T f_ip1j, f_im1j, f_ijp1, f_ijm1;
        DBL_T p0, p1, p2, p3, p4;

        for (COORD_T &crd: coord) {
            x = crd[0], y = crd[1];

            if (y == 0) {
                if (x == 0) {
                    f_ip1j = (*f)(x, y + 1);
                    f_im1j = 0;
                    f_ijp1 = 0;
                    f_ijm1 = (*f)(x + 1, y);

                    p1 = 0;
                    p2 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ip1j - f_im1j));
                    p3 = (DT * PARS->DN[IDX] / (H * H)) - (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ijp1 - f_ijm1));
                    p4 = 0;
                } else if (x == Y_LEN - 1) {
                    f_ip1j = (*f)(x, y + 1);
                    f_im1j = 0;
                    f_ijp1 = (*f)(x - 1, y);
                    f_ijm1 = 0;

                    p1 = 0;
                    p2 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ip1j - f_im1j));
                    p3 = 0;
                    p4 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ijp1 - f_ijm1));
                } else {
                    f_ip1j = (*f)(x, y + 1);
                    f_im1j = 0;
                    f_ijp1 = (*f)(x - 1, y);
                    f_ijm1 = (*f)(x + 1, y);

                    p1 = 0;
                    p2 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ip1j - f_im1j));
                    p3 = (DT * PARS->DN[IDX] / (H * H)) - (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ijp1 - f_ijm1));
                    p4 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ijp1 - f_ijm1));
                }

                p0 = 1.0 - (4.0 * DT * PARS->DN[IDX] / (H * H)) - (DT * PARS->GAMMA[IDX] / (H * H) * (
                        f_ip1j + f_im1j - 4 * (*f)(x, y) + f_ijp1 + f_ijm1
                )) + PARS->RN[IDX] * (1.0 - (*n)(x, y) - (*f)(x, y)) * DT;
            } else if (y == X_LEN - 1) {
                if (x == 0) {
                    f_ip1j = 0;
                    f_im1j = (*f)(x, y - 1);
                    f_ijp1 = 0;
                    f_ijm1 = (*f)(x + 1, y);

                    p1 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ip1j - f_im1j));
                    p2 = 0;
                    p3 = (DT * PARS->DN[IDX] / (H * H)) - (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ijp1 - f_ijm1));
                    p4 = 0;
                } else if (x == Y_LEN - 1) {
                    f_ip1j = 0;
                    f_im1j = (*f)(x, y - 1);
                    f_ijp1 = (*f)(x - 1, y);
                    f_ijm1 = 0;

                    p1 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ip1j - f_im1j));
                    p2 = 0;
                    p3 = 0;
                    p4 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ijp1 - f_ijm1));
                } else {
                    f_ip1j = 0;
                    f_im1j = (*f)(x, y - 1);
                    f_ijp1 = (*f)(x - 1, y);
                    f_ijm1 = (*f)(x + 1, y);

                    p1 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ip1j - f_im1j));
                    p2 = 0;
                    p3 = (DT * PARS->DN[IDX] / (H * H)) - (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ijp1 - f_ijm1));
                    p4 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ijp1 - f_ijm1));
                }

                p0 = 1.0 - (4.0 * DT * PARS->DN[IDX] / (H * H)) - (DT * PARS->GAMMA[IDX] / (H * H) * (
                        f_ip1j + f_im1j - 4 * (*f)(x, y) + f_ijp1 + f_ijm1
                )) + PARS->RN[IDX] * (1.0 - (*n)(x, y) - (*f)(x, y)) * DT;
            } else if (x == 0) {
                f_ip1j = (*f)(x, y + 1);
                f_im1j = (*f)(x, y - 1);
                f_ijp1 = 0;
                f_ijm1 = (*f)(x + 1, y);

                p0 = 1.0 - (4.0 * DT * PARS->DN[IDX] / (H * H)) - (DT * PARS->GAMMA[IDX] / (H * H) * (
                        f_ip1j + f_im1j - 4 * (*f)(x, y) + f_ijp1 + f_ijm1
                )) + PARS->RN[IDX] * (1.0 - (*n)(x, y) - (*f)(x, y)) * DT;
                p1 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ip1j - f_im1j));
                p2 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ip1j - f_im1j));
                p3 = (DT * PARS->DN[IDX] / (H * H)) - (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ijp1 - f_ijm1));
                p4 = 0;
            } else if (x == Y_LEN - 1) {
                f_ip1j = (*f)(x, y + 1);
                f_im1j = (*f)(x, y - 1);
                f_ijp1 = (*f)(x - 1, y);
                f_ijm1 = 0;

                p0 = 1.0 - (4.0 * DT * PARS->DN[IDX] / (H * H)) - (DT * PARS->GAMMA[IDX] / (H * H) * (
                        f_ip1j + f_im1j - 4 * (*f)(x, y) + f_ijp1 + f_ijm1
                )) + PARS->RN[IDX] * (1.0 - (*n)(x, y) - (*f)(x, y)) * DT;
                p1 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ip1j - f_im1j));
                p2 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ip1j - f_im1j));
                p3 = 0;
                p4 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ijp1 - f_ijm1));
            } else {
                f_ip1j = (*f)(x, y + 1);
                f_im1j = (*f)(x, y - 1);
                f_ijp1 = (*f)(x - 1, y);
                f_ijm1 = (*f)(x + 1, y);

                p0 = 1.0 - (4.0 * DT * PARS->DN[IDX] / (H * H)) - (DT * PARS->GAMMA[IDX] / (H * H) * (
                        f_ip1j + f_im1j - 4 * (*f)(x, y) + f_ijp1 + f_ijm1
                )) + PARS->RN[IDX] * (1.0 - (*n)(x, y) - (*f)(x, y)) * DT;
                p1 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ip1j - f_im1j));
                p2 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ip1j - f_im1j));
                p3 = (DT * PARS->DN[IDX] / (H * H)) - (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ijp1 - f_ijm1));
                p4 = (DT * PARS->DN[IDX] / (H * H)) + (DT * PARS->GAMMA[IDX] / (4.0 * (H * H)) * (f_ijp1 - f_ijm1));
            }

            DBL_T       p[5]  = {p0, p1, p2, p3, p4};
            DBL_T       p_sum = 0;
            for (double &i: p) {
                if (i < 0) { i = 0; } else if (i > 1) { i = 1; }
                p_sum += i;
            }

            const int movement = p_sum == 0 ? 0 : sample_int_index(5, p);
            assert(0 <= movement && movement <= 5);

            (*ind_pos)(x, y) = 0;

            if (movement == 2 && (*ind_pos)(x, y - 1) == 0) {
                y--;
            } else if (movement == 3 && (*ind_pos)(x, y + 1) == 0) {
                y++;
            } else if (movement == 4 && (*ind_pos)(x + 1, y) == 0) {
                x++;
            } else if (movement == 5 && (*ind_pos)(x - 1, y) == 0) {
                x--;
            }

            crd[0] = x;
            crd[1] = y;
            (*ind_pos)(x, y) = 1;
        }
    }
}

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::Dimension::density_matrix(const int time) {
    if ((time + 1) == (DAY_TIME_STEPS * 3)) {
        den_mat_out = new MatrixD<DBL_T>(Y_CUT_LEN, X_CUT_LEN, 0);
        ind_pos_out = new MatrixS<DBL_T, Y_LEN, X_LEN>(*ind_pos);

        for (int i = 0; i < Y_CUT_LEN; ++i) {
            for (int j = 0; j < X_CUT_LEN; ++j) {
                int ones = 0;
                (*ind_pos).iter_range(y_cut[i], x_cut[j], MAT_SIZE - 1, MAT_SIZE - 1, [&](DBL_T val) {
                    if (val == 1) { ones++; }
                });
                (*den_mat_out)(i, j) = ones / (MAT_SIZE * MAT_SIZE);
            }
        }

        n_out = new MatrixS<DBL_T, Y_LEN, X_LEN>(*n);
        f_out = new MatrixS<DBL_T, Y_LEN, X_LEN>(*f);
        m_out = new MatrixS<DBL_T, Y_LEN, X_LEN>(*m);
    }
}

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::Dimension::pde() {
    for (int time = 0; time < TIME_STEPS; ++time) {
        if (!end_of_day(time)) { return; }
        if (!solve_pde(time)) { return; }
        movement(time);
        density_matrix(time);
    }
}

#undef H
#undef DAY_TIME_STEPS
#undef PARS
#undef DT
#undef PDE_TIME_STEPS
#undef X_CUT_LEN
#undef Y_CUT_LEN
#undef MAT_SIZE
#undef TIME_STEPS
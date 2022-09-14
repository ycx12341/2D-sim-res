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

#define F               (*f)
#define M               (*m)
#define N               (*n)
#define IND_POS         (*ind_pos)

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
template<unsigned Nbr_Num>
void Sim_2D<N_DIMS, Y_LEN, X_LEN>::Dimension::cell_proliferate(
        const std::array<int, Nbr_Num> &nbr_temp,
        const std::array<COORD_T, Nbr_Num> &nghr_cord,
        COORD_T cell_pos
) {
    std::vector<unsigned> zeros = std_array_which_equals<int, Nbr_Num>(nbr_temp, 0);
    if ((unsigned) zeros.size() >= 2) {
        std::vector<unsigned> sample      = unif_index(2, (unsigned) zeros.size());
        assert(sample.size() == 2);

        const unsigned sample_a_idx = zeros.at(sample.at(0));
        const unsigned sample_b_idx = zeros.at(sample.at(1));
        assert(sample_a_idx != sample_b_idx);

        COORD_T               sample_a = nghr_cord[sample_a_idx];
        COORD_T               sample_b = nghr_cord[sample_b_idx];

        IND_POS(cell_pos[0], cell_pos[1]) = 0;
        IND_POS(sample_a[0], sample_a[1]) = 1;
        IND_POS(sample_b[0], sample_b[1]) = 1;
    }
}

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
void Sim_2D<N_DIMS, Y_LEN, X_LEN>::Dimension::proliferation(const unsigned PROF_CELLS_NUM, unsigned *prof_cells) {
    for (unsigned i = 0; i < PROF_CELLS_NUM; ++i) {
        COORD_T cell_pos = coord.at(prof_cells[i]);
        const unsigned x = cell_pos[0];
        const unsigned y = cell_pos[1];

        COORD_T        _r_pos = {x, y + 1};
        COORD_T        rd_pos = {x + 1, y + 1};
        COORD_T        _d_pos = {x + 1, y};
        COORD_T        ru_pos = {x - 1, y + 1};
        COORD_T        _u_pos = {x - 1, y};
        COORD_T        lu_pos = {x - 1, y - 1};
        COORD_T        _l_pos = {x, y - 1};
        COORD_T        ld_pos = {x + 1, y - 1};

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
                cell_proliferate<3>({R_, RD, D_},
                                    {_r_pos, rd_pos, _d_pos}, cell_pos);
            } else if (y == X_LEN - 1) {
                cell_proliferate<3>({L_, LD, D_},
                                    {_l_pos, ld_pos, _d_pos}, cell_pos);
            } else {
                cell_proliferate<5>({L_, R_, D_, LD, RD},
                                    {_l_pos, _r_pos, _d_pos, ld_pos, rd_pos}, cell_pos);
            }
        } else if (x == Y_LEN - 1) {
            if (y == 0) {
                cell_proliferate<3>({U_, R_, RU},
                                    {_u_pos, _r_pos, ru_pos}, cell_pos);
            } else if (y == X_LEN - 1) {
                cell_proliferate<3>({L_, U_, LU},
                                    {_l_pos, _u_pos, lu_pos}, cell_pos);
            } else {
                cell_proliferate<5>({L_, LU, U_, RU, R_},
                                    {_l_pos, lu_pos, _u_pos, ru_pos, _r_pos}, cell_pos);
            }
        } else if (y == 0) {
            cell_proliferate<5>({U_, RU, R_, RD, D_},
                                {_u_pos, ru_pos, _r_pos, rd_pos, _d_pos}, cell_pos);
        } else if (y == X_LEN - 1) {
            cell_proliferate<5>({U_, LU, L_, LD, D_},
                                {_u_pos, lu_pos, _l_pos, ld_pos, _d_pos}, cell_pos);
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

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
bool Sim_2D<N_DIMS, Y_LEN, X_LEN>::Dimension::end_of_day(unsigned int time) {
    if ((time + 1) % DAY_TIME_STEPS == 0) {
        std::vector<DBL_T> cell_den;
        for (COORD_T       c: coord) {
            cell_den.push_back(N(c[0], c[1]));
        }

        const unsigned DEAD_CELLS_NUM = round((double) ((DBL_T) coord.size() * PARS->PROB_DEATH[IDX]));
        unsigned       dead_cells[DEAD_CELLS_NUM];

        std::vector<unsigned int> mins;

        for (unsigned i = 0, mins_size, ind_idx, ind; i < DEAD_CELLS_NUM; ++i) {
            mins = vector_which_min<DBL_T>(cell_den);

            mins_size = (int) mins.size();
            ind_idx   = mins_size > 1 ? unif_index(mins_size) : 0;

            ind = mins.at(ind_idx);
            dead_cells[i] = ind;
            cell_den.at(ind) = (DBL_T) INFINITY;
        }

        if ((unsigned) coord.size() == 0) { return false; }

        coord    = vector_remove_many_by_index<COORD_T >(coord, DEAD_CELLS_NUM, dead_cells);
        cell_den = vector_remove_many_by_index<DBL_T>(cell_den, DEAD_CELLS_NUM, dead_cells);
        assert(coord.size() == cell_den.size());

        if ((unsigned) cell_den.size() == 0) { return false; }

        const unsigned PROF_CELLS_NUM = round((double) ((DBL_T) coord.size() * PARS->PROB_PROF[IDX]));
        unsigned       prof_cells[PROF_CELLS_NUM];

        std::vector<unsigned int> maxes;

        for (unsigned i = 0, maxes_size, ind_idx, ind; i < PROF_CELLS_NUM; ++i) {
            maxes      = vector_which_max<DBL_T>(cell_den);
            maxes_size = (unsigned) maxes.size();
            ind_idx    = maxes_size > 1 ? unif_index(maxes_size) : 0;

            ind = maxes.at(ind_idx);
            prof_cells[i] = ind;
            cell_den.at(ind) = (DBL_T) -INFINITY;
        }

        IND_POS.setAll(0);
        for (COORD_T &i: coord) { IND_POS(i[0], i[1]) = 1; }

        proliferation(PROF_CELLS_NUM, prof_cells);

        coord.clear();
        for (COORD_T c: IND_POS.matrix_which_equals(1)) { coord.push_back(c); }

        IND_POS.setAll(0);
        for (COORD_T c: coord) { IND_POS(c[0], c[1]) = 1; }
    }
    return true;
}

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
bool Sim_2D<N_DIMS, Y_LEN, X_LEN>::Dimension::solve_pde(unsigned int time) {
#define MARGIN_X_IDX    1
#define MARGIN_Y_IDX    1
#define MARGIN_X_LEN    (X_LEN - 2)
#define MARGIN_Y_LEN    (Y_LEN - 2)

#define ITER_RANGE(matrix, body) \
    ((matrix).iter_range_index(MARGIN_X_IDX, MARGIN_Y_IDX, MARGIN_Y_LEN, MARGIN_X_LEN, [&](unsigned i, unsigned j) {body} ))

#define NF(i, j)        (F((i), (j) + 1) + F((i), (j) - 1) - 4 * F((i), (j)) + F((i) - 1, (j)) + F((i) + 1, (j)))
#define NF1(i, j)       (F((i), (j) + 1) - F((i), (j) - 1))
#define NF2(i, j)       (F((i) - 1, (j)) - F((i) + 1, (j)))

    const DBL_T   F1 = DT * PARS->ETA[IDX];
    const DBL_T   N0 = DT * PARS->DN[IDX] / (H * H);
    const DBL_T   N1 = 1.0 - (4.0 * N0);
    const DBL_T   N2 = DT * PARS->RN[IDX];
    const DBL_T   N3 = DT * PARS->GAMMA[IDX] / (H * H);
    const DBL_T   N4 = DT * PARS->GAMMA[IDX] / (4.0 * (H * H));
    const DBL_T   M1 = 1.0 - (4.0 * DT * PARS->DM[IDX] / (H * H));
    const DBL_T   M2 = DT * PARS->ALPHA[IDX];
    const DBL_T   M3 = DT * PARS->DM[IDX] / (H * H);

    ITER_RANGE(F, { F(i, j) *= 1.0 - F1 * M(i, j); });

    if ((time + 1) > PDE_TIME_STEPS) {
        MatrixS<DBL_T, Y_LEN, X_LEN> n_cpy(N);
        ITER_RANGE(N, {
            N(i, j) = n_cpy(i, j) * (N1 - N3 * NF(i, j) + N2 * (1.0 - n_cpy(i, j) - F(i, j))) +
                      n_cpy(i, j + 1) * (N0 - N4 * NF1(i, j)) +
                      n_cpy(i, j - 1) * (N0 + N4 * NF1(i, j)) +
                      n_cpy(i - 1, j) * (N0 - N4 * NF2(i, j)) +
                      n_cpy(i + 1, j) * (N0 + N4 * NF2(i, j));
        });
        MatrixS<DBL_T, Y_LEN, X_LEN> m_cpy(M);
        ITER_RANGE(M, {
            M(i, j) = m_cpy(i, j) * M1 + M2 * N(i, j) + M3 * (
                    m_cpy(i, j + 1) + m_cpy(i, j - 1) + m_cpy(i - 1, j) + m_cpy(i + 1, j)
            );
        });
    } else {
        ITER_RANGE(N, { N(i, j) *= (1 + N2 * (1 - N(i, j) - F(i, j))); });
        ITER_RANGE(M, { M(i, j) += M2 * N(i, j); });
    }

    N.iter_cols([&](unsigned j) {
        N(0, j)         = N(1, j);
        F(0, j)         = F(1, j);
        M(0, j)         = M(1, j);
        N(Y_LEN - 1, j) = N(Y_LEN - 2, j);
        F(Y_LEN - 1, j) = F(Y_LEN - 2, j);
        M(Y_LEN - 1, j) = M(Y_LEN - 2, j);
    });
    N.iter_rows([&](unsigned i) {
        N(i, 0)         = N(i, 1);
        F(i, 0)         = F(i, 1);
        M(i, 0)         = M(i, 1);
        N(i, X_LEN - 1) = N(i, X_LEN - 2);
        F(i, X_LEN - 1) = F(i, X_LEN - 2);
        M(i, X_LEN - 1) = M(i, X_LEN - 2);
    });

    DBL_T nv, fv, mv;
    for (unsigned i  = 0; i < Y_LEN; ++i) {
        for (unsigned j = 0; j < X_LEN; ++j) {
            nv = N(i, j);
            fv = F(i, j);
            mv = M(i, j);
            if (std::isnan(nv) || nv < 0 ||
                std::isnan(fv) || fv < 0 ||
                std::isnan(mv) || mv < 0) {
                return false;
            }
        }
    }
    return true;

#undef MARGIN_X_IDX
#undef MARGIN_Y_IDX
#undef MARGIN_X_LEN
#undef MARGIN_Y_LEN
#undef ITER_RANGE
}

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
void Sim_2D<N_DIMS, Y_LEN, X_LEN>::Dimension::movement(unsigned int time) {
    if ((time + 1) > round(PDE_TIME_STEPS)) {
        unsigned x, y;
        DBL_T f_ip1j, f_im1j, f_ijp1, f_ijm1;
        DBL_T p0, p1, p2, p3, p4;

        const DBL_T PF1 = DT * PARS->DN[IDX] / (H * H);
        const DBL_T PF2 = DT * PARS->GAMMA[IDX] / (4.0 * (H * H));
        const DBL_T PF3 = 4.0 * DT * PARS->DN[IDX] / (H * H);
        const DBL_T PF4 = PARS->RN[IDX] * DT;
        const DBL_T PF5 = DT * PARS->GAMMA[IDX] / (H * H);

        for (COORD_T &crd: coord) {
            x = crd[0], y = crd[1];

            f_ip1j = NAN, f_im1j = NAN, f_ijp1 = NAN, f_ijm1 = NAN;
            p1 = NAN, p2 = NAN, p3 = NAN, p4 = NAN;

            if (y == 0) {
                f_im1j = 0;
                p1     = 0;
                if (x == 0) {
                    f_ijp1 = 0;
                    p4     = 0;
                } else if (x == Y_LEN - 1) {
                    f_ijm1 = 0;
                    p3     = 0;
                }
            } else if (y == X_LEN - 1) {
                f_ip1j = 0;
                p2     = 0;
                if (x == 0) {
                    f_ijp1 = 0;
                    p4     = 0;
                } else if (x == Y_LEN - 1) {
                    f_ijm1 = 0;
                    p3     = 0;
                }
            } else if (x == 0) {
                f_ijp1 = 0;
                p4     = 0;
            } else if (x == Y_LEN - 1) {
                f_ijm1 = 0;
                p3     = 0;
            }

            if (std::isnan(f_ip1j)) f_ip1j = F(x, y + 1);
            if (std::isnan(f_im1j)) f_im1j = F(x, y - 1);
            if (std::isnan(f_ijp1)) f_ijp1 = F(x - 1, y);
            if (std::isnan(f_ijm1)) f_ijm1 = F(x + 1, y);

            p0 = 1.0 - PF3 - (PF5 * (
                    f_ip1j + f_im1j - 4 * F(x, y) + f_ijp1 + f_ijm1
            )) + PF4 * (1.0 - N(x, y) - F(x, y));
            if (std::isnan(p1)) p1 = PF1 - PF2 * (f_ip1j - f_im1j);
            if (std::isnan(p2)) p2 = PF1 + PF2 * (f_ip1j - f_im1j);
            if (std::isnan(p3)) p3 = PF1 - PF2 * (f_ijp1 - f_ijm1);
            if (std::isnan(p4)) p4 = PF1 + PF2 * (f_ijp1 - f_ijm1);

            DBL_T      p[5]  = {p0, p1, p2, p3, p4};
            DBL_T      p_sum = 0;
            for (DBL_T &i: p) {
                if (i < 0) { i = 0; } else if (i > 1) { i = 1; }
                p_sum += i;
            }

            const unsigned movement = p_sum == 0 ? 0 : sample_int_index(5, p);
            assert(0 <= movement && movement <= 5);

            IND_POS(x, y) = 0;

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
            IND_POS(x, y) = 1;
        }
    }
}

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
void Sim_2D<N_DIMS, Y_LEN, X_LEN>::Dimension::density_matrix(unsigned int time) {
    if ((time + 1) == (DAY_TIME_STEPS * 3)) {
        den_mat_out = new MatrixD<DBL_T>(Y_CUT_LEN, X_CUT_LEN, 0);
        ind_pos_out = new MatrixS<DBL_T, Y_LEN, X_LEN>(IND_POS);

        for (unsigned i = 0; i < Y_CUT_LEN; ++i) {
            for (unsigned j = 0; j < X_CUT_LEN; ++j) {
                DBL_T sum = 0;
                IND_POS.iter_range(y_cut[i] - 1, x_cut[j] - 1, MAT_SIZE , MAT_SIZE , [&](DBL_T val) {
                    sum += val;
                });
                (*den_mat_out)(i, j) = sum / (MAT_SIZE * MAT_SIZE);
            }
        }

        n_out = new MatrixS<DBL_T, Y_LEN, X_LEN>N;
        f_out = new MatrixS<DBL_T, Y_LEN, X_LEN>F;
        m_out = new MatrixS<DBL_T, Y_LEN, X_LEN>M;
    }
}

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
void Sim_2D<N_DIMS, Y_LEN, X_LEN>::Dimension::pde() {
    for (unsigned time = 0; time < TIME_STEPS; ++time) {
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
#undef F
#undef M
#undef N
#undef IND_POS
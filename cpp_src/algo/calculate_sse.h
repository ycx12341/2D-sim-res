#ifndef CPP_SRC_2D_SIM_ALGO_H
#define CPP_SRC_2D_SIM_ALGO_H

#include <cmath>
#include <cfloat>
#include <vector>

#include "../config.h"
#include "../seed/seed.h"
#include "../collection/collection.h"

class Parameters {
public:
    const int N_DIMS;
    DBL_T *DN;
    DBL_T *GAMMA;
    DBL_T *RN;
    DBL_T *ETA;
    DBL_T *DM;
    DBL_T *ALPHA;
    DBL_T *INIT_CELLS_COLS;
    DBL_T *PROB_DEATH;
    DBL_T *PROB_PROF;

    explicit Parameters(const int n_dims) : N_DIMS(n_dims) {
        DN              = new DBL_T[n_dims];
        GAMMA           = new DBL_T[n_dims];
        RN              = new DBL_T[n_dims];
        ETA             = new DBL_T[n_dims];
        DM              = new DBL_T[n_dims];
        ALPHA           = new DBL_T[n_dims];
        INIT_CELLS_COLS = new DBL_T[n_dims];
        PROB_DEATH      = new DBL_T[n_dims];
        PROB_PROF       = new DBL_T[n_dims];

        runif_seq(DN, N_DIMS, DN_MIN, DN_MAX);
        runif_seq(GAMMA, N_DIMS, GAMMA_MIN, GAMMA_MAX);
        runif_seq(RN, N_DIMS, RN_MIN, RN_MAX);
        runif_seq(ETA, N_DIMS, ETA_MIN, ETA_MAX);
        runif_seq(DM, N_DIMS, DM_MIN, DM_MAX);
        runif_seq(ALPHA, N_DIMS, ALPHA_MIN, ALPHA_MAX);
        runif_seq(INIT_CELLS_COLS, N_DIMS, INIT_CELLS_COLS_MIN, INIT_CELLS_COLS_MAX);
        runif_seq(PROB_DEATH, N_DIMS, PROB_DEATH_MIN, PROB_DEATH_MAX);
        runif_seq(PROB_PROF, N_DIMS, PROB_PROF_MIN, PROB_PROF_MAX);
    }

    ~Parameters() {
        delete DN;
        delete GAMMA;
        delete RN;
        delete ETA;
        delete DM;
        delete ALPHA;
        delete INIT_CELLS_COLS;
        delete PROB_DEATH;
        delete PROB_PROF;
    }
};

template<int Y_LEN, int X_LEN>
class Sim_2D {
public:
    int N_DIMS;

    const DBL_T      H;
    const DBL_T      SPACE_LENGTH_Y;
    const DBL_T      SPACE_LENGTH_X;
    const DBL_T      T;
    const DBL_T      DT;
    const DBL_T      TIME_STEPS;
    const DBL_T      INT_TIME_STEPS;
    const int        DAY_TIME_STEPS;
    const DBL_T      PDE_TIME_STEPS;        // Diffusion starts having an impact after a certain amount of time
    const Parameters *pars;

    Sim_2D(
            const int n_dims,
            DBL_T H,
            DBL_T SPACE_LENGTH_Y,
            DBL_T SPACE_LENGTH_X,
            DBL_T T,
            DBL_T DT,
            DBL_T TIME_STEPS,
            DBL_T INT_TIME_STEPS,
            int DAY_TIME_STEPS,
            DBL_T PDE_TIME_STEPS)
            :
            N_DIMS(n_dims),
            H(H),
            SPACE_LENGTH_Y(SPACE_LENGTH_Y),
            SPACE_LENGTH_X(SPACE_LENGTH_X),
            T(T),
            DT(DT),
            TIME_STEPS(TIME_STEPS),
            INT_TIME_STEPS(INT_TIME_STEPS),
            DAY_TIME_STEPS(DAY_TIME_STEPS),
            PDE_TIME_STEPS(PDE_TIME_STEPS) {
        pars = new Parameters(DEFAULT_N_DIMS);
    }

    ~Sim_2D() {
        pars->~Parameters();
    }

    void calculate_sse(int idx);

private:
    int IDX = 0;

    Matrix<DBL_T, Y_LEN, X_LEN> n;
    Matrix<DBL_T, Y_LEN, X_LEN> f;
    Matrix<DBL_T, Y_LEN, X_LEN> m;
    Matrix<DBL_T, Y_LEN, X_LEN> ind_pos;

    std::vector<COORD_T > coord;

    void initial_condition();

    void *generate_pattern();

    void solve_pde(int t);

    bool pde();

    bool end_of_day(int t);

    void proliferation(int PROF_CELLS_NUM, int *prof_cells);

    template<int Nbr_Num>
    void cell_proliferate(
            std::array<int, Nbr_Num> nbr_temp,
            std::array<COORD_T, Nbr_Num> nghr_cord,
            COORD_T cell_pos
    );
};

class Sim_2D_Factory {
public:
    static auto SCC(const int n_dims) {
        constexpr static const DBL_T H                  = 1.0 / 59.0;
        constexpr static const DBL_T SCC_SPACE_LENGTH_Y = (1.0 / H) + 1.0;
        constexpr static const DBL_T SCC_SPACE_LENGTH_X = (DBL_T) SCC_SPACE_LENGTH_Y * (280.0 / 480.0);
        constexpr static const int   Y_LEN              = (int) SCC_SPACE_LENGTH_Y;
        constexpr static const int   X_LEN              = (int) SCC_SPACE_LENGTH_X;
        constexpr static const DBL_T T                  = 4.52;
        constexpr static const DBL_T DT                 = 0.0025;
        constexpr static const DBL_T TIME_STEPS         = T / DT;
        constexpr static const DBL_T INT_TIME_STEPS     = 1.0 / DT;
        constexpr static const int   DAY_TIME_STEPS     = 600;


        return new Sim_2D<Y_LEN, X_LEN>(
                n_dims, H,
                SCC_SPACE_LENGTH_Y, SCC_SPACE_LENGTH_X,
                T, DT,
                TIME_STEPS, INT_TIME_STEPS, DAY_TIME_STEPS,
                (DBL_T) round((double) DAY_TIME_STEPS * (95.0 / 96.0))
        );
    }
};

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::calculate_sse(const int idx) {
    this->IDX = idx;
    generate_pattern();
}

#endif //CPP_SRC_2D_SIM_ALGO_H

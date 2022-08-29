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

class Sim_2D {
public:
    int              N_DIMS;
    const DBL_T       H;
    const DBL_T       SPACE_LENGTH_Y;
    const DBL_T       SPACE_LENGTH_X;
    const int        X_LEN;
    const int        Y_LEN;
    const DBL_T       T;
    const DBL_T       DT;
    const DBL_T       TIME_STEPS;
    const DBL_T       INT_TIME_STEPS;
    const int        DAY_TIME_STEPS;
    const double     PDE_TIME_STEPS;        // Diffusion starts having an impact after a certain amount of time
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
            double PDE_TIME_STEPS)
            :
            N_DIMS(n_dims),
            H(H),
            SPACE_LENGTH_Y(SPACE_LENGTH_Y),
            SPACE_LENGTH_X(SPACE_LENGTH_X),
            X_LEN((int) SPACE_LENGTH_X),
            Y_LEN((int) SPACE_LENGTH_Y),
            T(T),
            DT(DT),
            TIME_STEPS(TIME_STEPS),
            INT_TIME_STEPS(INT_TIME_STEPS),
            DAY_TIME_STEPS(DAY_TIME_STEPS),
            PDE_TIME_STEPS(PDE_TIME_STEPS) {
        pars = new Parameters(DEFAULT_N_DIMS);
    }

    static Sim_2D sim_scc(const int n_dims) {
        const DBL_T h              = 1.0 / 59.0;
        const DBL_T space_length_y = (1.0 / h) + 1.0;
        const DBL_T t              = 4.52;
        const DBL_T dt             = 0.0025;
        const DBL_T day_time_steps = 600.0;

        return *new Sim_2D(
                n_dims,
                h,
                space_length_y,
                round((double) space_length_y * (280.0 / 480.0)),
                4.52,
                0.0025,
                t / dt,
                1 / dt,
                6000,
                round((double) day_time_steps * (95.0 / 96.0))
        );
    }

    ~Sim_2D() {
        pars->~Parameters();
    }

    void calculate_sse(int idx);

private:
    int IDX;

    Matrix<DBL_T> n;
    Matrix<DBL_T> f;
    Matrix<DBL_T> m;
    Matrix<DBL_T> ind_pos;

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

#endif //CPP_SRC_2D_SIM_ALGO_H

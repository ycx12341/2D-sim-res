#ifndef CPP_SRC_2D_SIM_ALGO_H
#define CPP_SRC_2D_SIM_ALGO_H

#define LDBL    long double

#include <cmath>
#include <cfloat>

#include "../config.h"
#include "../seed/seed.h"
#include "../collection/collection.h"

class Parameters {
public:
    const int N_DIMS;
    LDBL *DN;
    LDBL *GAMMA;
    LDBL *RN;
    LDBL *ETA;
    LDBL *DM;
    LDBL *ALPHA;
    LDBL *INIT_CELLS_COLS;
    LDBL *PROB_DEATH;
    LDBL *PROB_PROF;

    explicit Parameters(const int n_dims) : N_DIMS(n_dims) {
        DN              = new LDBL[n_dims];
        GAMMA           = new LDBL[n_dims];
        RN              = new LDBL[n_dims];
        ETA             = new LDBL[n_dims];
        DM              = new LDBL[n_dims];
        ALPHA           = new LDBL[n_dims];
        INIT_CELLS_COLS = new LDBL[n_dims];
        PROB_DEATH      = new LDBL[n_dims];
        PROB_PROF       = new LDBL[n_dims];

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
    const LDBL       H;
    const LDBL       SPACE_LENGTH_Y;
    const LDBL       SPACE_LENGTH_X;
    const int        X_LEN;
    const int        Y_LEN;
    const LDBL       T;
    const LDBL       DT;
    const LDBL       TIME_STEPS;
    const LDBL       INT_TIME_STEPS;
    const int        DAY_TIME_STEPS;
    const Parameters *pars;

    Sim_2D(
            const int n_dims,
            LDBL H,
            LDBL SPACE_LENGTH_Y,
            LDBL SPACE_LENGTH_X,
            LDBL T,
            LDBL DT,
            LDBL TIME_STEPS,
            LDBL INT_TIME_STEPS,
            int DAY_TIME_STEPS)
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
            DAY_TIME_STEPS(DAY_TIME_STEPS) {
        pars = new Parameters(DEFAULT_N_DIMS);
    }

    static Sim_2D sim_scc(const int n_dims) {
        const LDBL h              = 1.0 / 59.0;
        const LDBL space_length_y = (1.0 / h) + 1.0;
        const LDBL space_length_x = round((double) space_length_y * (280.0 / 480.0));
        const LDBL t              = 4.52;
        const LDBL dt             = 0.0025;
        const LDBL time_steps     = t / dt;
        const LDBL int_time_steps = 1 / dt;
        const LDBL day_time_steps = 600.0;

        return *new Sim_2D(
                n_dims,
                h,
                space_length_y,
                space_length_x,
                t,
                dt,
                time_steps,
                int_time_steps,
                (int) day_time_steps
        );
    }

    ~Sim_2D() {
        pars->~Parameters();
    }

    void calculate_sse(int idx);

private:
    int IDX;

    MATRIX_T(LDBL) n;
    MATRIX_T(LDBL) f;
    MATRIX_T(LDBL) m;
    MATRIX_T(LDBL) ind_pos;

    std::vector<COORD_T > coord;

    void initial_condition();

    void *generate_pattern();

    bool solve_pde();

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

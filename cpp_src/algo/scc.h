#ifndef SIM_2D_CPP_SCC_H
#define SIM_2D_CPP_SCC_H

#include "../config.h"
#include "../seed/seed.h"
#include "../collection/collection.h"

#include <map>

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
private:
    typedef struct {
        DBL_T diff;                             // excluding NAN value
        DBL_T wt;                               // wt_obj
        DBL_T resample;                         // resamp.prob

    } Info_T;
public:
    const int        N_DIMS;
    const DBL_T      h;
    const DBL_T      space_length_y;
    const DBL_T      space_length_x;
    const DBL_T      t;
    const DBL_T      dt;
    const DBL_T      time_steps;
    const DBL_T      int_time_steps;
    const int        day_time_steps;
    const DBL_T      pde_time_steps;            // Diffusion starts having an impact after a certain amount of time
    const DBL_T      mat_size;
    const Parameters *pars;
    int              y_cut_len;
    int              x_cut_len;
    const int        power_len = (int) ceil((POWER_MAX - POWER_MIN) / POWER_STEP);

    std::map<unsigned, DBL_T>  diffs;           // <idx, diff>
    std::map<unsigned, Info_T> infos;           // <idx, { non-NAN diff, ess, ... } >
    std::map<DBL_T, DBL_T>     ess_map;
    DBL_T ess_obj = NAN;
    DBL_T bw_obj  = NAN;

    Sim_2D(const unsigned int seed,
           const int n_dims,
           DBL_T h,
           DBL_T space_length_y,
           DBL_T space_length_x,
           DBL_T t,
           DBL_T dt,
           DBL_T time_steps,
           DBL_T int_time_steps,
           int day_time_steps,
           DBL_T pde_time_steps,
           DBL_T mat_size)
            :
            N_DIMS(n_dims),
            h(h),
            space_length_y(space_length_y),
            space_length_x(space_length_x),
            t(t),
            dt(dt),
            time_steps(time_steps),
            int_time_steps(int_time_steps),
            day_time_steps(day_time_steps),
            pde_time_steps(pde_time_steps),
            mat_size(mat_size) {
        set_seed(seed);
        pars      = new Parameters(DEFAULT_N_DIMS);
        y_cut_len = (int) ceil((double) ((space_length_y - 1.0) / mat_size));
        x_cut_len = (int) ceil((double) ((space_length_x - 1.0) / mat_size));
    }

    ~Sim_2D() {
        pars->~Parameters();
    }

    void simulate();

private:
    void calculate_sse();

    void calculate_bw();

    class Dimension {
    private:
        Sim_2D<Y_LEN, X_LEN> *parent;

        int IDX = 0;
        DBL_T diff = NAN;

        std::vector<COORD_T > coord;

        Matrix<DBL_T> *n            = nullptr;
        Matrix<DBL_T> *f            = nullptr;
        Matrix<DBL_T> *m            = nullptr;
        Matrix<DBL_T> *ind_pos      = nullptr;
        Matrix<DBL_T> *n_out        = nullptr;
        Matrix<DBL_T> *f_out        = nullptr;
        Matrix<DBL_T> *m_out        = nullptr;
        Matrix<DBL_T> *ind_pos_out  = nullptr;
        Matrix<DBL_T> *ind_pos_init = nullptr;
        Matrix<DBL_T> *den_mat_out  = nullptr;

        DBL_T *y_cut = nullptr;
        DBL_T *x_cut = nullptr;

    public:
        explicit Dimension(Sim_2D<Y_LEN, X_LEN> *parent, const int idx) : parent(parent), IDX(idx) {
        };

        ~Dimension() {
            delete[] x_cut;
            delete[] y_cut;
            delete n;
            delete f;
            delete m;
            delete ind_pos;
            delete n_out;
            delete f_out;
            delete m_out;
            delete ind_pos_out;
            delete ind_pos_init;
            delete den_mat_out;
        }

        void calculate();

        DBL_T get_diff() {
            return diff;
        }

    private:
        void initial_condition();

        void generate_pattern();

        void pde();

        bool solve_pde(int time);

        void movement(int time);

        bool end_of_day(int time);

        void proliferation(int PROF_CELLS_NUM, int *prof_cells);

        template<int Nbr_Num>
        void cell_proliferate(
                const std::array<int, Nbr_Num> &nbr_temp,
                const std::array<COORD_T, Nbr_Num> &nghr_cord,
                COORD_T cell_pos
        );

        void density_matrix(int time);
    };
};

class Sim_2D_Factory {
public:
    static auto SCC(const int n_dims, const unsigned int seed) {
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
                seed,
                n_dims, H,
                SCC_SPACE_LENGTH_Y, SCC_SPACE_LENGTH_X,
                T, DT,
                TIME_STEPS, INT_TIME_STEPS, DAY_TIME_STEPS,
                (DBL_T) round((double) DAY_TIME_STEPS * (95.0 / 96.0)),
                SCC_SPACE_LENGTH_Y / 12.0
        );
    }
};

#endif //SIM_2D_CPP_SCC_H

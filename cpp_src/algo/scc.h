#ifndef SIM_2D_CPP_SCC_H
#define SIM_2D_CPP_SCC_H

#include "../config.h"
#include "../seed/seed.h"
#include "../collection/collection.h"

#include <unordered_map>
#include <iomanip>

class Parameters {
public:
    constexpr static const int FEATURES_NUM = 9;
    enum FEATURE_T {
        DN_E, GAMMA_E, RN_E,
        ETA_E, DM_E, ALPHA_E,
        INIT_CELLS_COLS_E, PROB_DEATH_E, PROB_PROF_E
    };

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
        all_zero();
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

    DBL_T &operator()(const FEATURE_T i, const int j) {
        if (i == DN_E) { return DN[j]; }
        if (i == GAMMA_E) { return GAMMA[j]; }
        if (i == RN_E) { return RN[j]; }
        if (i == ETA_E) { return ETA[j]; }
        if (i == DM_E) { return DM[j]; }
        if (i == ALPHA_E) { return ALPHA[j]; }
        if (i == INIT_CELLS_COLS_E) { return INIT_CELLS_COLS[j]; }
        if (i == PROB_DEATH_E) { return PROB_DEATH[j]; }
        if (i == PROB_PROF_E) { return PROB_PROF[j]; }
    }

    friend std::ostream &operator<<(std::ostream &os, Parameters &that) {
        os << "        DN      GAMMA         RN        ETA         DM      ALPHA        ICC         PD         PP"
           << std::endl;
        for (int j = 0; j < that.N_DIMS; ++j) {
            for (int i = 0; i < FEATURES_NUM; ++i) {
                os << std::setw(10) << that((FEATURE_T) i, j) << " ";
            }
            os << std::endl;
        }
        return os;
    }

    void init(bool random_seed = false) const {
        if (random_seed) { set_seed(0); }
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

    Parameters resample(const std::vector<int> &idxes) {
        Parameters p((int) idxes.size());

        for (int i = 0; i < FEATURES_NUM; ++i) {
            for (int j = 0, js = (int) idxes.size(); j < js; ++j) {
                p((FEATURE_T) i, j) = this->operator()((FEATURE_T) i, idxes.at(j));
            }
        }
        return p;
    }

    DBL_T feature_mean(const FEATURE_T i) {
        DBL_T    mean = 0;
        for (int j    = 0; j < N_DIMS; ++j) {
            mean += this->operator()(i, j);
        }
        mean /= N_DIMS;
        return mean;
    }

    DBL_T feature_sd(const FEATURE_T i) {
        DBL_T    var = 0, mean = feature_mean(i);
        for (int j   = 0; j < N_DIMS; ++j) {
            var += pow(this->operator()(i, j) - mean, 2);
        }
        var /= N_DIMS - 1;
        return sqrt(var);
    }

private:
    void all_zero() const {
        for (int i = 0; i < N_DIMS; ++i) {
            DN[i]              = 0;
            GAMMA[i]           = 0;
            RN[i]              = 0;
            ETA[i]             = 0;
            DM[i]              = 0;
            ALPHA[i]           = 0;
            INIT_CELLS_COLS[i] = 0;
            PROB_DEATH[i]      = 0;
            PROB_PROF[i]       = 0;
        }
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
    const int   N_DIMS;
    const DBL_T h;
    const DBL_T space_length_y;
    const DBL_T space_length_x;
    const DBL_T t;
    const DBL_T dt;
    const DBL_T time_steps;
    const DBL_T int_time_steps;
    const int   day_time_steps;
    const DBL_T pde_time_steps;            // Diffusion starts having an impact after a certain amount of time
    const DBL_T mat_size;
    const int   power_len = (int) ceil((POWER_MAX - POWER_MIN) / POWER_STEP);
    int         y_cut_len;
    int         x_cut_len;
    Parameters  *pars;

    std::vector<int>                     nnan_idxs;       // IDXes of which has non-NAN diff
    std::unordered_map<unsigned, DBL_T>  diffs;           // <idx, diff>
    std::unordered_map<unsigned, Info_T> infos;           // <idx, { non-NAN diff, ess, ... } >
    std::unordered_map<DBL_T, DBL_T>     ess_map;         // <power, ess>
    DBL_T ess_obj = NAN;
    DBL_T bw_obj  = NAN;

    const DBL_T ABC_BCD_LB[ABC_BCD_PAR_NUM] = ABC_BCD_PAR_LB;
    const DBL_T ABC_BCD_UB[ABC_BCD_PAR_NUM] = ABC_BCD_PAR_UB;

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
        pars = new Parameters(DEFAULT_N_DIMS);
        pars->init();
        y_cut_len = (int) ceil((double) ((space_length_y - 1.0) / mat_size));
        x_cut_len = (int) ceil((double) ((space_length_x - 1.0) / mat_size));
    }

    ~Sim_2D() {
        pars->~Parameters();
    }

    void simulate();

    Parameters abc_bcd();

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
        constexpr static const int Y_LEN = SCC_Y_LEN;
        constexpr static const int X_LEN = SCC_X_LEN;

        return new Sim_2D<Y_LEN, X_LEN>(
                seed, n_dims, SCC_H,
                SCC_SPACE_LENGTH_Y, SCC_SPACE_LENGTH_X,
                SCC_T, SCC_DT,
                SCC_TIME_STEPS, SCC_INT_TIME_STEPS, SCC_DAY_TIME_STEPS,
                (DBL_T) round((double) SCC_DAY_TIME_STEPS * (95.0 / 96.0)),
                SCC_SPACE_LENGTH_Y / 12.0
        );
    }
};

#endif //SIM_2D_CPP_SCC_H

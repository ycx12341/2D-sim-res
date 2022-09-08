#ifndef SIM_2D_CPP_SCC_H
#define SIM_2D_CPP_SCC_H

#include "pars.h"
#include "../collection/collection.h"

#include <unordered_map>

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
class Sim_2D {
private:
    typedef struct {
        DBL_T least_square;                             // excluding NAN value
        DBL_T wt;                                       // wt.obj
        DBL_T resample;                                 // resamp.prob

    } Info_T;
public:
    const DBL_T        h;
    const DBL_T        space_length_y;
    const DBL_T        space_length_x;
    const DBL_T        t;
    const DBL_T        dt;
    const DBL_T        time_steps;
    const DBL_T        int_time_steps;
    const int          day_time_steps;
    const DBL_T        pde_time_steps;                         // Diffusion starts having an impact after a certain amount of time
    const DBL_T        mat_size;
    const int          power_len = (int) ceil((POWER_MAX - POWER_MIN) / POWER_STEP);
    int                y_cut_len;
    int                x_cut_len;
    Parameters<N_DIMS> pars;

    std::vector<int>                     nnan_idxs;     // IDXes of which has non-NAN least_square
    std::unordered_map<unsigned, DBL_T>  least_square;  // <idx, least_square>
    std::unordered_map<unsigned, Info_T> infos;         // <idx, { non-NAN least_square, ess, ... } >
    std::unordered_map<DBL_T, DBL_T>     ess_map;       // <power, ess>
    DBL_T ess_obj  = NAN;
    DBL_T bw_obj   = NAN;
    DBL_T sum_diff = NAN;

    Sim_2D(const unsigned int seed,
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
        y_cut_len = (int) ceil((double) ((space_length_y - 1.0) / mat_size));
        x_cut_len = (int) ceil((double) ((space_length_x - 1.0) / mat_size));
    }

    Parameters<N_DIMS> simulate(bool multithreading = false);

    void reset() {
        least_square.clear();
        infos.clear();
        nnan_idxs.clear();
        ess_map.clear();
        ess_obj  = NAN;
        bw_obj   = NAN;
        sum_diff = NAN;
    }

    void export_least_square(const std::string &fn = CSV_DIFF_FNAME);

    void export_summary(const std::string &fn = CSV_SMRY_FNAME);

private:
    void calculate_sse(bool multithreading = false);

    void calculate_bw();

    Parameters<N_DIMS> abc_bcd();

    class Dimension {
    private:
        Sim_2D<N_DIMS, Y_LEN, X_LEN> *parent;

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
        explicit Dimension(Sim_2D<N_DIMS, Y_LEN, X_LEN> *parent, const int idx) : parent(parent), IDX(idx) {
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

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
class Sim_2D_Factory {
public:
    static auto SCC(const unsigned int seed) {
        return new Sim_2D<N_DIMS, Y_LEN, X_LEN>(
                seed, SCC_H,
                SCC_SPACE_LENGTH_Y, SCC_SPACE_LENGTH_X,
                SCC_T, SCC_DT,
                SCC_TIME_STEPS, SCC_INT_TIME_STEPS, SCC_DAY_TIME_STEPS,
                (DBL_T) round((double) SCC_DAY_TIME_STEPS * (95.0 / 96.0)),
                SCC_SPACE_LENGTH_Y / 12.0
        );
    }

    static auto SCC_375(const unsigned int seed) {
        return new Sim_2D<N_DIMS, Y_LEN, X_LEN>(
                seed, SCC_H_375,
                SCC_SPACE_LENGTH_Y_375, SCC_SPACE_LENGTH_X_375,
                SCC_T_375, SCC_DT_375,
                SCC_TIME_STEPS_375, SCC_INT_TIME_STEPS_375, SCC_DAY_TIME_STEPS_375,
                (DBL_T) round((double) SCC_DAY_TIME_STEPS_375 * (95.0 / 96.0)),
                SCC_SPACE_LENGTH_Y_375 / 12.0
        );
    }
};

#endif //SIM_2D_CPP_SCC_H

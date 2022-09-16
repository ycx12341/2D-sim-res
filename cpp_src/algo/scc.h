/**
 * 2D simulation.
 */

#ifndef SIM_2D_CPP_SCC_H
#define SIM_2D_CPP_SCC_H

#include "pars.h"
#include "../collection/collection.h"

#include <unordered_map>
#include <filesystem>

/**
 * A 2D simulation task.
 * Simulate cells evolution in X_LEN * Y_LEN grids.
 * @tparam N_DIMS Number of dimensions.
 * @tparam Y_LEN  Number of rows in one grid.
 * @tparam X_LEN  Number of columns in one grid.
 */
template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
class Sim_2D {
private:
    typedef struct {
        DBL_T least_square;                             // least square         (diff)
        DBL_T wt;                                       // weight               (wt.obj)
        DBL_T resample;                                 // resample probability (resamp.prob)
    }     Info_T;

    bool        need_export = false;                    // export csv?
    std::string name;                                   // name of this simulation: for directory name

    unsigned power_len = (unsigned) NAN;                // number of powers
    DBL_T ess_target = NAN;                             // desirable ESS
    DBL_T power_min  = NAN;                             // lower bound of the decimal number sequence to be searched for the correct bandwidth factor.
    DBL_T power_max  = NAN;                             // upper bound of the decimal number sequence to be searched for the correct bandwidth factor.
    DBL_T step_size  = NAN;                             // step of increment for the search of correct bandwidth factor in the sequence.

    std::chrono::time_point<std::chrono::system_clock> start_time;  // when the simulation started
    std::chrono::time_point<std::chrono::system_clock> end_time;    // when the simulation ended

    Matrix<DBL_T> *ref_den = nullptr;                   // cells density
public:
    /* Space discretization */
    const DBL_T h;                                      // 1 / (Y_LEN - 1)
    const DBL_T space_length_y;                         // Space length of grids in Y dimension
    const DBL_T space_length_x;                         // Space length of grids in X dimension

    /* Time discretization */
    const DBL_T    t;
    const DBL_T    dt;                                  // dt chosen as 0.0025 so the stability condition ( dn < ( ( h^2 ) / 4dt ) ) can be maintained
    const unsigned time_steps;                          // Total dimensionless time steps and the dimensionless time steps for one single day.
    const unsigned int_time_steps;                      // initial time steps
    const unsigned day_time_steps;                      // time step
    const unsigned pde_time_steps;                      // diffusion starts having an impact after a certain amount of time

    /* Cut Y_LEN * X_LEN grids into smaller fragments (density matrix). */
    const DBL_T mat_size;                               // size of density matrix
    unsigned    y_cut_len;                              // Y length of density matrix
    unsigned    x_cut_len;                              // X length of density matrix

    Parameters<N_DIMS> *pars;                           // Parameters

    std::vector<unsigned>                nnan_idxs;     // IDXes of which has non-NAN least_square
    std::unordered_map<unsigned, DBL_T>  least_square;  // <idx, least_square>
    std::unordered_map<unsigned, Info_T> infos;         // <idx, { non-NAN least_square, ess, ... } >
    std::unordered_map<DBL_T, DBL_T>     ess_map;       // <power, ess>

    DBL_T ess_obj  = NAN;                               // ess.obj
    DBL_T bw_obj   = NAN;                               // bw.obj
    DBL_T sum_diff = NAN;                               // sum of least squares in all dimensions

    Sim_2D(const unsigned int seed,
           DBL_T h,
           DBL_T space_length_y,
           DBL_T space_length_x,
           DBL_T t,
           DBL_T dt,
           unsigned time_steps,
           unsigned int_time_steps,
           unsigned day_time_steps,
           unsigned pde_time_steps,
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
        y_cut_len = ceil((double) ((space_length_y - 1.0) / mat_size));
        x_cut_len = ceil((double) ((space_length_x - 1.0) / mat_size));
        pars      = new Parameters<N_DIMS>();
    }

    ~Sim_2D() {
        delete pars;
        delete ref_den;
    }

    /**
     * Reset simulation status to initial state.
     */
    void reset() {
        least_square.clear();
        infos.clear();
        nnan_idxs.clear();
        ess_map.clear();

        ess_obj  = NAN;
        bw_obj   = NAN;
        sum_diff = NAN;
    }

    /**
     * Export simulation results to CSV file.
     * @param simulation_name Name of this simulation. (for file/directory name)
     */
    void export_csv(const std::string &simulation_name) {
        std::string dir = CSV_DIR(simulation_name);
        if (!std::filesystem::exists(dir) && !std::filesystem::create_directory(dir)) {
            std::cerr << "[SYSTEM] Unable to create directory for details: " << dir << std::endl;
            return;
        }

        need_export = true;
        name        = simulation_name;
    }

    /**
     * Set conditions for calculating bw that would yield the desirable ESS.
     * @param target desirable ESS.
     * @param lb_bw  lower bound of the decimal number sequence to be searched for the correct bandwidth factor.
     * @param ub_bw  upper bound of the decimal number sequence to be searched for the correct bandwidth factor.
     * @param step   step of increment for the search of correct bandwidth factor in the sequence.
     */
    void set_bw(const DBL_T target, const DBL_T lb_bw, const DBL_T ub_bw, const DBL_T step) {
        if (std::isnan(target) || std::isnan(lb_bw) || std::isnan(ub_bw) || std::isnan(step)) {
            return;
        }
        assert(lb_bw <= ub_bw);
        ess_target = target;
        power_min  = lb_bw;
        power_max  = ub_bw;
        step_size  = step;
        power_len  = (int) ceil((power_max - power_min) / step_size);
    }

    /**
     * Set cells density at a particular day.
     * @param ref_den_matrix Reference Density Matrix.
     */
    void set_ref_den(Matrix<DBL_T> *ref_den_matrix) {
        ref_den = ref_den_matrix;
    }

    /**
     * Simulate and generate parameters for the next round.
     * Basic workflow:
     *      1. Reset simulation state.
     *      2. Generate a SCC invasion pattern with the given parameters.
     *      3. Calculates the least square difference between the simulated invasion pattern and the reference invasion pattern.
     *      4. Calculate the correct bw that would yield the desirable ESS.
     *      5. Given the current set of parameters and the corresponding results of summary statistics, generate a new
     *         set of parameters to be evaluated in the next round.
     * @param multithreading True to use multi threading.
     * @return parameters for the next round.
     */
    Parameters<N_DIMS> *simulate(bool multithreading = false);

    /**
     * @return the number of seconds spent in the simulation.
     */
    long long get_time() const {
        return (std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time)).count();
    }

    /**
     * Calculates the least square difference between the simulated invasion pattern and the reference invasion pattern.
     * @param multithreading True to use multi threading.
     */
    void calculate_sse(bool multithreading = false);

private:
    /**
     * Automatically correct the range of searching the correct bandwidth factor.
     */
    void correct_bw_boundary();

    /**
     * Calculate the correct bandwidth factor which can yield the desirable effective sample size (ESS).
     */
    void calculate_bw();

    /**
     * Given the current set of parameters and the corresponding results of summary statistics, generate a new
     * set of parameters to be evaluated in the next round.
     * @return parameters for the next round.
     */
    Parameters<N_DIMS> *abc_bcd();

    /**
     * Export least squares for all dimensions.
     * @param fn File name.
     */
    void export_least_square(const std::string &fn);

    /**
     * Export summary of this simulation.
     * @param fn File name.
     */
    void export_summary(const std::string &fn);

    /**
     * Dimension.
     * A standard simulation task requires N_DIMS Dimensions.
     */
    class Dimension {
    private:
        Sim_2D<N_DIMS, Y_LEN, X_LEN> *parent;           // Simulation task

        unsigned IDX = 0;                               // Index of the dimension
        DBL_T diff = NAN;                               // Least square

        std::vector<COORD_T > coord;                    // cbind(x.coord, y.corrd)

        Matrix<DBL_T> *n            = nullptr;          // n
        Matrix<DBL_T> *f            = nullptr;          // f
        Matrix<DBL_T> *m            = nullptr;          // m
        Matrix<DBL_T> *ind_pos      = nullptr;          // ind.position
        Matrix<DBL_T> *n_out        = nullptr;          // n for output
        Matrix<DBL_T> *f_out        = nullptr;          // f for output
        Matrix<DBL_T> *m_out        = nullptr;          // m for output
        Matrix<DBL_T> *ind_pos_out  = nullptr;          // ind.position for output
        Matrix<DBL_T> *ind_pos_init = nullptr;          // initial ind.position
        Matrix<DBL_T> *den_mat_out  = nullptr;          // density matrix for output

        DBL_T *y_cut = nullptr;
        DBL_T *x_cut = nullptr;

        std::chrono::time_point<std::chrono::system_clock> start_time;
        std::chrono::time_point<std::chrono::system_clock> end_time;

    public:
        explicit Dimension(Sim_2D<N_DIMS, Y_LEN, X_LEN> *parent, const unsigned int idx) : parent(parent), IDX(idx) {
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

#ifdef EXPORT_CSV

        /**
         * Export the details of a Dimension.
         * @param fn File name.
         */
        void export_details(const std::string &fn) {
            if (std::isnan(get_diff())) { return; }

            std::ofstream csv;

            std::time_t       now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::stringstream fn_s;
            fn_s << std::put_time(std::localtime(&now), fn.c_str());

            csv.open(fn_s.str());

            /* Time */
            if (end_time > start_time) {
                std::chrono::duration period = end_time - start_time;
                csv << "time" << CSV_SEPARATOR
                    << (std::chrono::duration_cast<std::chrono::seconds>(period)).count()
                    << std::endl;
            }

            /* Difference / Least Square */
            csv << "diff" << CSV_SEPARATOR << get_diff() << std::endl << std::endl;

            /* Density Matrix */
            csv << "den.mat";
            den_mat_out->iter_index([&](int i, int j) {
                if (j == 0) { csv << std::endl; }
                csv << std::fixed << std::setprecision(CSV_DBL_PRECISION)
                    << (*den_mat_out)(i, j) << CSV_SEPARATOR;
            });

            /* Ind Positions */
            csv << std::endl << std::endl << "ind.pos";
            ind_pos_out->iter_index([&](int i, int j) {
                if (j == 0) { csv << std::endl; }
                csv << std::fixed << std::setprecision(CSV_DBL_PRECISION)
                    << (*ind_pos_out)(i, j) << CSV_SEPARATOR;
            });

            /* Initial Ind Positions */
            csv << std::endl << std::endl << "ind.pos.init";
            ind_pos_init->iter_index([&](int i, int j) {
                if (j == 0) { csv << std::endl; }
                csv << std::fixed << std::setprecision(CSV_DBL_PRECISION)
                    << (*ind_pos_init)(i, j) << CSV_SEPARATOR;
            });

            /* N */
            csv << std::endl << std::endl << "n";
            n_out->iter_index([&](int i, int j) {
                if (j == 0) { csv << std::endl; }
                csv << std::fixed << std::setprecision(CSV_DBL_PRECISION)
                    << (*n_out)(i, j) << CSV_SEPARATOR;
            });

            /* F */
            csv << std::endl << std::endl << "f";
            f_out->iter_index([&](int i, int j) {
                if (j == 0) { csv << std::endl; }
                csv << std::fixed << std::setprecision(CSV_DBL_PRECISION)
                    << (*f_out)(i, j) << CSV_SEPARATOR;
            });

            /* M */
            csv << std::endl << std::endl << "m";
            m_out->iter_index([&](int i, int j) {
                if (j == 0) { csv << std::endl; }
                csv << std::fixed << std::setprecision(CSV_DBL_PRECISION)
                    << (*m_out)(i, j) << CSV_SEPARATOR;
            });
            csv << std::endl;

            csv.close();
        }

#endif

    private:
        /**
         * Initial condition.
         */
        void initial_condition();

        /**
         * Generate a SCC invasion pattern with the given parameters.
         */
        void generate_pattern();

        /**
         * Numerical scheme which solves the PDE system.
         */
        void pde();

        /**
         * Solve PDE system at a specific time unit.
         * @param time Time in range [0:time_step].
         * @return True if the specific time unit does not yield a NAN value.
         */
        bool solve_pde(unsigned int time);

        /**
         * Roll the dice based on the calculated probabilities, see which direction the cell will choose to move to.
         * Once the direction of movement is confirmed, if the designated position is not occupied, the cell will move to that location.
         * @param time Time in range [0:time_step].
         */
        void movement(unsigned int time);

        /**
         * At the end of everyday, some of the current cells in the domain will undergo extinction or mitosis.
         * @param time Time in range [0:time_step].
         * @return True if the specific time unit does not yield a NAN value.
         */
        bool end_of_day(unsigned int time);

        /**
         * Proliferation mechanism.
         * @param PROF_CELLS_NUM Number of cells in prof_cells.
         * @param prof_cells Cell mitosis: some of the current cells at the locations with the highest densities will undergo mitosis.
         */
        void proliferation(unsigned PROF_CELLS_NUM, unsigned *prof_cells);

        /**
         * If the cell has more than two neighbouring positions which are not occupied, it will proliferate.
         * The original cell will vanish and split into two daughter cells,
         * which will be randomly distributed into two unoccupied neighbouring locations.
         * @tparam Nbr_Num  Number of neighbouring cells.
         * @param nbr_temp  The position of neighbouring cells.
         * @param nghr_cord The value of neighbouring cells.
         * @param cell_pos  The position of current cell.
         */
        template<unsigned Nbr_Num>
        void cell_proliferate(
                const std::array<int, Nbr_Num> &nbr_temp,
                const std::array<COORD_T, Nbr_Num> &nghr_cord,
                COORD_T cell_pos
        );

        /**
         * Compute the density matrix at the end of day.
         * @param time Time in range [0:time_step].
         */
        void density_matrix(unsigned int time);
    };
};

/**
 * Builder for simulation task.
 * @tparam N_DIMS Number of dimensions.
 * @tparam Y_LEN  Number of rows in one grid.
 * @tparam X_LEN  Number of columns in one grid.
 */
template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
class Sim_2D_Builder {
private:
    unsigned           _seed           = (unsigned) NAN;
    DBL_T              _h              = NAN;
    DBL_T              _space_length_y = NAN;
    DBL_T              _space_length_x = NAN;
    DBL_T              _t              = NAN;
    DBL_T              _dt             = NAN;
    unsigned           _time_steps     = (unsigned) NAN;
    unsigned           _int_time_steps = (unsigned) NAN;
    unsigned           _day_time_steps = (unsigned) NAN;
    unsigned           _pde_time_steps = (unsigned) NAN;
    DBL_T              _mat_size       = NAN;
    bool               _need_export    = false;
    std::string        _name;
    Parameters<N_DIMS> *_pars          = nullptr;
    DBL_T              _ess_target     = NAN;
    DBL_T              _power_min      = NAN;
    DBL_T              _power_max      = NAN;
    DBL_T              _step_size      = NAN;
    Matrix<DBL_T>      *_ref_den       = nullptr;
public:
    Sim_2D_Builder seed(unsigned int seed) {
        _seed = seed;
        return *this;
    }

    Sim_2D_Builder h(DBL_T h) {
        _h = h;
        return *this;
    }

    Sim_2D_Builder space_length_y(DBL_T space_length_y) {
        _space_length_y = space_length_y;
        return *this;
    }

    Sim_2D_Builder space_length_x(DBL_T space_length_x) {
        _space_length_x = space_length_x;
        return *this;
    }

    Sim_2D_Builder t(DBL_T t) {
        _t = t;
        return *this;
    }

    Sim_2D_Builder dt(DBL_T dt) {
        _dt = dt;
        return *this;
    }

    Sim_2D_Builder time_step(unsigned time_step) {
        _time_steps = time_step;
        return *this;
    }

    Sim_2D_Builder int_time_steps(unsigned int_time_steps) {
        _int_time_steps = int_time_steps;
        return *this;
    }

    Sim_2D_Builder day_time_steps(unsigned day_time_steps) {
        _day_time_steps = day_time_steps;
        return *this;
    }

    Sim_2D_Builder pde_time_steps(unsigned pde_time_steps) {
        _pde_time_steps = pde_time_steps;
        return *this;
    }

    Sim_2D_Builder mat_size(DBL_T mat_size) {
        _mat_size = mat_size;
        return *this;
    }

    Sim_2D_Builder export_csv(const std::string &simulation_name) {
        _need_export = true;
        _name        = simulation_name;
        return *this;
    }

    Sim_2D_Builder load_pars(Parameters<N_DIMS> *pars = nullptr) {
        _pars = pars;
        return *this;
    }

    Sim_2D_Builder bw(DBL_T ess_target, DBL_T power_min, DBL_T power_max, DBL_T step_size) {
        _ess_target = ess_target;
        _power_min  = power_min;
        _power_max  = power_max;
        _step_size  = step_size;
        return *this;
    }

    Sim_2D_Builder ref_den(Matrix<DBL_T> *ref_den_matrix) {
        _ref_den = ref_den_matrix;
        return *this;
    }

    auto build() {
        auto scc = new Sim_2D<N_DIMS, Y_LEN, X_LEN>(
                _seed, _h, _space_length_y, _space_length_x,
                _t, _dt,
                _time_steps, _int_time_steps, _day_time_steps, _pde_time_steps,
                _mat_size
        );
        if (_need_export) { scc->export_csv(_name); }
        if (_pars == nullptr) { scc->pars->init(); } else { scc->pars = _pars; }
        scc->set_bw(_ess_target, _power_min, _power_max, _step_size);
        scc->set_ref_den(_ref_den);
        return scc;
    }

    static auto SCC_Builder(const unsigned int seed) {
        return Sim_2D_Builder<N_DIMS, Y_LEN, X_LEN>()
                .seed(seed)
                .h(SCC_H)
                .space_length_y(SCC_SPACE_LENGTH_Y)
                .space_length_x(SCC_SPACE_LENGTH_X)
                .t(SCC_T)
                .dt(SCC_DT)
                .time_step(SCC_TIME_STEPS)
                .int_time_steps(SCC_INT_TIME_STEPS)
                .day_time_steps(SCC_DAY_TIME_STEPS)
                .pde_time_steps((DBL_T) round((double) SCC_DAY_TIME_STEPS * (95.0 / 96.0)))
                .mat_size(SCC_SPACE_LENGTH_Y / 12.0);
    }

    static auto SCC_375_Builder(const unsigned int seed) {
        return Sim_2D_Builder<N_DIMS, Y_LEN, X_LEN>()
                .seed(seed)
                .h(SCC_H_375)
                .space_length_y(SCC_SPACE_LENGTH_Y_375)
                .space_length_x(SCC_SPACE_LENGTH_X_375)
                .t(SCC_T_375)
                .dt(SCC_DT_375)
                .time_step(SCC_TIME_STEPS_375)
                .int_time_steps(SCC_INT_TIME_STEPS_375)
                .day_time_steps(SCC_DAY_TIME_STEPS_375)
                .pde_time_steps((DBL_T) round((double) SCC_DAY_TIME_STEPS_375 * (95.0 / 96.0)))
                .mat_size(SCC_SPACE_LENGTH_Y_375 / 12.0);
    }
};

#endif //SIM_2D_CPP_SCC_H

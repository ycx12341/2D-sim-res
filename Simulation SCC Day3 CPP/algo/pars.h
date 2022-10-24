/**
 * Parameters to be evaluated for the simulation of SCC pattern and returns the parameters to be evaluated.
 */

#ifndef SIM_2D_CPP_PARS_H
#define SIM_2D_CPP_PARS_H

#include "../config.h"
#include "../seed/seed.h"

#include <iomanip>
#include <cassert>
#include <iostream>
#include <fstream>

/* All parameters needed */
enum FEATURE_T {
    DN_E, GAMMA_E, RN_E,
    ETA_E, DM_E, ALPHA_E,
    INIT_CELLS_COLS_E, PROB_DEATH_E, PROB_PROF_E
};

/**
 * Parameters.
 * @tparam N_DIMS Number of dimensions.
 */
template<unsigned N_DIMS>
class Parameters {
public:
    constexpr static const DBL_T    NaN          = (DBL_T) NAN;
    constexpr static const unsigned FEATURES_NUM = 9;

    DBL_T                           DN[N_DIMS]              = {0};
    DBL_T                           GAMMA[N_DIMS]           = {0};
    DBL_T                           RN[N_DIMS]              = {0};
    DBL_T                           ETA[N_DIMS]             = {0};
    DBL_T                           DM[N_DIMS]              = {0};
    DBL_T                           ALPHA[N_DIMS]           = {0};
    DBL_T                           INIT_CELLS_COLS[N_DIMS] = {0};
    DBL_T                           PROB_DEATH[N_DIMS]      = {0};
    DBL_T                           PROB_PROF[N_DIMS]       = {0};

    explicit Parameters() = default;

    /**
     * Copy constructor.
     * @param that The other Parameter object.
     */
    Parameters(Parameters<N_DIMS> &that) {
        for (unsigned i = 0; i < FEATURES_NUM; ++i) {
            for (unsigned j = 0; j < N_DIMS; ++j) {
                this->operator()((FEATURE_T) i, j) = that((FEATURE_T) i, j);
            }
        }
    }

    /**
     * Get the reference of a value of a specified parameter type and a specified dimension.
     * @param i Feature/Parameter type.
     * @param j The index of the dimension is being simulated.
     * @return Feature[Dimension]
     */
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
        return const_cast<DBL_T &>(NaN);
    }

    friend std::ostream &operator<<(std::ostream &os, Parameters &that) {
        os << "        DN      GAMMA         RN        ETA         DM      ALPHA        ICC         PD         PP"
           << std::endl;
        for (unsigned j = 0; j < N_DIMS; ++j) {
            for (unsigned i = 0; i < FEATURES_NUM; ++i) {
                os << std::fixed << std::setw(10) << std::setprecision(5)
                   << that((FEATURE_T) i, j) << " ";
            }
            os << std::endl;
        }
        return os;
    }

    /**
     * Initialize parameters.
     * @param random_seed Seed for randomness generator.
     */
    void init(bool random_seed = false) const {
        if (random_seed) { set_seed(0); }
        runif_seq(N_DIMS, (DBL_T *) DN, DN_MIN, DN_MAX);
        runif_seq(N_DIMS, (DBL_T *) GAMMA, GAMMA_MIN, GAMMA_MAX);
        runif_seq(N_DIMS, (DBL_T *) RN, RN_MIN, RN_MAX);
        runif_seq(N_DIMS, (DBL_T *) ETA, ETA_MIN, ETA_MAX);
        runif_seq(N_DIMS, (DBL_T *) DM, DM_MIN, DM_MAX);
        runif_seq(N_DIMS, (DBL_T *) ALPHA, ALPHA_MIN, ALPHA_MAX);
        runif_seq(N_DIMS, (DBL_T *) INIT_CELLS_COLS, INIT_CELLS_COLS_MIN, INIT_CELLS_COLS_MAX);
        runif_seq(N_DIMS, (DBL_T *) PROB_DEATH, PROB_DEATH_MIN, PROB_DEATH_MAX);
        runif_seq(N_DIMS, (DBL_T *) PROB_PROF, PROB_PROF_MIN, PROB_PROF_MAX);
    }

    /**
     * Resample parameters based on the resampling probabilities.
     * @param idxes Resampled indexes of parameters.
     * @param nnan_idxs Indexes of dimensions which have Non-NAN difference.
     * @return Resampled parameters.
     */
    Parameters<N_DIMS> *resample(const std::vector<unsigned> &idxes, const std::vector<unsigned> &nnan_idxs) {
        assert(idxes.size() == N_DIMS);
        auto *p = new Parameters<N_DIMS>;

        for (unsigned i = 0; i < FEATURES_NUM; ++i) {
            for (unsigned j = 0, js = (int) idxes.size(); j < js; ++j) {
                (*p)((FEATURE_T) i, j) = this->operator()((FEATURE_T) i, nnan_idxs.at(idxes.at(j)));
            }
        }
        return p;
    }

    /**
     * Mean value of the specified feature for all dimensions.
     * @param i Feature.
     * @return mean(Feature[:])
     */
    DBL_T feature_mean(const FEATURE_T i) {
        DBL_T    mean = 0;
        for (int j    = 0; j < N_DIMS; ++j) {
            mean += this->operator()(i, j);
        }
        mean /= N_DIMS;
        return mean;
    }

    /**
     * Standard deviation of the specified feature for all dimensions.
     * @param i Feature.
     * @return sd(Feature[:])
     */
    DBL_T feature_sd(const FEATURE_T i) {
        DBL_T    var = 0, mean = feature_mean(i);
        for (int j   = 0; j < N_DIMS; ++j) {
            var += pow(this->operator()(i, j) - mean, 2);
        }
        var /= N_DIMS - 1;
        return sqrt(var);
    }

#ifdef EXPORT_CSV

#define CSV_TITLE_LINE "DN,GAMMA,RN,ETA,DM,ALPHA,ICC,PD,PP"

    /**
     * Export parameters to a csv file.
     * Create file when it does not exist.
     * @param fn File name.
     */
    void export_csv(const std::string &fn) {
        std::ofstream csv;

        std::time_t       now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::stringstream fn_s;
        fn_s << std::put_time(std::localtime(&now), fn.c_str());

        csv.open(fn_s.str());
        std::cout << "[SYSTEM] Exporting Parameters results to " << fn_s.str()
                  << " DIMs: " << N_DIMS << std::endl;
        csv << CSV_TITLE_LINE << std::endl;

        for (unsigned j = 0; j < N_DIMS; ++j) {
            for (unsigned i = 0; i < FEATURES_NUM; ++i) {
                csv << std::fixed << std::setprecision(CSV_DBL_PRECISION)
                    << this->operator()((FEATURE_T) i, j) << CSV_SEPARATOR;
            }
            csv << std::endl;
        }

        csv.close();
    }

    /**
     * Load parameters from a CSV file.
     * @param fn File name.
     */
    void load_csv(const std::string &fn) {
        std::cout << "[SYSTEM] Loading Parameters from " << fn
                  << " DIMs: " << N_DIMS << std::endl;

        std::ifstream csv(fn);
        std::string   line;
        std::string   value;

        std::getline(csv, line);
        assert(line == CSV_TITLE_LINE);

        unsigned int i, j = 0;
        while (std::getline(csv, line)) {
            if (j >= N_DIMS) { break; }
            i = 0;
            std::istringstream line_s(line);
            while (std::getline(line_s, value, CSV_SEPARATOR)) {
                assert(i < PARS_NUM);
                this->operator()((FEATURE_T) i, (int) j) = (DBL_T) std::stold(value);
                i++;
            }
            assert(i == PARS_NUM);
            j++;
        }
        csv.close();

        assert(j == N_DIMS);
    }

#undef CSV_TITLE_LINE

#endif

    /**
     * Mean values of all features for all dimensions.
     * @return 1 dimension Parameters.
     */
    Parameters<1U> features_mean() {
        Parameters<1U> mean;
        for (unsigned j = 0; j < N_DIMS; ++j) {
            for (unsigned i = 0; i < FEATURES_NUM; ++i) {
                mean.operator()((FEATURE_T) i, 0) += this->operator()((FEATURE_T) i, j);
            }
        }
        for (unsigned i = 0; i < FEATURES_NUM; ++i) {
            mean.operator()((FEATURE_T) i, 0) /= N_DIMS;
        }
        return mean;
    }
};

#endif //SIM_2D_CPP_PARS_H

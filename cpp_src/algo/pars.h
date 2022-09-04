#ifndef SIM_2D_CPP_PARS_H
#define SIM_2D_CPP_PARS_H

#include "../config.h"
#include "../seed/seed.h"

#include <iomanip>
#include <cassert>
#include <iostream>
#include <fstream>

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
                os << std::fixed << std::setw(10) << std::setprecision(5)
                   << that((FEATURE_T) i, j) << " ";
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

    void load(Parameters &that) {
        assert(that.N_DIMS == N_DIMS);
        for (int i = 0; i < FEATURES_NUM; ++i) {
            for (int j = 0; j < N_DIMS; ++j) {
                this->operator()((FEATURE_T) i, j) = that((FEATURE_T) i, j);
            }
        }
    }

    Parameters resample(const std::vector<int> &idxes, const std::vector<int> &nnan_idxs) {
        Parameters p((int) idxes.size());

        for (int i = 0; i < FEATURES_NUM; ++i) {
            for (int j = 0, js = (int) idxes.size(); j < js; ++j) {
                p((FEATURE_T) i, j) = this->operator()((FEATURE_T) i, nnan_idxs.at(idxes.at(j)));
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

#ifdef EXPORT_CSV

#define CSV_TITLE_LINE "DN,GAMMA,RN,ETA,DM,ALPHA,ICC,PD,PP"

    void export_csv(const std::string &fn = CSV_PARS_FNAME) {
        std::ofstream csv;

        std::time_t       now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::stringstream fn_s;
        fn_s << std::put_time(std::localtime(&now), fn.c_str());

        csv.open(fn_s.str());
        std::cout << "[SYSTEM] Exporting Parameters results to " << fn_s.str()
                  << " DIMs: " << this->N_DIMS << std::endl;
        csv << CSV_TITLE_LINE << std::endl;

        for (int j = 0; j < N_DIMS; ++j) {
            for (int i = 0; i < FEATURES_NUM; ++i) {
                csv << std::fixed << std::setprecision(CSV_DBL_PRECISION)
                    << this->operator()((FEATURE_T) i, j) << CSV_SEPARATOR;
            }
            csv << std::endl;
        }

        csv.close();
    }

    void load_csv(const std::string &fn) {
        std::cout << "[SYSTEM] Loading Parameters from " << fn
                  << " DIMs: " << this->N_DIMS << std::endl;

        std::ifstream csv(fn);
        std::string   line;
        std::string   value;

        std::getline(csv, line);
        assert(line == CSV_TITLE_LINE);

        unsigned int i, j = 0;
        while (std::getline(csv, line)) {
            assert(j < N_DIMS);
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

#endif //SIM_2D_CPP_PARS_H

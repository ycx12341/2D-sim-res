#ifndef CPP_SRC_2D_SIM_ALGO_H
#define CPP_SRC_2D_SIM_ALGO_H

#include <cmath>
#include <cfloat>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>

#include "scc.h"
#include "../ref_den/t3_ref_den.h"

#define MINIMUM_MULTI_THREAD_DIMS   30
#define NAN_EXIT_CODE               (-1)

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
void Sim_2D<N_DIMS, Y_LEN, X_LEN>::Dimension::calculate() {
    generate_pattern();

    if (n_out == nullptr || f_out == nullptr || m_out == nullptr ||
        den_mat_out == nullptr && ind_pos_out == nullptr) {
        diff = NAN;
    } else {
        DBL_T    sum = 0;
        for (int i   = 0; i < parent->y_cut_len; ++i) {
            for (int j = 0; j < parent->x_cut_len; ++j) {
                sum += pow((*den_mat_out)(i, j) - T3_REF_DEN[i][j], 2);
            }
        }

        diff = sum;
    }
}

DBL_T calculate_ess(const std::vector<DBL_T> &resamp_prob) {
    DBL_T      sum = 0, square_sum = 0;
    for (DBL_T v: resamp_prob) {
        sum += v;
        square_sum += pow(v, 2);
    }
    return pow(sum, 2) / square_sum;
}

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
void Sim_2D<N_DIMS, Y_LEN, X_LEN>::calculate_bw() {
    assert(!infos.empty());
    DBL_T power[power_len];
    assert(seq_by<DBL_T>(power, POWER_MIN, POWER_MAX, POWER_STEP) == power_len);

    std::map<unsigned, DBL_T> wt;
    DBL_T w_min, w_max, ess;

    for (int i = 0; i < power_len; ++i) {
        for (auto &[idx, info]: infos) {
            wt[idx] = pow(info.least_square, -power[i]);
        }

        w_min = map_values_min<unsigned>(wt);
        w_max = map_values_max<unsigned>(wt);
        std::vector<DBL_T> resamp_prob;
        for (auto const &[idx, w]: wt) {
            if (w == w_min) {
                resamp_prob.push_back(0);
            } else if (w == w_max) {
                resamp_prob.push_back(1);
            } else {
                resamp_prob.push_back((w - w_min) / (w_max - w_min));
            }
        }
        assert(resamp_prob.size() == wt.size());

        ess = calculate_ess(resamp_prob);
        if (!std::isnan(ess)) { ess_map.insert({power[i], ess}); }
    }
    assert(ess_map.size() <= power_len);

    DBL_T ess_diff_min = INFINITY, ess_diff;
    for (auto const &[p, e]: ess_map) {
        ess_diff = abs(e - ESS_TARGET);
        if (ess_diff < ess_diff_min) {
            ess_obj      = e;
            ess_diff_min = ess_diff;
            bw_obj       = p;
        }
    }
#ifdef CONSOLE_REPORT
    if (std::isnan(bw_obj) || std::isnan(ess_obj)) {
        std::cerr << "[SYSTEM] No valid ess calculated: Use a different seed or increase the number of dimensions."
                  << std::endl << "[SYSTEM] Program Terminated."
                  << std::endl;
        exit(NAN_EXIT_CODE);
    }
#endif
    assert(!std::isnan(bw_obj) && !std::isnan(ess_obj));

    DBL_T wt_min = INFINITY, wt_max = (DBL_T) -INFINITY, info_wt;
    for (auto &[_, info]: infos) {
        info_wt = pow(info.least_square, -bw_obj);
        info.wt = info_wt;
        if (info_wt < wt_min) { wt_min = info_wt; }
        if (info_wt > wt_max) { wt_max = info_wt; }
    }

    for (auto &[_, info]: infos) {
        info_wt = info.wt;
        if (info_wt == wt_min) {
            info.resample = 0;
        } else if (info_wt == wt_max) {
            info.resample = 1;
        } else {
            info.resample = (info_wt - wt_min) / (wt_max - wt_min);
        }
    }
}

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
Parameters<N_DIMS> Sim_2D<N_DIMS, Y_LEN, X_LEN>::simulate(bool multithreading) {
    reset();
    calculate_sse(multithreading);
    calculate_bw();
    return abc_bcd();
}

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
Parameters<N_DIMS> Sim_2D<N_DIMS, Y_LEN, X_LEN>::abc_bcd() {
    assert(Parameters<N_DIMS>::FEATURES_NUM == PARS_NUM);

    std::vector<DBL_T> probs;
    for (const int     idx: nnan_idxs) {
        probs.push_back(infos[idx].resample);
        assert(!std::isnan(infos[idx].least_square));
    }

    std::vector<int> resamp_idx = sample_indices(N_DIMS, probs, true);
    assert(resamp_idx.size() == N_DIMS);

    Parameters<N_DIMS> paras_nr_unperturbed = pars.resample(resamp_idx, nnan_idxs);
    Parameters<N_DIMS> paras_nr_perturbed;

#define FT(x) (FEATURE_T) (x)
    const DBL_T H            = ABC_BCD_H;
    const DBL_T LB[PARS_NUM] = ABC_BCD_PAR_LB;
    const DBL_T UB[PARS_NUM] = ABC_BCD_PAR_UB;
    DBL_T p;

    for (int i = 0; i < Parameters<N_DIMS>::FEATURES_NUM; ++i) {
        for (int j = 0; j < N_DIMS; ++j) {
            do {
                paras_nr_perturbed(FT(i), j) = rnorm(
                        H * paras_nr_unperturbed(FT(i), j) + (1 - H) * paras_nr_unperturbed.feature_mean(FT(i)),
                        0.05 * paras_nr_unperturbed.feature_sd(FT(i))
                );
                p = paras_nr_perturbed(FT(i), j);
            } while (p > UB[i] || p < LB[i]);
        }
    }
#undef FT
    return paras_nr_perturbed;
}

#ifdef CONSOLE_REPORT

static unsigned FINISHED_TASK = 0;
static auto     TASK_START    = std::chrono::system_clock::now();
static auto     TASK_END      = std::chrono::system_clock::now();

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
void progress_report(const Sim_2D<N_DIMS, Y_LEN, X_LEN> *simulation) {
    FINISHED_TASK++;
    if (FINISHED_TASK == 1) { TASK_START = std::chrono::system_clock::now(); }
    if (FINISHED_TASK %
        (N_DIMS > MINIMUM_MULTI_THREAD_DIMS
         ? (int) round((DBL_T) N_DIMS / MINIMUM_MULTI_THREAD_DIMS)
         : 1) == 0) {
        std::cerr << "-";
    }
    if (FINISHED_TASK == N_DIMS) {
        TASK_END                     = std::chrono::system_clock::now();
        std::chrono::duration period = TASK_END - TASK_START;
        std::cout << std::endl
                  << "Estimated spending time: "
                  << (std::chrono::duration_cast<std::chrono::seconds>(period)).count()
                  << " secs" << std::endl
                  << "Valid results:           " << simulation->infos.size() << std::endl
                  << "Mean Square Differences: " << simulation->sum_diff / simulation->infos.size()
                  << std::endl << std::endl;
        FINISHED_TASK = 0;
    }
}

#else

void progress_report() {};

#endif

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
void Sim_2D<N_DIMS, Y_LEN, X_LEN>::calculate_sse(bool multithreading) {
    sum_diff = 0;

#ifdef CONSOLE_REPORT
    FINISHED_TASK = 0;
#endif

    if (multithreading) {
        const unsigned int suggestion = std::thread::hardware_concurrency();
        if (suggestion <= 1 || N_DIMS <= MINIMUM_MULTI_THREAD_DIMS) {
            std::cerr << "[SYSTEM] MULTITHREADING DISABLED: Hardware concurrency is not well defined or not needed"
                      << std::endl;
            goto single_thread;
        }

        std::mutex               lock;
        std::vector<std::thread> threads;

        const int threads_num = N_DIMS > suggestion ? suggestion : N_DIMS;
        const int batch_size  = ceil((double) N_DIMS / suggestion);

        for (int batch = 0; batch < threads_num; ++batch) {
            threads.push_back(std::thread([&lock, this, batch, batch_size]() {
                for (int d = batch * batch_size; d < (batch + 1) * batch_size && d < N_DIMS; d++) {
                    Dimension dimension(this, d);
                    dimension.calculate();
                    DBL_T df = dimension.get_diff();

                    lock.lock();
                    least_square.insert({d, df});
                    if (!std::isnan(df)) {
                        infos.insert({d, {df, NAN, NAN}});
                        nnan_idxs.push_back(d);
                        sum_diff += df;
                    }
                    lock.unlock();

                    progress_report(this);
                }
            }));
        }

        std::cout << "[SYSTEM] MULTITHREADING ENABLED" << std::endl;
        std::cout << "[SYSTEM] TASK REGISTERED: " << threads_num << " threads in use" << std::endl;
        for (std::thread &thread: threads) {
            thread.join();
        }
    } else {
        single_thread:
        DBL_T    df;
        for (int i = 0; i < N_DIMS; ++i) {
            Dimension dimension(this, i);
            dimension.calculate();

            df = dimension.get_diff();
            least_square.insert({i, df});
            if (!std::isnan(df)) {
                infos.insert({i, {df, NAN, NAN}});
                nnan_idxs.push_back(i);
                sum_diff += df;
            }

            progress_report(this);
        }
    }

    if (infos.empty()) {
        std::cerr << "[SYSTEM] No valid dimension found: Use a different seed or increase the number of dimensions."
                  << std::endl << "[SYSTEM] Program Terminated."
                  << std::endl;
        exit(NAN_EXIT_CODE);
    }
}

#ifdef EXPORT_CSV

#define CSV_TITLE_LINE "Index,LeastSquare"

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
void Sim_2D<N_DIMS, Y_LEN, X_LEN>::export_least_square(const std::string &fn) {
    std::ofstream csv;

    std::time_t       now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::stringstream fn_s;
    fn_s << std::put_time(std::localtime(&now), fn.c_str());

    csv.open(fn_s.str());
    std::cout << "[SYSTEM] Exporting Least Square results to " << fn_s.str()
              << " DIMs: " << this->N_DIMS << std::endl;
    csv << CSV_TITLE_LINE << std::endl;

    for (const auto &[idx, ls]: least_square) {
        csv << idx << CSV_SEPARATOR
            << std::fixed << std::setprecision(CSV_DBL_PRECISION)
            << ls << std::endl;
    }

    csv.close();
}

#undef CSV_TITLE_LINE

#endif

#ifdef EXPORT_CSV

template<unsigned N_DIMS, unsigned Y_LEN, unsigned X_LEN>
void Sim_2D<N_DIMS, Y_LEN, X_LEN>::export_summary(const std::string &fn) {
    std::ofstream csv;

    std::time_t       now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::stringstream fn_s;
    fn_s << std::put_time(std::localtime(&now), fn.c_str());

    csv.open(fn_s.str());
    std::cout << "[SYSTEM] Exporting Summary to " << fn_s.str()
              << " DIMs: " << this->N_DIMS << std::endl;

    csv << std::put_time(std::localtime(&now), CSV_TIME_FORMAT)
        << std::endl << std::endl;

    csv << "non-NAN results" << CSV_SEPARATOR << nnan_idxs.size() << std::endl;
    csv << "ess.obj" << CSV_SEPARATOR << ess_obj << std::endl;
    csv << "bw.obj" << CSV_SEPARATOR << bw_obj << std::endl;

    csv << std::endl
        << "Valid Index" << CSV_SEPARATOR
        << "Least Square" << CSV_SEPARATOR
        << "wt.obj" << CSV_SEPARATOR
        << "resamp.prob.obj" << std::endl;
    for (const auto &[idx, info]: infos) {
        csv << idx << CSV_SEPARATOR
            << std::fixed << std::setprecision(CSV_DBL_PRECISION)
            << info.least_square << CSV_SEPARATOR
            << info.wt << CSV_SEPARATOR
            << info.resample << std::endl;
    }

    csv << std::endl
        << "power" << CSV_DBL_PRECISION
        << "ess" << std::endl;
    for (const auto &[p, ess]: ess_map) {
        csv << std::fixed << std::setprecision(CSV_DBL_PRECISION)
            << p << CSV_SEPARATOR
            << ess << std::endl;
    }

    csv.close();
}

#undef CSV_TITLE_LINE

#endif

#endif //CPP_SRC_2D_SIM_ALGO_H

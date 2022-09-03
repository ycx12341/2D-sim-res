#ifndef CPP_SRC_2D_SIM_ALGO_H
#define CPP_SRC_2D_SIM_ALGO_H

#include <cmath>
#include <cfloat>
#include <vector>

#include "scc.h"

static const DBL_T T3_REF_DEN[12][7] = {
        {0.26500, 0.00000, 0.00000, 0, 0, 0, 0.00000},
        {0.35125, 0.00562, 0.00500, 0, 0, 0, 0.00000},
        {0.46812, 0.00000, 0.00000, 0, 0, 0, 0.00000},
        {0.48000, 0.00000, 0.00000, 0, 0, 0, 0.00000},
        {0.47938, 0.00000, 0.00000, 0, 0, 0, 0.00375},
        {0.48812, 0.00938, 0.00000, 0, 0, 0, 0.00000},
        {0.51750, 0.00000, 0.00000, 0, 0, 0, 0.00000},
        {0.38750, 0.00000, 0.00000, 0, 0, 0, 0.00000},
        {0.41812, 0.00000, 0.00000, 0, 0, 0, 0.00000},
        {0.44875, 0.00000, 0.00125, 0, 0, 0, 0.00000},
        {0.31875, 0.00000, 0.00875, 0, 0, 0, 0.00000},
        {0.44625, 0.00000, 0.00000, 0, 0, 0, 0.00000},
};

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::Dimension::calculate() {
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
        diff         = sum;
    }
}

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::calculate_sse() {
    DBL_T diff;
    for (int i = 0; i < N_DIMS; ++i) {
        Dimension dimension(this, i);
        dimension.calculate();

        diff = dimension.get_diff();
        diffs.insert({i, diff});
        if (!std::isnan(diff)) {
            infos.insert({i, {diff, NAN, NAN}});
            nnan_idxs.push_back(i);
        }
        std::cout << i << " -> " << diff << std::endl;  // TODO
    }

    // TODO output diffs and print mean (?)

    DBL_T    mean = 0;
    for (auto const &[_, info]: infos) { mean += info.diff; }
    mean /= infos.size();
    std::cout << "Mean: " << mean << std::endl;
}

DBL_T calculate_ess(const std::vector<DBL_T> &resamp_prob) {
    DBL_T      sum = 0, square_sum = 0;
    for (DBL_T v: resamp_prob) {
        sum += v;
        square_sum += pow(v, 2);
    }
    return pow(sum, 2) / square_sum;
}

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::calculate_bw() {
    assert(!infos.empty());
    DBL_T power[power_len];
    assert(seq_by<DBL_T>(power, POWER_MIN, POWER_MAX, POWER_STEP) == power_len);

    std::map<unsigned, DBL_T> wt;

    for (int i = 0; i < power_len; ++i) {
        for (auto          &[idx, info]: infos) {
            wt[idx] = pow(info.diff, -power[i]);
        }

        DBL_T w_min = map_values_min<unsigned>(wt);
        DBL_T w_max = map_values_max<unsigned>(wt);
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

        DBL_T ess = calculate_ess(resamp_prob);
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
    assert(!std::isnan(bw_obj) && !std::isnan(ess_obj));

    DBL_T wt_min = INFINITY, wt_max = (DBL_T) -INFINITY, info_wt;
    for (auto &[_, info]: infos) {
        info_wt = pow(info.diff, -bw_obj);
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

template<int Y_LEN, int X_LEN>
Parameters Sim_2D<Y_LEN, X_LEN>::simulate() {
    calculate_sse();
    calculate_bw();
    return abc_bcd();
}

template<int Y_LEN, int X_LEN>
Parameters Sim_2D<Y_LEN, X_LEN>::abc_bcd() {
    assert(Parameters::FEATURES_NUM == ABC_BCD_PAR_NUM);

    std::vector<DBL_T> probs;
    for (const int     idx: nnan_idxs) {
        probs.push_back(infos[idx].resample);
        assert(!std::isnan(infos[idx].diff));
    }

    std::vector<int> resamp_idx = sample_indices(N_DIMS, probs, true);
    assert(resamp_idx.size() == N_DIMS);

    Parameters paras_nr_unperturbed = pars->resample(resamp_idx, nnan_idxs);
    Parameters paras_nr_perturbed(N_DIMS);

#define FT(x) (Parameters::FEATURE_T) x
    const DBL_T H                   = ABC_BCD_H;
    const DBL_T LB[ABC_BCD_PAR_NUM] = ABC_BCD_PAR_LB;
    const DBL_T UB[ABC_BCD_PAR_NUM] = ABC_BCD_PAR_UB;
    DBL_T p;

    for (int i = 0; i < Parameters::FEATURES_NUM; ++i) {
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

#endif //CPP_SRC_2D_SIM_ALGO_H

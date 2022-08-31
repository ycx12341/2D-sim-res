#ifndef CPP_SRC_2D_SIM_ALGO_H
#define CPP_SRC_2D_SIM_ALGO_H

#include <cmath>
#include <cfloat>
#include <vector>

#include "scc.h"

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::Dimension::calculate() {
    generate_pattern();

    if (n_out == nullptr || f_out == nullptr || m_out == nullptr ||
        den_mat_out == nullptr && ind_pos_out == nullptr) {
        diff = NAN;
    } else {
        diff = den_mat_out->sum<DBL_T>();
    }
}

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::calculate_sse() {
    DBL_T diff;
    std::pair<unsigned, DBL_T> pair;

    for (int i = 0; i < N_DIMS; ++i) {
        std::cout << i << std::endl;    // TODO

        Dimension dimension(this, i);
        dimension.calculate();

        diff = dimension.get_diff();
        pair = {i, diff};
        diffs.insert(pair);
        if (!std::isnan(diff)) { diffs_valid.insert(pair); }
    }

    // TODO output diffs and print mean (?)

    DBL_T                      mean = 0;
    for (auto const &[_, d]: diffs_valid) { mean += d; }
    mean /= diffs_valid.size();
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
    DBL_T power[power_len];
    assert(seq_by<DBL_T>(power, POWER_MIN, POWER_MAX, POWER_STEP) == power_len);

    std::vector<DBL_T>        ess_vec;
    std::map<unsigned, DBL_T> wt(diffs_valid);

    for (int i = 0; i < power_len; ++i) {
        for (auto const    &[idx, d]: diffs_valid) {
            wt[idx] = pow(d, -power[i]);
        }

        DBL_T min = map_values_min<unsigned>(wt);
        DBL_T max = map_values_max<unsigned>(wt);
        std::vector<DBL_T> resamp_prob;
        for (auto const &[idx, d]: wt) {
            if (d == min) {
                resamp_prob.push_back(0);
            } else if (d == max) {
                resamp_prob.push_back(1);
            } else {
                resamp_prob.push_back((d - min) / (max - min));
            }
        }
        assert(resamp_prob.size() == wt.size());

        ess_vec.push_back(calculate_ess(resamp_prob));
    }
    assert(ess_vec.size() == power_len);
}

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::simulate() {
    calculate_sse();
    calculate_bw();
}

#endif //CPP_SRC_2D_SIM_ALGO_H
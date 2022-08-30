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
        if (!std::isnan(diff)) { diffs_nNAN.insert(pair); }
    }

    // TODO output diffs and print mean (?)

    DBL_T                      mean = 0;
    for (auto const &[_, d]: diffs_nNAN) { mean += d; }
    mean /= diffs_nNAN.size();
}

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::calculate_bw() {

}

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::simulate() {
    calculate_sse();
    calculate_bw();
}

#endif //CPP_SRC_2D_SIM_ALGO_H

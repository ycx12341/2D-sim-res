#ifndef CPP_SRC_2D_SIM_ALGO_H
#define CPP_SRC_2D_SIM_ALGO_H

#include <cmath>
#include <cfloat>
#include <vector>

#include "scc.h"

template<int Y_LEN, int X_LEN>
void Sim_2D<Y_LEN, X_LEN>::calculate_sse(const int idx) {
    this->IDX = idx;
    generate_pattern();

//    assert(n_out != nullptr && f_out != nullptr && m_out != nullptr &&
//           dsy_mat_out != nullptr && ind_pos_out != nullptr);
}

#endif //CPP_SRC_2D_SIM_ALGO_H

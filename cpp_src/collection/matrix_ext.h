#ifndef SIM_2D_CPP_MATRIX_EXT_H
#define SIM_2D_CPP_MATRIX_EXT_H

#include <Eigen/Eigen>

#define MATRIX_T(type)          Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic>

#define MATRIX_ZERO(type, r, c) Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic>::Zero(r, c)

#define MATRIX_ONES(type, r, c) Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic>::Ones(r, c)

#define COORD_T                 std::array<int, 2>

template<typename T>
std::vector<COORD_T > matrix_which_max(const MATRIX_T(T) *matrix) {
    std::vector<COORD_T > maxes;

    T max = -INFINITY;
    T v;

    if (matrix->size() <= 0) { return maxes; }

    for (int i = 0, r = (int) matrix->rows(); i < r; ++i) {
        for (int j = 0, c = (int) matrix->cols(); j < c; ++j) {
            v = (*matrix)(i, j);
            if (v < max) { continue; }
            if (v > max) {
                maxes.clear();
                max = v;
            }
            maxes.push_back({i, j});
        }
    }
    return maxes;
}

#endif //SIM_2D_CPP_MATRIX_EXT_H

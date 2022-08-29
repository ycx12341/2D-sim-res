#ifndef SIM_2D_CPP_MATRIX_H
#define SIM_2D_CPP_MATRIX_H

#include <cassert>
#include <array>
#include <cmath>
#include <vector>
#include <iostream>

#define COORD_T std::array<int, 2>

template<typename T, int ROWS, int COLS>
class Matrix {

private:
    T MATRIX[ROWS][COLS];

public:
    int Rows = 0;
    int Cols = 0;

    Matrix() = default;

    explicit Matrix(const T val) {
        setAll(val);
    }

    Matrix(Matrix<T, ROWS, COLS> &that) {
        for (int i = 0; i < ROWS; ++i) {
            for (int j = 0; j < COLS; ++j) {
                MATRIX[i][j] = that.MATRIX[i][j];
            }
        }
    }

    T &operator()(const int i, const int j) {
        assert(0 <= i && i < ROWS && 0 <= j && j < COLS);
        return MATRIX[i][j];
    }

    void setAll(T val) {
        for (int i = 0; i < ROWS; ++i) {
            for (int j = 0; j < COLS; ++j) {
                MATRIX[i][j] = val;
            }
        }
    }

    template<typename F>
    void iter(F f) {
        iter_by_index(0, 0, ROWS, COLS, f);
    }

    template<typename F>
    void iter_range(const int i, const int j, const int r, const int c, F f) {
        iter_range_index(0, 0, ROWS, COLS, f);
    }

    template<typename F>
    void iter_by_index(F f) {
        for (int i = 0; i < ROWS; ++i) {
            for (int j = 0; j < COLS; ++j) {
                f(i, j);
            }
        }
    }

    template<typename F>
    void iter_range_index(const int i, const int j, const int r, const int c, F f) {
        for (int ii = i, il = i + r; ii < il; ++ii) {
            for (int jj = j, jl = j + c; jj < jl; ++jj) {
                f(ii, jj);
            }
        }
    }

    template<typename F>
    void iter_cols(F f) {
        for (int j = 0; j < COLS; ++j) {
            f(j);
        }
    }

    template<typename F>
    void iter_rows(F f) {
        for (int i = 0; i < ROWS; ++i) {
            f(i);
        }
    }

    long long size() {
        return ROWS * COLS;
    }

    void print(const char *format) {
        for (int i = 0; i < ROWS; ++i) {
            for (int j = 0; j < COLS; ++j) {
                printf(format, MATRIX[i][j]);
            }
            printf("\n");
        }
    }

    std::vector<COORD_T > matrix_which_max() {
        std::vector<COORD_T > maxes;

        T max = -INFINITY;
        T v;

        if (size() <= 0) { return maxes; }

        for (int i = 0; i < ROWS; ++i) {
            for (int j = 0; j < COLS; ++j) {
                v = MATRIX[i][j];
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

    std::vector<COORD_T > matrix_which_equals(const T val) {
        std::vector<COORD_T > res;
        if (size() <= 0) { return res; }

        for (int j = 0; j < COLS; ++j) {
            for (int i = 0; i < ROWS; ++i) {
                if (MATRIX[i][j] == val) { res.push_back({i, j}); }
            }
        }
        return res;
    }
};

#undef MATRIX

#endif //SIM_2D_CPP_MATRIX_H

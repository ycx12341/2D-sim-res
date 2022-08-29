#ifndef SIM_2D_CPP_MATRIX_H
#define SIM_2D_CPP_MATRIX_H

#include <cassert>
#include <array>
#include <cmath>
#include <vector>
#include <iostream>

#define COORD_T std::array<int, 2>

#define MATRIX(i, j) MATRIX[i * Cols + j]

template<typename T>
class Matrix {

private:
    T *MATRIX = nullptr;

public:
    int Rows = 0;
    int Cols = 0;

    Matrix() = default;

    Matrix(const int rows, const int cols) : Rows(rows), Cols(cols) {
        MATRIX = new T[rows * cols];
    }

    Matrix(const int rows, const int cols, const T val) : Matrix(rows, cols) {
        setAll(val);
    }

    Matrix(Matrix<T> &that) {
        Rows       = that.Rows;
        Cols       = that.Cols;
        MATRIX     = new T[that.Rows * that.Cols];
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                MATRIX(i, j) = that.MATRIX(i, j);
            }
        }
    }

    ~Matrix() {
        delete[] MATRIX;
    }

    T &operator()(const int i, const int j) {
        assert(0 <= i && i < Rows && 0 <= j && j < Cols && MATRIX != nullptr);
        return MATRIX(i, j);
    }

    void setAll(T val) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                MATRIX(i, j) = val;
            }
        }
    }

    template<typename F>
    void iterate(F f) {
        iterate_by_index(0, 0, Rows, Cols, f);
    }

    template<typename F>
    void iterate_range(const int i, const int j, const int r, const int c, F f) {
        iterate_range_index(0, 0, Rows, Cols, f);
    }

    template<typename F>
    void iterate_by_index(F f) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                MATRIX(i, j) = f(i, j);
            }
        }
    }

    template<typename F>
    void iterate_range_index(const int i, const int j, const int r, const int c, F f) {
        for (int ii = i, il = i + r; ii < il; ++ii) {
            for (int jj = j, jl = j + c; jj < jl; ++jj) {
                MATRIX(ii, jj) = f(ii, jj);
            }
        }
    }

    long long size() {
        return Rows * Cols;
    }

    void print(const char *format) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                printf(format, MATRIX(i, j));
            }
            printf("\n");
        }
    }

    std::vector<COORD_T > matrix_which_max() {
        std::vector<COORD_T > maxes;

        T max = -INFINITY;
        T v;

        if (size() <= 0) { return maxes; }

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                v = MATRIX(i, j);
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

        for (int j = 0; j < Cols; ++j) {
            for (int i = 0; i < Rows; ++i) {
                if (MATRIX(i, j) == val) { res.push_back({i, j}); }
            }
        }
        return res;
    }
};

#undef MATRIX

#endif //SIM_2D_CPP_MATRIX_H

/**
 * Highly customized matrix which is designed for Sim2D use.
 */

#ifndef SIM_2D_CPP_MATRIX_H
#define SIM_2D_CPP_MATRIX_H

#include <cassert>
#include <array>
#include <cmath>
#include <vector>
#include <iostream>
#include <map>

/* Coordinator Type: [X, Y] */
#define COORD_T std::array<int, 2>

/**
 * Abstract Matrix class.
 * @tparam T Type of the elements in matrix.
 */
template<typename T>
class Matrix {
public:
    virtual ~Matrix() = default;

    /**
     * @return The number of columns of the matrix.
     */
    [[nodiscard]] virtual int cols() const = 0;

    /**
     * @return The number of rows of the matrix.
     */
    [[nodiscard]] virtual int rows() const = 0;

    /**
     * Access an element in the matrix.
     * @param i The index of row in which the element presents.
     * @param j The index of column in which the element presents.
     * @return The element presents in the i-th row and j-th column.
     */
    virtual T &operator()(int i, int j) = 0;

    friend std::ostream &operator<<(std::ostream &os, Matrix<T> &that) {
        for (int i = 0, r = that.rows(); i < r; ++i) {
            for (int j = 0, c = that.cols(); j < c; ++j) {
                os << that.operator()(i, j) << " ";
            }
            os << std::endl;
        }
        return os;
    }

    /**
     * Set all elements in the matrix to specified value.
     * @param val Value of all elements to set.
     */
    void setAll(T val) {
        for (int i = 0, r = rows(); i < r; ++i) {
            for (int j = 0, c = cols(); j < c; ++j) {
                this->operator()(i, j) = val;
            }
        }
    }

    /**
     * Compare elements in two matrices.
     * @param that The second matrix to compare with.
     * @return true if two matrices have same elements at same positions.
     */
    bool equals(Matrix<T> *that) {
        if (this == that) { return true; }
        if (this->cols() != that->cols() || this->rows() != that->rows()) { return false; }

        for (int i = 0, r = rows(); i < r; ++i) {
            for (int j = 0, c = cols(); j < c; ++j) {
                if (this->operator()(i, j) != that->operator()(i, j)) { return false; }
            }
        }
        return true;
    }

    /**
     * Iterate through all elements in matrix and apply a function to each of them.
     * Equivalent to:
     *      for (int i = 0; i < rows; ++i) {
     *          for (int j = 0; j < cols; ++j) {
     *              f( &MATRIX[i][j] )
     *          }
     *      }
     * @tparam F Function type.
     * @param f  Function in form of: F f(T& element) { ... }
     *           The reference of each element will be passed to this function.
     */
    template<typename F>
    void iter(F f) {
        for (int i = 0, r = rows(); i < r; ++i) {
            for (int j = 0, c = cols(); j < c; ++j) {
                f(this->operator()(i, j));
            }
        }
    }

    /**
     * Iterate through all elements in matrix and apply a function to indexes of each element.
     * Equivalent to:
     *      for (int i = 0; i < rows; ++i) {
     *          for (int j = 0; j < cols; ++j) {
     *              f( i, j )
     *          }
     *      }
     * @tparam F
     * @param f
     */
    template<typename F>
    void iter_index(F f) {
        for (int i = 0, r = rows(); i < r; ++i) {
            for (int j = 0, c = cols(); j < c; ++j) {
                f(i, j);
            }
        }
    }

    template<typename F>
    void iter_range(const int i, const int j, const int r, const int c, F f) {
        for (int ii = i, il = i + r; ii < il; ++ii) {
            for (int jj = j, jl = j + c; jj < jl; ++jj) {
                f(this->operator()(ii, jj));
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
        for (int j = 0, c = cols(); j < c; ++j) {
            f(j);
        }
    }

    template<typename F>
    void iter_rows(F f) {
        for (int i = 0, r = rows(); i < r; ++i) {
            f(i);
        }
    }

    template<typename F>
    bool any(F f) {
        for (int i = 0, r = rows(); i < r; ++i) {
            for (int j = 0, c = cols(); j < c; ++j) {
                if (f(this->operator()(i, j))) { return true; }
            }
        }
        return false;
    }

    template<typename F>
    DBL_T sum() {
        DBL_T sum = 0;
        for (int i = 0, r = rows(); i < r; ++i) {
            for (int j = 0, c = cols(); j < c; ++j) {
                sum += this->operator()(i, j);
            }
        }
        return sum;
    }

    [[nodiscard]] long long size() const {
        return rows() * cols();
    }

    std::vector<COORD_T > matrix_which_max() {
        std::vector<COORD_T > maxes;

        T max = -INFINITY, v;

        if (size() <= 0) { return maxes; }

        for (int i = 0, r = rows(); i < r; ++i) {
            for (int j = 0, c = cols(); j < c; ++j) {
                v = this->operator()(i, j);
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

        for (int i = 0, r = rows(); i < r; ++i) {
            for (int j = 0, c = cols(); j < c; ++j) {
                if (this->operator()(i, j) == val) { res.push_back({i, j}); }
            }
        }
        return res;
    }
};

template<typename T, int ROWS, int COLS>
class MatrixS : public Matrix<T> {

private:
    T MATRIX[ROWS][COLS];

public:

    MatrixS() = default;

    explicit MatrixS(const T val) {
        Matrix<T>::setAll(val);
    }

    explicit MatrixS(Matrix<T> &that) {
        assert(that.rows() == ROWS && that.cols() == COLS);
        for (int i = 0; i < ROWS; ++i) {
            for (int j = 0; j < COLS; ++j) {
                MATRIX[i][j] = that.operator()(i, j);
            }
        }
    }

    T &operator()(const int i, const int j) {
        assert(0 <= i && i < ROWS && 0 <= j && j < COLS);
        return MATRIX[i][j];
    }

    [[nodiscard]] int cols() const {
        return COLS;
    }

    [[nodiscard]] int rows() const {
        return ROWS;
    }
};

template<typename T>
class MatrixD : public Matrix<T> {

private:
    T **MATRIX = nullptr;
    int ROWS = 0;
    int COLS = 0;

public:
    MatrixD(const int r, const int c) : ROWS(r), COLS(c) {
        MATRIX = new T *[ROWS];
        for (int i = 0; i < ROWS; ++i) {
            MATRIX[i] = new T[COLS];
        }
    }

    MatrixD(const int r, const int c, const T val) : MatrixD(r, c) {
        Matrix<T>::setAll(val);
    }

    ~MatrixD() {
        for (int i = 0; i < ROWS; ++i) {
            delete[] MATRIX[i];
        }
        delete[] MATRIX;
        MATRIX = nullptr;
        ROWS = 0;
        COLS = 0;
    }

    T &operator()(const int i, const int j) {
        assert(0 <= i && i < ROWS && 0 <= j && j < COLS);
        return MATRIX[i][j];
    }

    [[nodiscard]] int cols() const {
        return COLS;
    }

    [[nodiscard]] int rows() const {
        return ROWS;
    }
};

#endif //SIM_2D_CPP_MATRIX_H

#ifndef C_SRC_MATRIX_H
#define C_SRC_MATRIX_H

#include "../collection/collection.h"

#define MATRIX_ITR(r, c, body)                  \
    for (int _i_ = 0; _i_ < r; ++_i_) {         \
        for (int _j_ = 0; _j_ < c; ++_j_) {     \
            body                                \
        }                                       \
    }

#define MATRIX_ITR2(r, c, body)                 \
    for (int _j_ = 0; _j_ < c; ++_j_) {         \
        for (int _i_ = 0; _i_ < r; ++_i_) {     \
            body                                \
        }                                       \
    }

#define MATRIX_MAP(src, dst, r, c, f)           \
    for (int _i_ = 0; _i_ < r; ++_i_) {         \
        for (int _j_ = 0; _j_ < c; ++_j_) {     \
            dst[_i_][_j_] = f(src[_i_][_j_]);   \
        }                                       \
    }

#define MATRIX_INIT(src, r, c, val)             \
    for (int _i_ = 0; _i_ < r; ++_i_) {         \
        for (int _j_ = 0; _j_ < c; ++_j_) {     \
            src[_i_][_j_] = val;                \
        }                                       \
    }

#define MATRIX_PRINT(xs, r, c, pattern)         \
    for (int _i_ = 0; _i_ < r; ++_i_) {         \
        for (int _j_ = 0; _j_ < c; ++_j_) {     \
            printf(pattern, xs[_i_][_j_]);      \
        }                                       \
        printf("\n");                           \
    }

#define MATRIX_COPY(src, dst, r, c)             \
    for (int _i_ = 0; _i_ < r; ++_i_) {         \
        for (int _j_ = 0; _j_ < c; ++_j_) {     \
            dst[_i_][_j_] = src[_i_][_j_];      \
        }                                       \
    }

/**
 * Find the maximum value in the matrix and the count of them.
 * @param r number of rows
 * @param c number of columns
 * @return Pair(max_value, count_of_max_values)
 */
pair_t matrix_max(int r, int c, const double matrix[r][c]);

/**
 * Find the minimum value in the matrix and the count of them.
 * @param r number of rows
 * @param c number of columns
 * @return Pair(min_value, count_of_min_values)
 */
pair_t matrix_min(int r, int c, const double matrix[r][c]);

/**
 * Find the n_th element in the matrix with `element == value`
 * @param r number of rows
 * @param c number of columns
 * @return Pair(x, y) where matrix[x][y] indicates the location of found element
 */
pair_t matrix_find(int r, int c, const double matrix[r][c], double value, int n);

int matrix_count(int r, int c, const double matrix[r][c], double value);

#endif //C_SRC_MATRIX_H

#ifndef C_SRC_MATRIX_H
#define C_SRC_MATRIX_H

#include "../collection/collection.h"

#define MATRIX_ITR(i, r, j, c, body)            \
    for (i = 0; i < r; ++i) {                   \
        for (j = 0; j < c; ++j) {               \
            body                                \
        }                                       \
    }

#define MATRIX_MAP(src, dst, i, r, j, c, f)     \
    for (i = 0; i < r; ++i) {                   \
        for (j = 0; j < c; ++j) {               \
            dst[i][j] = f(src[i][j]);           \
        }                                       \
    }

#define MATRIX_PRINT(xs, i, r, j, c, pattern)   \
    for (i = 0; i < r; ++i) {                   \
        for (j = 0; j < c; ++j) {               \
            printf(pattern, xs[i][j]);          \
        }                                       \
        printf("\n");                           \
    }

#define MATRIX_COPY(src, dst, i, r, j, c)       \
    for (i = 0; i < r; ++i) {                   \
        for (j = 0; j < c; ++j) {               \
            dst[i][j] = src[i][j];              \
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
 * Find the n_th element in the matrix with `element == value`
 * @param r number of rows
 * @param c number of columns
 * @return Pair(x, y) where matrix[x][y] indicates the location of found element
 */
pair_t matrix_find(int r, int c, const double matrix[r][c], double value, int n);

#endif //C_SRC_MATRIX_H

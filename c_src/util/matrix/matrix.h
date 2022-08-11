#ifndef C_SRC_MATRIX_H
#define C_SRC_MATRIX_H

#define MATRIX_ITR(r, c, body)              \
    for (int i = 0; i < r; ++i) {           \
        for (int j = 0; j < c; ++j) {       \
            body                            \
        }                                   \
    }

#define MATRIX_MAP(src, dst, r, c, f)       \
    for (int i = 0; i < r; ++i) {           \
        for (int j = 0; j < c; ++j) {       \
            dst[i][j] = f(src[i][j]);       \
        }                                   \
    }

#define MATRIX_PRINT(x, r, c, pattern)      \
    for (int i = 0; i < r; ++i) {           \
        for (int j = 0; j < c; ++j) {       \
            printf(pattern, x[i][j]);       \
        }                                   \
        printf("\n");                       \
    }

#define MATRIX_COPY(src, dst, r, c)         \
    for (int i = 0; i < r; ++i) {           \
        for (int j = 0; j < c; ++j) {       \
            dst[i][j] = src[i][j];          \
        }                                   \
    }

#endif //C_SRC_MATRIX_H

#include "matrix.h"
#include <float.h>
#include <math.h>
#include <assert.h>

pair_t matrix_max(const int r, const int c, const double matrix[r][c]) {
    double   max   = -DBL_MAX, v;
    int      count = 0;
    for (int i     = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            v = matrix[i][j];
            if (v < max) { continue; }
            if (v == max) { count++; }
            else {
                count = 1;
                max   = v;
            }
        }
    }

    node_t _max, _count;
    _max._double = max;
    _count._int  = count;
    pair_t res;
    res.x = _max;
    res.y = _count;
    return res;
}

pair_t matrix_min(const int r, const int c, const double matrix[r][c]) {
    double   min   = DBL_MAX, v;
    int      count = 0;
    for (int i     = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            v = matrix[i][j];
            if (v > min) { continue; }
            if (v == min) { count++; }
            else {
                count = 1;
                min   = v;
            }
        }
    }

    node_t _min, _count;
    _min._double = min;
    _count._int  = count;
    pair_t res;
    res.x = _min;
    res.y = _count;
    return res;
}

pair_t matrix_find(const int r, const int c, const double matrix[r][c], const double value, int n) {
    assert(n > 0);
    pair_t res;
    res.x._int = (int) NAN;
    res.y._int = (int) NAN;
    for (int j = 0; j < c; ++j) {
        for (int i = 0; i < r; ++i) {
            if (matrix[i][j] == value) { n--; }
            if (n == 0) {
                res.x._int = i;
                res.y._int = j;
                return res;
            }
        }
    }
    return res;
}

int matrix_count(const int r, const int c, const double matrix[r][c], const double value) {
    int count = 0;
    double v;
    for (int i     = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            if (matrix[i][j] == value) { count++; }
        }
    }
    return count;
}
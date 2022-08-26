#include "collection.h"

#include <float.h>
#include <math.h>
#include <assert.h>

int seq_length_out(double *buf, double from, double to, const int length_out) {
    if (length_out <= 0) { return 0; }
    if (from > to) {
        double tmp = from;
        from = to;
        to   = tmp;
    }
    if (from < -DBL_MAX || to > DBL_MAX) { return 0; }

    double       value = from;
    const double by    = (to - from) / (length_out - 1);

    for (int i = 0; i < length_out; ++i) {
        buf[i] = value;
        value += by;
    }
    return length_out;
}

int seq_by(double *buf, double from, double to, double by) {
    if (from > to) {
        double tmp = from;
        from = to;
        to   = tmp;
    }
    if (from < -DBL_MAX || to > DBL_MAX) { return 0; }

    double value = from;
    int    i;
    for (i = 0; value < to; ++i) {
        buf[i] = value;
        value += by;
    }
    return i;
}

pair_t array_min(const int len, const double arr[len]) {
    double   min   = DBL_MAX;
    int      count = 0;
    for (int i     = 0; i < len; ++i) {
        double v = arr[i];
        if (v > min) { continue; }
        if (v == min) { count++; }
        else {
            count = 1;
            min   = v;
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

pair_t array_max(const int len, const double arr[len]) {
    double   max   = -DBL_MAX;
    int      count = 0;
    for (int i     = 0; i < len; ++i) {
        double v = arr[i];
        if (v < max) { continue; }
        if (v == max) { count++; }
        else {
            count = 1;
            max   = v;
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

int int_array_count(const int len, const int arr[len], const int value) {
    int      count = 0;
    for (int i     = 0; i < len; i++) {
        if (arr[i] == value) { count++; }
    }
    return count;
}

int double_array_find(const int len, const double arr[len], const double val, int n) {
    assert(n > 0);
    for (int i = 0; i < len; ++i) {
        if (arr[i] == val) { n--; }
        if (n == 0) {
            return i;
        }
    }
    return -1;
}

int int_array_find(const int len, const int arr[len], const int val, int n) {
    assert(n > 0);
    for (int i = 0; i < len; ++i) {
        if (arr[i] == val) { n--; }
        if (n == 0) {
            return i;
        }
    }
    return -1;
}

bool int_array_contain(const int len, const int arr[len], int value) {
    for (int i = 0; i < len; ++i) {
        if (value == arr[i]) { return true; }
    }
    return false;
}

int double_array_delete_many(const int len, double arr[len], int idx_len, const int indices[idx_len]) {
    for (int i = 0; i < idx_len; ++i) {
        arr[indices[i]] = NAN;
    }

    int next_pos = 0, deleted = 0;

    for (int i = 0; i < len; i++) {
        if (!isnan(arr[i])) {
            arr[next_pos] = arr[i];
            if (i != next_pos) { arr[i] = NAN; }
            ++next_pos;
        } else { ++deleted; }
    }

    return len - deleted;
}
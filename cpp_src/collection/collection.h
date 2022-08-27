#ifndef SIM_2D_CPP_COLLECTION_H
#define SIM_2D_CPP_COLLECTION_H

#include "matrix_ext.h"

template<typename T>
int seq_length_out(T *buf, T from, T to, const int length_out) {
    if (length_out <= 0) { return 0; }
    if (from > to) {
        long double tmp = from;
        from = to;
        to   = tmp;
    }
    assert(to > -INFINITY && from < INFINITY);

    T       value = from;
    const T by    = (to - from) / (length_out - 1.0);

    for (int i = 0; i < length_out; ++i) {
        buf[i] = value;
        value += by;
    }
    return length_out;
}

template<typename T>
int seq_by(T *buf, T from, T to, T by) {
    if (from > to) {
        T tmp = from;
        from = to;
        to   = tmp;
    }
    if (from < -INFINITY || to > INFINITY) { return 0; }

    double value = from;
    int    i;
    for (i = 0; value < to; ++i) {
        buf[i] = value;
        value += by;
    }
    return i;
}

template<typename T>
std::vector<int> vector_which_min(std::vector<T> vector) {
    std::vector<int> mins;

    T min = INFINITY;
    T v;

    if ((int) vector.size() <= 0) { return mins; }

    for (int i = 0, len = (int) vector.size(); i < len; ++i) {
        v = vector.at(i);
        if (v > min) { continue; }
        if (v < min) {
            mins.clear();
            min = v;
        }
        mins.push_back(i);
    }
    return mins;
}

template<typename T>
std::vector<int> vector_which_max(std::vector<T> vector) {
    std::vector<int> maxes;

    T max = -INFINITY;
    T v;

    if ((int) vector.size() <= 0) { return maxes; }

    for (int i = 0, len = (int) vector.size(); i < len; ++i) {
        v = vector.at(i);
        if (v < max) { continue; }
        if (v > max) {
            maxes.clear();
            max = v;
        }
        maxes.push_back(i);
    }
    return maxes;
}

template<typename T>
bool array_contains(T val, const int len, T *arr) {
    for (int i = 0; i < len; ++i) {
        if (val == arr[i]) { return true; }
    }
    return false;
}

template<typename T>
std::vector<T> vector_remove_many_by_index(std::vector<T> vector, const int len, int *indexes) {
    std::vector<T> res;

    for (int i = 0, l = vector.size(); i < l; ++i) {
        if (array_contains<int>(i, len, indexes)) { continue; }
        res.push_back(vector.at(i));
    }
    return res;
}

#endif //SIM_2D_CPP_COLLECTION_H

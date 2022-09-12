/**
 * Extended functions for collections.
 */

#ifndef SIM_2D_CPP_COLLECTION_H
#define SIM_2D_CPP_COLLECTION_H

#include "matrix.h"

/**
 * Generate a sequence in an array.
 * Equivalent to:
 *      seq(from = from, to = to, length.out = length_out)
 * @tparam T         type of elements.
 * @param buf        buffer: a T type array.
 * @param from       the starting value of the sequence.
 * @param to         the (maximal) end value of the sequence.
 * @param length_out desired length of the sequence.
 * @attention        desired length is not always equal to the size of array.
 * @return           The number of elements generated in the array.
 */
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

/**
 * Generate a sequence in an array.
 * Equivalent to:
 *      seq(from = from, to = to, by = by)
 * @tparam T   type of elements.
 * @param buf  buffer: a T type array.
 *             by = ((to - from)/(length.out - 1))
 * @param from the starting value of the sequence.
 * @param to   the (maximal) end value of the sequence.
 * @param by   number: increment of the sequence.
 * @return     The number of elements generated in the array.
 */
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

/**
 * Return indexes of those elements which are equal to specified value.
 * @tparam T   element type.
 * @tparam Len length of the std::array.
 * @param arr  std::array.
 * @param val  T type value.
 * @return     A vector of indexes.
 */
template<typename T, int Len>
std::vector<int> std_array_which_equals(const std::array<int, Len> &arr, const T val) {
    std::vector<int> res;

    if (Len <= 0) { return res; }

    for (int i = 0; i < Len; ++i) {
        if (arr[i] == val) {
            res.push_back(i);
        }
    }
    return res;
}

/**
 * Check if an array contains specified element.
 * @tparam T   element type.
 * @param val  T type value.
 * @param len  length of the array.
 * @param arr  T type array.
 * @return     True if this array contains at least one val.
 */
template<typename T>
bool array_contains(T val, const unsigned int len, T *arr) {
    for (int i = 0; i < len; ++i) {
        if (val == arr[i]) { return true; }
    }
    return false;
}

/**
 * Return indexes of those smallest elements in a vector.
 * @tparam T     element type.
 * @param vector T type vector.
 * @return       A vector of indexes.
 */
template<typename T>
std::vector<unsigned int> vector_which_min(const std::vector<T> &vector) {
    std::vector<unsigned int> mins;

    T min = INFINITY, v;

    if ((unsigned) vector.size() <= 0) { return mins; }

    for (unsigned int i = 0, len = (int) vector.size(); i < len; ++i) {
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

/**
 * Return indexes of those largest elements in a vector.
 * @tparam T     element type.
 * @param vector T type vector.
 * @return       A vector of indexes.
 */
template<typename T>
std::vector<unsigned int> vector_which_max(const std::vector<T> &vector) {
    std::vector<unsigned int> maxes;

    T max = -INFINITY, v;

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

/**
 * Remove multiple elements from a vector by their indexes.
 * @tparam T      element type.
 * @param vector  T type vector.
 * @param len     Number of indexes.
 * @param indexes Int array holding indexes of elements to remove.
 * @return A vector with specified elements removed.
 */
template<typename T>
std::vector<T> vector_remove_many_by_index(const std::vector<T> &vector, const unsigned int len, unsigned int *indexes) {
    std::vector<T> res;

    for (int i = 0, l = vector.size(); i < l; ++i) {
        if (array_contains<unsigned>(i, len, indexes)) { continue; }
        res.push_back(vector.at(i));
    }
    return res;
}

/**
 * Get the minimum value from a <K, DBL_T> std::map.
 * @tparam K  Key type.
 * @param map std::map.
 * @return minimum value in the std::map.
 *         INFINITY if map is empty.
 */
template<typename K>
DBL_T map_values_min(const std::map<K, DBL_T> &map) {
    DBL_T min = INFINITY;
    for (auto const &[_, d]: map) {
        min = d < min ? d : min;
    }
    return min;
}

/**Get the maximum value from a <K, DBL_T> std::map.
 * @tparam K  Key type.
 * @param map std::map.
 * @return maximum value in the std::map.
 *         -INFINITY if map is empty.
 */
template<typename K>
DBL_T map_values_max(const std::map<K, DBL_T> &map) {
    DBL_T max = (DBL_T) -INFINITY;
    for (auto const &[_, d]: map) {
        max = d > max ? d : max;
    }
    return max;
}

#endif //SIM_2D_CPP_COLLECTION_H

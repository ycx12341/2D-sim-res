#include "collection.h"

#include <float.h>

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
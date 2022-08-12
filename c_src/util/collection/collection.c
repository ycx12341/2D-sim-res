#include "collection.h"

#include <float.h>

int seq(double *buf, double from, double to, const int num) {
    if (num <= 0) { return 0; }
    if (from > to) {
        double tmp = from;
        from = to;
        to   = tmp;
    }
    if (from < -DBL_MAX || to > DBL_MAX) { return 0; }

    double       value  = from;
    const double offset = (to - from) / (num - 1);

    for (int i = 0; i < num; ++i) {
        buf[i] = value;
        value += offset;
    }
    return num;
}
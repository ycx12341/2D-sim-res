#include <stdio.h>

#include "util/math/math_ext.h"

int main() {
    /* testing code */
    set_seed(874513);

    printf("runif: %f\n", runif(1, 2));
    printf("runif: %f\n", runif(1, 2));
    printf("runif: %f\n", runif(1, 2));
    printf("runif: %f\n", runif(1, 2));
    printf("runif: %f\n", runif(1, 2));

    return 0;
}

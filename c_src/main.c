#include <stdio.h>

#include "util/math/math_ext.h"
#include "algorithm/scc.h"

int main() {
    /* testing code */
    set_seed(874513);

    printf("runif: %f\n", runif(1, 2));
    printf("runif: %f\n", runif(1, 2));
    printf("runif: %f\n", runif(1, 2));
    printf("runif: %f\n", runif(1, 2));
    printf("runif: %f\n", runif(1, 2));

    calculate_scc(NULL);

    return 0;
}

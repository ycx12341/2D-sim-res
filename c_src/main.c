#include <stdio.h>

#include "util/math/math_ext.h"
#include "util/matrix/matrix.h"
#include "algorithm/scc.h"
#include "util/collection/collection.h"

int main() {
    set_seed(SEED);

    double prob[6] = {0.1, 0.2, 0.1, 0.2, 0.1, 0.3};

    printf("%d \n", sample_prob1(6, prob));
    printf("%d \n", sample_prob1(6, prob));
    printf("%d \n", sample_prob1(6, prob));
    printf("%d \n", sample_prob1(6, prob));
    printf("%d \n", sample_prob1(6, prob));

//    sse_pars_t pars;
//    init_pars(&pars);
//    for (int i = 0; i < 1; ++i) {
//        calculate_scc(&pars, i);
//    }

    return 0;
}

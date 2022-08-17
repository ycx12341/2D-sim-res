#include <stdio.h>

#include "util/math/math_ext.h"
#include "util/matrix/matrix.h"
#include "algorithm/scc.h"
#include "util/collection/collection.h"

int main() {
    set_seed(SEED);

    sse_pars_t pars;
    init_pars(&pars);
    for (int i = 0; i < 2; ++i) {
        calculate_scc(&pars, i);
    }

    return 0;
}

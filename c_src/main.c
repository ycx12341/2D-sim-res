#include <stdio.h>

#include "util/math/math_ext.h"
#include "algorithm/scc.h"

int main() {
    set_seed(SEED);

    sse_pars_t pars;
    init_pars(&pars);
    calculate_scc(&pars, 0);

    return 0;
}

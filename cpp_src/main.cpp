#include <iostream>

#include "algo/sim2d.h"

int main() {
    set_seed(SEED);
    auto scc = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS);
    scc.calculate_sse();

    return 0;
}

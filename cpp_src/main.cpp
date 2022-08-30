#include <iostream>

#include "algo/sim2d.h"

int main() {
    auto scc = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
    scc.simulate();
    set_seed(SEED);
    return 0;
}

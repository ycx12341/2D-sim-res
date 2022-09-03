#include <iostream>

#include "algo/sim2d.h"

int main() {
    auto scc = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
    scc.pars->init();
    Parameters p_r2 = scc.simulate();

    auto scc_r2 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
    scc_r2.pars->load(p_r2);
    Parameters p_r3 = scc_r2.simulate();

    return 0;
}

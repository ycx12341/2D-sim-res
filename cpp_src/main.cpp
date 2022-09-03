#include <iostream>

#include "algo/sim2d.h"

int main() {
    auto scc = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
    scc.pars->init();
    Parameters p_r2 = scc.simulate(true);

    auto scc_r2 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
    scc_r2.pars->load(p_r2);
    Parameters p_r3 = scc_r2.simulate(true);

    auto scc_r3 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
    scc_r3.pars->load(p_r3);
    Parameters p_r4 = scc_r3.simulate(true);

    return 0;
}

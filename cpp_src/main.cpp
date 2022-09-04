#include <iostream>

#include "algo/sim2d.h"

int main() {
    auto scc = *Sim_2D_Factory::SCC_375(DEFAULT_N_DIMS, SEED);
    scc.pars->init();
    Parameters p_r2 = scc.simulate(MULTI_THREADING);

    auto scc_r2 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
    scc_r2.pars->load(p_r2);
    Parameters p_r3 = scc_r2.simulate(MULTI_THREADING);
//
//    auto scc_r3 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
//    scc_r3.pars->load(p_r3);
//    Parameters p_r4 = scc_r3.simulate(MULTI_THREADING);
//
//    auto scc_r4 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
//    scc_r4.pars->load(p_r4);
//    Parameters p_r5 = scc_r4.simulate(MULTI_THREADING);
//
//    auto scc_r5 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
//    scc_r5.pars->load(p_r5);
//    Parameters p_r6 = scc_r5.simulate(MULTI_THREADING);
//
//    auto scc_r6 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
//    scc_r6.pars->load(p_r6);
//    Parameters p_r7 = scc_r6.simulate(MULTI_THREADING);
//
//    auto scc_r7 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
//    scc_r7.pars->load(p_r7);
//    Parameters p_r8 = scc_r7.simulate(MULTI_THREADING);
//
//    auto scc_r8 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
//    scc_r8.pars->load(p_r8);
//    Parameters p_r9 = scc_r8.simulate(MULTI_THREADING);
//
//    auto scc_r9 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
//    scc_r9.pars->load(p_r9);
//    Parameters p_r10 = scc_r9.simulate(MULTI_THREADING);
//
//    auto scc_r10 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
//    scc_r10.pars->load(p_r10);
//    Parameters p_r11 = scc_r10.simulate(MULTI_THREADING);
//
//    auto scc_r11 = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
//    scc_r11.pars->load(p_r11);
//    Parameters p_r12 = scc_r11.simulate(MULTI_THREADING);

    return 0;
}

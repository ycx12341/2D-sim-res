#include <iostream>

#include "algo/sim2d.h"

int main() {
//    auto scc = *Sim_2D_Factory::SCC_375(DEFAULT_N_DIMS, SEED);
//    scc.pars->load_csv("Parameters_2022-09-04_19-11-29.csv");
//    Parameters p_r2 = scc.simulate(MULTI_THREADING);
//    p_r2.export_csv();
    constexpr static const unsigned X = SCC_X_LEN;
    constexpr static const unsigned Y = SCC_Y_LEN;
    constexpr static const unsigned N = DEFAULT_N_DIMS;

    auto scc = *Sim_2D_Factory<N, Y, X>::SCC_375(SEED);
    scc.pars.init();
    Parameters<N> p_r2 = scc.simulate(MULTI_THREADING);
//    p_r2.export_csv();
//
    auto scc2 = *Sim_2D_Factory<N, Y, X>::SCC_375(SEED);
    scc2.pars = p_r2;
    Parameters<N> p_r3 = scc2.simulate(MULTI_THREADING);

    auto scc3 = *Sim_2D_Factory<N, Y, X>::SCC_375(SEED);
    scc3.pars = p_r3;
    Parameters<N> p_r4 = scc3.simulate(MULTI_THREADING);


//    scc_r2.pars->load(p_r2);
//    Parameters p_r3 = scc_r2.simulate(MULTI_THREADING);
//    scc_r2.export_least_square();
//    scc_r2.export_summary();
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

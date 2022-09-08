#include "../algo/sim2d.h"

/* Terminate simulation if the decrease in averaged Least Square Differences has dropped below 5%. */
#define TERMINATE_CONDITION 0.05

void scc_simple() {
    Parameters p(DEFAULT_N_DIMS);
    DBL_T      mean_diff;

    /* Round 1
     *
     */
    auto scc_r1 = *Sim_2D_Factory::SCC_375(DEFAULT_N_DIMS, SEED);
    scc_r1.pars->init();
    p         = scc_r1.simulate(MULTI_THREADING);
    std::cout << p;
    auto scc = *Sim_2D_Factory::SCC_375(DEFAULT_N_DIMS, SEED);
    scc.pars->load(p);
    p = scc.simulate(MULTI_THREADING);
    std::cout << "A";
//    auto scc2 = *Sim_2D_Factory::SCC_375(DEFAULT_N_DIMS, SEED);
//    *scc2.pars = p;
//    p = scc2.simulate(MULTI_THREADING);
//
//    /* Round 2 ~ Round 3
//     *
//     */
//    for (int i = 2; i <= 3; ++i) {
//        std::cout << "A";
//        auto scc = *Sim_2D_Factory::SCC_375(DEFAULT_N_DIMS, SEED);
//        *scc.pars = p;
//        p = scc.simulate(MULTI_THREADING);
//    mean_diff = scc_r1.sum_diff / (DBL_T) scc_r1.nnan_idxs.size();
//    }
//
//    /* Round 4 ~ Round X
//     *
//     */
//    for (int i = 4;; ++i) {
//        std::cout << "B";
//        auto scc = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS, SEED);
//        *scc.pars = p;
//        p = scc.simulate(MULTI_THREADING);
//        if ((mean_diff - scc.sum_diff / (DBL_T) scc.nnan_idxs.size())
//            / mean_diff < TERMINATE_CONDITION) {
//            break;
//        }
//    mean_diff = scc_r1.sum_diff / (DBL_T) scc_r1.nnan_idxs.size();
//    }
}
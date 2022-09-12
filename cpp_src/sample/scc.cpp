/**
 * Possible SCC sample implementations using Sim-2D algorithms.
 */

#include "../algo/sim2d.h"

/* Terminate simulation if the decrease in averaged Least Square Differences has dropped below 5%. */
#define TERMINATE_CONDITION 0.05

constexpr static const unsigned X     = SCC_X_LEN;
constexpr static const unsigned X_375 = SCC_X_LEN_375;
constexpr static const unsigned Y     = SCC_Y_LEN;
constexpr static const unsigned Y_375 = SCC_Y_LEN_375;
constexpr static const unsigned N     = DEFAULT_N_DIMS;

typedef struct {
    DBL_T lb_bw;
    DBL_T ub_bw;
    DBL_T ess_target;
    DBL_T step_size;
}                               ROUND_SETTING_T;

static const unsigned        MAX_ROUND        = 10;
static const ROUND_SETTING_T ROUND[MAX_ROUND] = {
        {0.0, 2.0, 722,     0.01},
        {0.5, 2.5, 2500,    0.01},
        {1.0, 3.0, 2250,    0.01},
        {1.5, 3.5, 2025,    0.01},
        {2.5, 3.5, 1822.5,  0.01},
        {2.5, 3.5, 1640.25, 0.01},
        {2.5, 3.5, 2500,    0.01},
        {2.5, 3.5, 2250,    0.01},
        {2.5, 3.5, 2025,    0.01},
        {2.5, 3.5, 1822.5,  0.01}
};


void scc_simple_d3() {
    Parameters<N> *p;
    DBL_T         mean_before;

    /* Round 1 */
    auto scc_r1 = *Sim_2D_Factory<N, Y_375, X_375>::SCC_375(SEED);
    scc_r1.export_csv("ROUND1");
    scc_r1.pars->init();
    scc_r1.set_bw(ROUND[0].ess_target, ROUND[0].lb_bw, ROUND[0].ub_bw, ROUND[0].step_size);
    p = scc_r1.simulate(MULTI_THREADING);

    /* Round 2 ~ Round 3 */
    for (int i = 1; i <= 2; ++i) {
        auto scc = *Sim_2D_Factory<N, Y_375, X_375>::SCC_375(SEED);
        scc.export_csv("ROUND" + std::to_string(i + 1));
        scc.pars = p;
        scc.set_bw(ROUND[i].ess_target, ROUND[i].lb_bw, ROUND[i].ub_bw, ROUND[i].step_size);
        p           = scc.simulate(MULTI_THREADING);
        mean_before = scc.sum_diff / (DBL_T) scc.nnan_idxs.size();
    }

    /* Round 4 ~ Round X */
    for (int i = 3;; ++i) {
        auto scc = *Sim_2D_Factory<N, Y, X>::SCC(SEED);
        scc.export_csv("ROUND" + std::to_string(i + 1));
        scc.pars = p;

        if (i >= MAX_ROUND) {
            scc.set_bw(ROUND[MAX_ROUND - 1].ess_target, ROUND[MAX_ROUND - 1].lb_bw,
                       ROUND[MAX_ROUND - 1].ub_bw, ROUND[MAX_ROUND - 1].step_size);
        } else {
            scc.set_bw(ROUND[i].ess_target, ROUND[i].lb_bw, ROUND[i].ub_bw, ROUND[i].step_size);
        }

        p = scc.simulate(MULTI_THREADING);

        DBL_T mean_after = scc.sum_diff / (DBL_T) scc.nnan_idxs.size();
        if (mean_before - mean_after < TERMINATE_CONDITION * mean_before) {
            std::cout << "Terminate simulation: the decrease in averaged Least Square Differences has dropped below 5%."
                      << std::endl;
            break;
        }
        mean_before = mean_after;
    }
    delete p;
}
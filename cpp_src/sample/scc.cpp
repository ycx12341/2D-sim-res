#include "../algo/sim2d.h"

/* Terminate simulation if the decrease in averaged Least Square Differences has dropped below 5%. */
#define TERMINATE_CONDITION 0.05

constexpr static const unsigned X = SCC_X_LEN;
constexpr static const unsigned Y = SCC_Y_LEN;
constexpr static const unsigned N = DEFAULT_N_DIMS;

/**
 * Example: TODO
 */
void scc_simple() {
    Parameters<N> p;
    DBL_T         mean_before;

    /* Round 1
     *
     */
    auto scc_r1 = *Sim_2D_Factory<N, Y, X>::SCC_375(SEED);
    scc_r1.pars.init();
    p = scc_r1.simulate(MULTI_THREADING);

    p.export_csv();
    scc_r1.export_least_square();
    scc_r1.export_summary();

    /* Round 2 ~ Round 3
     *
     */
    for (int i = 2; i <= 3; ++i) {
        auto scc = *Sim_2D_Factory<N, Y, X>::SCC_375(SEED);
        scc.pars = p;
        p           = scc.simulate(MULTI_THREADING);
        mean_before = scc.sum_diff / (DBL_T) scc.nnan_idxs.size();

        p.export_csv();
        scc_r1.export_least_square();
        scc_r1.export_summary();
    }

    /* Round 4 ~ Round X
     *
     */
    for (int i = 4;; ++i) {
        auto scc = *Sim_2D_Factory<N, Y, X>::SCC(SEED);
        scc.pars = p;
        p = scc.simulate(MULTI_THREADING);

        p.export_csv();
        scc_r1.export_least_square();
        scc_r1.export_summary();

        DBL_T mean_after = scc.sum_diff / (DBL_T) scc.nnan_idxs.size();
        if (mean_before - mean_after < TERMINATE_CONDITION * mean_before) {
            std::cout << "Terminate simulation: the decrease in averaged Least Square Differences has dropped below 5%."
                      << std::endl;
            break;
        }
        mean_before = mean_after;
    }
}
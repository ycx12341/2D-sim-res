/**
 * Possible SCC sample implementations using Sim-2D algorithms.
 *
 * Sample bw Settings:
 * {0.0, 2.0, 722.000, 0.01},      // Round 1 375
 * {1.5, 3.0, 2500.00, 0.01},      // Round 2 375  2.22
 * {2.5, 4.0, 2250.00, 0.01},      // Round 3 375  3.08
 * {2.0, 3.5, 2025.00, 0.01},      // Round 4      2.80
 * {2.5, 4.0, 1822.50, 0.01},      // Round 5      3.06
 * {3.5, 5.0, 1640.25, 0.01},      // Round 6      4.16
 * {3.5, 5.0, 2500.00, 0.01},      // Round 7      4.00
 * {3.5, 5.0, 2250.00, 0.01},      // Round 8      4.02
 * {4.0, 5.5, 2025.00, 0.01},      // Round 9      4.44
 * {3.5, 5.5, 1822.50, 0.01}       // Round 10
 */

#include "../algo/sim2d.h"

/* Terminate simulation if the decrease in averaged Least Square Differences has dropped below 5%. */
#define TERMINATE_CONDITION 0.05

#define ESS_DECREMENT       0.9
#define DEFAULT_STEP        0.01

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
}                               BW_T;
static const BW_T               BW_r1 = {0.0, 2.0, 722.000, DEFAULT_STEP};
static const BW_T               BW_r2 = {1.5, 3.0, 2500.00, DEFAULT_STEP};
static const BW_T               BW_r7 = {3.5, 5.0, 2500.00, DEFAULT_STEP};
static BW_T                     G_bw  = BW_r1;

void step_bw(const int round) {
    if (round == 0) { G_bw = BW_r1; }
    else if (round == 1) { G_bw = BW_r2; }
    else if (round == 6) { G_bw = BW_r7; }
    else {
        G_bw.ess_target *= ESS_DECREMENT;
        G_bw.step_size = DEFAULT_STEP;
    }
    assert(G_bw.ess_target >= 0);
}

auto scc_48_28(const std::string &name, Parameters<N> *p) {
    return Sim_2D_Builder<N, Y_375, X_375>::SCC_375_Builder(SEED)
            .export_csv(name)
            .bw(G_bw.ess_target, G_bw.lb_bw, G_bw.ub_bw, G_bw.step_size)
#ifdef USE_PRELOAD_REF
            .ref_den(MatrixS<DBL_T, REF_DEN_ROWS, REF_DEN_COLS>::of(D3_REF_DEN))
#endif
            .load_pars(p)
            .build();
}

auto scc_60_35(const std::string &name, Parameters<N> *p) {
    return Sim_2D_Builder<N, Y, X>::SCC_Builder(SEED)
            .export_csv(name)
            .bw(G_bw.ess_target, G_bw.lb_bw, G_bw.ub_bw, G_bw.step_size)
#ifdef USE_PRELOAD_REF
            .ref_den(MatrixS<DBL_T, REF_DEN_ROWS, REF_DEN_COLS>::of(D3_REF_DEN))
#endif
            .load_pars(p)
            .build();
}

void scc_simple_d3() {
    Parameters<N> *p = nullptr;
    DBL_T         mean_before;

    /* Round 1 */
    step_bw(0);
    p = scc_48_28("ROUND1", p)->simulate(MULTI_THREADING);

    /* Round 2 ~ Round 3 */
    for (int i = 1; i <= 2; ++i) {
        step_bw(i);
        auto scc = scc_48_28("ROUND" + std::to_string(i + 1), p);
        p           = scc->simulate(MULTI_THREADING);
        mean_before = scc->sum_diff / (DBL_T) scc->nnan_idxs.size();
        delete scc;
    }

    /* Round 4 ~ Round X */
    for (int i = 3;; ++i) {
        step_bw(i);
        auto scc = scc_60_35("ROUND" + std::to_string(i + 1), p);
        p = scc->simulate(MULTI_THREADING);

        DBL_T mean_after = scc->sum_diff / (DBL_T) scc->nnan_idxs.size();
        delete scc;
        if (mean_before - mean_after < TERMINATE_CONDITION * mean_before) {
            std::cout << "Terminate simulation: the decrease in averaged Least Square Differences has dropped below "
                      << TERMINATE_CONDITION * 100 << "%."
                      << std::endl;
            break;
        }
        mean_before = mean_after;
    }

    /* Final Simulation Output */
    std::cout << "Final Simulation Output" << std::endl;
    Parameters<1U> pp = p->features_mean();
    std::cout << "Mean Parameters Value" << std::endl
              << pp << std::endl;

    auto scc = Sim_2D_Builder<1U, Y, X>::SCC_Builder(SEED)
            .export_csv("OUTPUT ROUND")
#ifdef USE_PRELOAD_REF
            .ref_den(MatrixS<DBL_T, REF_DEN_ROWS, REF_DEN_COLS>::of(D3_REF_DEN))
#endif
            .load_pars(&pp)
            .build();
    scc->calculate_sse();

    delete p;
}
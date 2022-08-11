#include "../util/math/math_ext.h"
#include "scc.h"

/* Lowerbound and upperbound of initial parameters */
#define DN_MIN                  0.000069
#define DN_MAX                  0.02
#define GAMMA_MIN               0.005
#define GAMMA_MAX               0.26
#define RN_MIN                  0.0008
#define RN_MAX                  0.08
#define ETA_MIN                 7.0
#define ETA_MAX                 18.0
#define DM_MIN                  0.0001
#define DM_MAX                  0.033
#define ALPHA_MIN               0.07
#define ALPHA_MAX               0.18
#define INIT_CELLS_COLS_MIN     1.0
#define INIT_CELLS_COLS_MAX     5.0
#define PROB_DEATH_MIN          0.01
#define PROB_DEATH_MAX          0.1
#define PROB_PROF_MIN           0.2
#define PROB_PROF_MAX           1.0

void init_pars(sse_pars_t* pars_p) {
    runif_seq(pars_p->dn, N_DIMS, DN_MIN, DN_MAX);
    runif_seq(pars_p->gamma, N_DIMS, GAMMA_MIN, GAMMA_MAX);
    runif_seq(pars_p->rn, N_DIMS, RN_MIN, RN_MAX);
    runif_seq(pars_p->eta, N_DIMS, ETA_MIN, ETA_MAX);
    runif_seq(pars_p->dm, N_DIMS, DM_MIN, DM_MAX);
    runif_seq(pars_p->alpha, N_DIMS, ALPHA_MIN, ALPHA_MAX);
    runif_seq(pars_p->init_cells_cols, N_DIMS, INIT_CELLS_COLS_MIN, INIT_CELLS_COLS_MAX);
    runif_seq(pars_p->prob_death, N_DIMS, PROB_DEATH_MIN, PROB_DEATH_MAX);
    runif_seq(pars_p->prob_prof, N_DIMS, PROB_PROF_MIN, PROB_PROF_MAX);
}
#ifndef C_SRC_ALGO_H
#define C_SRC_ALGO_H

#include "../config.h"
#include "../collection/collection.h"

#define SPACE_LENGTH_Y  ((1.0 / H) + 1.0)
#define SPACE_LENGTH_X  round(SPACE_LENGTH_Y * (280.0 / 480.0))

static const double H        = 1.0 / 59.0;
static const int    Y_LEN    = 60;                      // (1 / H) + 1
static const int    X_LEN    = 35;                      // round( Y_LEN * (280 / 480) )
static const int    MAT_SIZE = Y_LEN / 12;              // 5

static const double T              = 4.52;
static const double DT             = 0.0025;
static const int    TIME_STEPS     = (int) (T / DT);    // 1808
static const int    INT_TIME_STEPS = (int) (1 / DT);    // 400
static const int    DAY_TIME_STEPS = (int) 600;

typedef struct {
    /* parameters of the PDE model */
    double dn[N_DIMS];
    double gamma[N_DIMS];
    double rn[N_DIMS];
    double eta[N_DIMS];
    double dm[N_DIMS];
    double alpha[N_DIMS];

    double init_cells_cols[N_DIMS];     // initial columns of cells being set at the left boundary of the domain.
    double prob_death[N_DIMS];          // proportion of cells that will be undergo extinction at the end of every day.
    double prob_prof[N_DIMS];           // proportion of cells that will undergo mitosis at the end of every day.
}                   sse_pars_t;

double generate_pattern(const sse_pars_t *pars, int idx);

void calculate_scc(sse_pars_t *pars, int idx);

void init_pars(sse_pars_t *pars_p);

void proliferation(int PROF_CELLS_NUM, double prof_cells[PROF_CELLS_NUM], double ind_position[Y_LEN][X_LEN],
                   arraylist_t *coord);

void pde(int t, int idx, const sse_pars_t *pars, double f[Y_LEN][X_LEN], double m[Y_LEN][X_LEN]);

#endif //C_SRC_ALGO_H

#ifndef C_SRC_SCC_H
#define C_SRC_SCC_H

#include "../config.h"

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
} sse_pars_t;

void calculate_scc(sse_pars_t *pars, int idx);

void init_pars(sse_pars_t* pars_p);

#endif //C_SRC_SCC_H

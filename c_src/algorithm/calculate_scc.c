#include "scc.h"

/**
 * Calculates the least square difference between the simulated
 * invasion pattern and the reference invasion pattern.
 * @param pars Parameters
 */
void calculate_scc(sse_pars_t *pars, const int idx) {
    generate_pattern(pars, idx);
}
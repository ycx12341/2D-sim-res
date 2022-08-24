#include "scc.h"

#include <stdio.h>
#include <math.h>

/**
 * Calculates the least square difference between the simulated
 * invasion pattern and the reference invasion pattern.
 * @param pars Parameters
 */
void calculate_scc(sse_pars_t *pars, const int idx) {
    double v = generate_pattern(pars, idx);
    if (isnan(v)) printf("NAN\n");
    else printf(">> %f \n", v);
}
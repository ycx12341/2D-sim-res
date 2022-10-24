/**
 * Possible sample implementations using Sim-2D algorithms.
 */

#ifndef SIM_2D_CPP_SAMPLE_H
#define SIM_2D_CPP_SAMPLE_H

/**
 * Day 3 simulation of SCC invasion.
 *
 * Simulation Terminating Condition:
 *      If the decrease in averaged Least Square Differences has dropped below 5%.
 *
 * Simulation Outline:
 *      [Round 1] ----------------------------------------------------------------------------------
 *          Using 375 (48 * 28) grids for simulation.
 *          Generate the initial parameters and simulate.
 *          Lots of the parameters could lead to NaN Least Square Difference.
 *          Resample parameters to be used in next round from those valid parameters.
 *      [Round 2] [Round 3] ------------------------------------------------------------------------
 *          Using (48 * 28) grids for simulation.
 *          Load parameters resampled from previous round and simulate.
 *          Resample parameters to be used in next round from those valid parameters.
 *      [Round 4] ... [Round ?] --------------------------------------------------------------------
 *          Using (60 * 35) grids for simulation.
 *          Load parameters resampled from previous round and simulate.
 *          Keep simulating until terminated.
 */
void scc_simple_d3();

void test();

#endif //SIM_2D_CPP_SAMPLE_H

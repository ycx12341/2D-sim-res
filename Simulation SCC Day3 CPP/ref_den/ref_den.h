/**
 * Pre-loaded cell density.
 */

#ifndef SIM_2D_CPP_REF_DEN_H
#define SIM_2D_CPP_REF_DEN_H

#include "../config.h"

#ifdef USE_PRELOAD_REF

#define REF_DEN_ROWS 12
#define REF_DEN_COLS 7

/**
 * Pre-loaded cell density at day 3.
 *
 * t3.dat.p1 <- read.csv("Summary pic day 3 p1.csv")
 * t3.dat.p2 <- read.csv("Summary pic day 3 p2.csv")
 * t3.dat.p3 <- read.csv("Summary pic day 3 p3.csv")
 * t3.dat.p4 <- read.csv("Summary pic day 3 p4.csv")
 * t3.dat.p5 <- read.csv("Summary pic day 3 p5.csv")
 * t3.dat.p6 <- read.csv("Summary pic day 3 p6.csv")
 * t3.dat.p7 <- read.csv("Summary pic day 3 p7.csv")
 *
 * t3.ref.den <- cbind(
 *      rev(t3.dat.p1[, 5]), rev(t3.dat.p2[, 5]), rev(t3.dat.p3[, 5]),
 *      rev(t3.dat.p4[, 5]), rev(t3.dat.p5[, 5]), rev(t3.dat.p6[, 5]),
 *      rev(t3.dat.p7[, 5])
 * ) / 100
 */
static const DBL_T D3_REF_DEN[REF_DEN_ROWS][REF_DEN_COLS] = {
        {0.26500, 0.00000, 0.00000, 0, 0, 0, 0.00000},
        {0.35125, 0.00562, 0.00500, 0, 0, 0, 0.00000},
        {0.46812, 0.00000, 0.00000, 0, 0, 0, 0.00000},
        {0.48000, 0.00000, 0.00000, 0, 0, 0, 0.00000},
        {0.47938, 0.00000, 0.00000, 0, 0, 0, 0.00375},
        {0.48812, 0.00938, 0.00000, 0, 0, 0, 0.00000},
        {0.51750, 0.00000, 0.00000, 0, 0, 0, 0.00000},
        {0.38750, 0.00000, 0.00000, 0, 0, 0, 0.00000},
        {0.41812, 0.00000, 0.00000, 0, 0, 0, 0.00000},
        {0.44875, 0.00000, 0.00125, 0, 0, 0, 0.00000},
        {0.31875, 0.00000, 0.00875, 0, 0, 0, 0.00000},
        {0.44625, 0.00000, 0.00000, 0, 0, 0, 0.00000},
};

#endif

#endif //SIM_2D_CPP_REF_DEN_H

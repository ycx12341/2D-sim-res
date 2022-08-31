#ifndef SIM_2D_CPP_CONFIG_H
#define SIM_2D_CPP_CONFIG_H

#define DBL_T                   double

#define DEFAULT_N_DIMS          50
#define SEED                    234567          // R-style 6 digits unsigned integer

#define DN_MIN                  0.000069L
#define DN_MAX                  0.02L
#define GAMMA_MIN               0.005L
#define GAMMA_MAX               0.26L
#define RN_MIN                  0.0008L
#define RN_MAX                  0.08L
#define ETA_MIN                 7.0L
#define ETA_MAX                 18.0L
#define DM_MIN                  0.0001L
#define DM_MAX                  0.033L
#define ALPHA_MIN               0.07L
#define ALPHA_MAX               0.18L
#define INIT_CELLS_COLS_MIN     1.0L
#define INIT_CELLS_COLS_MAX     5.0L
#define PROB_DEATH_MIN          0.01L
#define PROB_DEATH_MAX          0.1L
#define PROB_PROF_MIN           0.2L
#define PROB_PROF_MAX           1.0L

#define POWER_MIN               0.0
#define POWER_MAX               2.0
#define POWER_STEP              0.01

#endif //SIM_2D_CPP_CONFIG_H

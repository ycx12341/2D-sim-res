#ifndef SIM_2D_CPP_CONFIG_H
#define SIM_2D_CPP_CONFIG_H

#define DBL_T                   double

#define DEFAULT_N_DIMS          50
#define SEED                    123456          // use system timestamp as see by setting to 0
#define USE_PRELOAD_REF                         // using preloaded T3_REF_DEN csv data

/* Initial Parameters */
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

/* SCC Default Settings */
#define SCC_H                   (1.0 / 59.0)
#define SCC_SPACE_LENGTH_Y      ((1.0 / SCC_H) + 1.0)
#define SCC_SPACE_LENGTH_X      (DBL_T) SCC_SPACE_LENGTH_Y * (280.0 / 480.0)
#define SCC_Y_LEN               (int) SCC_SPACE_LENGTH_Y
#define SCC_X_LEN               (int) SCC_SPACE_LENGTH_X
#define SCC_T                   4.52
#define SCC_DT                  0.0025
#define SCC_TIME_STEPS          (SCC_T / SCC_DT)
#define SCC_INT_TIME_STEPS      (1.0 / SCC_DT)
#define SCC_DAY_TIME_STEPS      600

#define POWER_MIN               0.0
#define POWER_MAX               2.0
#define POWER_STEP              0.01

#define ESS_TARGET              722

#define ABC_BCD_PAR_NUM         9
#define ABC_BCD_PAR_LB          { DN_MIN, GAMMA_MIN, RN_MIN, ETA_MIN, DM_MIN, ALPHA_MIN, INIT_CELLS_COLS_MIN, PROB_DEATH_MIN, PROB_PROF_MIN }
#define ABC_BCD_PAR_UB          { DN_MAX, GAMMA_MAX, RN_MAX, ETA_MAX, DM_MAX, ALPHA_MAX, INIT_CELLS_COLS_MAX, PROB_DEATH_MAX, PROB_PROF_MAX }
#define ABC_BCD_H               sqrt(1 - pow(0.05, 2))

#endif //SIM_2D_CPP_CONFIG_H

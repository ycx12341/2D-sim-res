#ifndef SIM_2D_CPP_CONFIG_H
#define SIM_2D_CPP_CONFIG_H

#define DBL_T                   double

#define DEFAULT_N_DIMS          20
#define SEED                    123456

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
#define ABC_BCD_PAR_LB          { 0.000069, 0.005, 0.0008,  7, 0.0001, 0.07, 1, 0.01, 0.2 }
#define ABC_BCD_PAR_UB          { 0.02    , 0.26 , 0.08  , 18, 0.033 , 0.18, 5, 0.1 , 1   }
#define ABC_BCD_H               sqrt(1 - pow(0.05, 2))

#endif //SIM_2D_CPP_CONFIG_H

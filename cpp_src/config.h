#ifndef SIM_2D_CPP_CONFIG_H
#define SIM_2D_CPP_CONFIG_H

#define DBL_T                       double      // { long double, double }

#define MULTI_THREADING             true        // { true, false }
#define DEFAULT_N_DIMS              10000
#define SEED                        0           // use system timestamp as see by setting to 0
#define USE_PRELOAD_REF                         // using preloaded T3_REF_DEN csv data
#define CONSOLE_REPORT                          // enable console report

/* Initial Parameters */
#define DN_MIN                      0.000069L
#define DN_MAX                      0.02L
#define GAMMA_MIN                   0.005L
#define GAMMA_MAX                   0.26L
#define RN_MIN                      0.0008L
#define RN_MAX                      0.08L
#define ETA_MIN                     7.0L
#define ETA_MAX                     18.0L
#define DM_MIN                      0.0001L
#define DM_MAX                      0.033L
#define ALPHA_MIN                   0.07L
#define ALPHA_MAX                   0.18L
#define INIT_CELLS_COLS_MIN         1.0L
#define INIT_CELLS_COLS_MAX         5.0L
#define PROB_DEATH_MIN              0.01L
#define PROB_DEATH_MAX              0.1L
#define PROB_PROF_MIN               0.2L
#define PROB_PROF_MAX               1.0L
#define ABC_BCD_PAR_LB              { DN_MIN, GAMMA_MIN, RN_MIN, ETA_MIN, DM_MIN, ALPHA_MIN, INIT_CELLS_COLS_MIN, PROB_DEATH_MIN, PROB_PROF_MIN }
#define ABC_BCD_PAR_UB              { DN_MAX, GAMMA_MAX, RN_MAX, ETA_MAX, DM_MAX, ALPHA_MAX, INIT_CELLS_COLS_MAX, PROB_DEATH_MAX, PROB_PROF_MAX }
#define ABC_BCD_H                   sqrt(1 - pow(0.05, 2))
#define PARS_NUM                    9

/* SCC (60 * 35) Settings */
#define SCC_H                       (1.0 / 59.0)
#define SCC_SPACE_LENGTH_Y          ((1.0 / SCC_H) + 1.0)
#define SCC_SPACE_LENGTH_X          (DBL_T) SCC_SPACE_LENGTH_Y * (280.0 / 480.0)
#define SCC_Y_LEN                   (int) SCC_SPACE_LENGTH_Y
#define SCC_X_LEN                   (int) SCC_SPACE_LENGTH_X
#define SCC_T                       4.52
#define SCC_DT                      0.0025
#define SCC_TIME_STEPS              (SCC_T / SCC_DT)
#define SCC_INT_TIME_STEPS          (1.0 / SCC_DT)
#define SCC_DAY_TIME_STEPS          600

/* SCC (48 * 28) Settings */
#define SCC_H_375                   (1.0 / 47.0)
#define SCC_SPACE_LENGTH_Y_375      ((1.0 / SCC_H_375) + 1.0)
#define SCC_SPACE_LENGTH_X_375      (DBL_T) SCC_SPACE_LENGTH_Y_375 * (280.0 / 480.0)
#define SCC_Y_LEN_375               (int) SCC_SPACE_LENGTH_Y_375
#define SCC_X_LEN_375               (int) SCC_SPACE_LENGTH_X_375
#define SCC_T_375                   4.52
#define SCC_DT_375                  0.004
#define SCC_TIME_STEPS_375          (SCC_T_375 / SCC_DT_375)
#define SCC_INT_TIME_STEPS_375      (1.0 / SCC_DT_375)
#define SCC_DAY_TIME_STEPS_375      375

#define POWER_MIN                   0.0
#define POWER_MAX                   2.0
#define POWER_STEP                  0.01

#define ESS_TARGET                  722

/* CSV settings */
#define EXPORT_CSV
#define CSV_SEPARATOR               ','
#define CSV_DBL_PRECISION           10
#define CSV_TIME_FORMAT             "%Y-%m-%d_%H-%M-%S"

#define CSV_DIR(name)               ((name) + std::string("_Results/"))
#define CSV_PARS_FNAME(name)        (CSV_DIR(name)  + std::string("Parameters_%Y-%m-%d_%H-%M-%S.csv"))
#define CSV_DIFF_FNAME(name)        (CSV_DIR(name)  + std::string("Differences_%Y-%m-%d_%H-%M-%S.csv"))
#define CSV_SMRY_FNAME(name)        (CSV_DIR(name)  + std::string("Summary_%Y-%m-%d_%H-%M-%S.csv"))
#define CSV_DETL_FNAME(name,idx)    (CSV_DIR(name)  + std::string("Idx_") + (idx) + std::string("_Details_%Y-%m-%d_%H-%M-%S.csv"))

#endif //SIM_2D_CPP_CONFIG_H

## CT on run 3 ##

This folder contains the detailed results of coverage test on the final parameter samples of run 3 for the simulations of post-day 3 pattern.

File **Reference densities SCC.rds** contains the reference dataset derived from the SCC invasion patterns in binary form. 

File **Round 9 parameters.txt** and **Round 9 Least Square.txt** contain the final parameter samples and their corresponding discrepancy measurement with the reference dataset. 

File **Round_9_paras15251_res.rds** contains the best-fitted output among all the outputs generated from the final parameter vectors. The standard deviation of the discrepancies between this best-fitted output and the reference dataset is used to add perturbations to the original simulation outputs in the coverage test.

File **Coverage test r9.R** and **PDE 2D ABC functions adjusted SCC (time-dependent parameters).R** carry out the coverage test on the final parameter samples of run 3. The detailed coverage test results are stored in **Coverage_test_disp_wt_prob_r11_pert.rds** and **Coverage_test_disp_wt_prob_r11_pert_non_zero.rds** were stored in folder **Cov_results_r9** and **Cov_results_r9_pert_all**. The corresponding coverage probabilities are stored in **Coverage probabilities r9 pert all.txt** and **Coverage probabilities glioma r9 non zero only**. At the end, the histogram plots of the coverage probabilities for each parameter are given in the .png files in the folder. 

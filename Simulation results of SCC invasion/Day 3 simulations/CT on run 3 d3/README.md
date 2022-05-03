## CT on run 3 d3 ##
This folder contains the detailed results of coverage test on the final parameter samples of run 3 for the simulations of post-day 3 pattern.

File **Summary pic day 3 p1-7.csv** contains the reference dataset derived from the post-day 3 SCC invasion patterns in binary form. 

File **Round 11 parameters.txt** and **Round 11 Least Square.txt** contain the final parameter samples and their corresponding discrepancy measurement with the reference dataset. The detailed simulation outputs (including the discrepancy measurements) are stored in folder **LS_results_r11**. 

File **Full output paras r11 run 1.rds** stores all the original simulation output of parameter vectors in **Round 11 parameters.txt**. 

File **Round_11_paras9900_res.rds** contains the best-fitted output among all the outputs generated from the final parameter vectors. The standard deviation of the discrepancies between this best-fitted output and the reference dataset is used to add perturbations to the original simulation outputs in the coverage test.

File **Coverage test r11.R** and **PDE 2D ABC functions adjusted.R** carry out the coverage test on the final parameter samples of run 3. The full coverage test results are stored in **Coverage_test_disp_wt_prob_r11_pert.rds** and **Coverage_test_disp_wt_prob_r11_pert_non_zero.rds**, with the detailed ones stored in folder **Cov_results_r11** and **Cov_results_r11_pert_all**. The corresponding coverage probabilities are stored in **Coverage probabilities r11.txt** and **Coverage probabilities glioma r11 non zero only**. At the end, the histogram plots of the coverage probabilities for each parameter are given in the .png files in the folder. 

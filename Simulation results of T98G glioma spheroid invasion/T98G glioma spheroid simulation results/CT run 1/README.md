## CT run 1 ##
This folder contains the detailed results of coverage test on the final parameter samples of run 1.

File **Glioma simulation full reference densities.rds** contains the reference dataset derived from the T98G glioma invasion patterns in binary form. 

File **Round 7 parameters.txt** and **Round 7 Least Square.txt** contain the final parameter samples and their corresponding discrepancy measurement with the reference dataset. The detailed simulation outputs (including the discrepancy measurements) are stored in folder **LS_results_r7**. 

File **Full output paras r7 run 1.rds** stores all the original simulation output of parameter vectors in **Round 7 parameters.txt**. 

File **Round_7_paras307_res.rds** contains the best-fitted output among all the outputs generated from the final parameter vectors. The standard deviation of the discrepancies between this best-fitted output and the reference dataset is used to add perturbations to the original simulation outputs in the coverage test.

File **Coverage test r7.R** and **PDE 2D ABC functions adjusted.R** carry out the coverage test on the final parameter samples of run 1. The full coverage test results are stored in **Coverage_test_disp_wt_prob_r7_d1_d3_pert_sep.rds** and **Coverage_test_disp_wt_prob_r7_d1_d3_pert_sep_non_zero.rds**, with the detailed ones stored in folder **Cov_results_r7** and **Cov_results_r7_pert_all**. The corresponding coverage probabilities are stored in **Coverage probabilities glioma r7 d1 d3 perturbed separately.txt** and **Coverage probabilities glioma r7 d1 d3 perturbed separately non zero only**. At the end, the histogram plots of the coverage probabilities for each parameter are given in the .png files in the folder. 

Folder **CT run 2** and **CT run 3** share similar structure with this folder. 

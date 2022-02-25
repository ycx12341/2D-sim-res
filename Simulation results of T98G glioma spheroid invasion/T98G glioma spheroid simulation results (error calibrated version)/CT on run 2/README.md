## CT on run 2 ##
This folder contains the detailed results of coverage test on the final parameter samples of run 2 - the run which produced the outputs chosen as the final simulation results.

File **Glioma simulation full reference densities.rds** contains the reference dataset derived from the T98G glioma invasion patterns in binary form. 

File **Round 6 parameters log transform.txt** and **Round 6 Least Square.txt** contain the final parameter samples and their corresponding discrepancy measurement with the reference dataset. 

File **Full output paras r6 run 2.rds** stores all the original simulation output of parameter vectors in **Round 6 parameters log transform.txt**. 

File **Round_6_paras2518_res.rds** contains the best-fitted output among all the outputs generated from the final parameter vectors. The standard deviation of the discrepancies between this best-fitted output and the reference dataset is used to generate the perturbations added to the original simulation outputs in the coverage test.

File **Coverage test r6.R** and **PDE 2D ABC functions adjusted.R** carry out the coverage test on the final parameter samples of run 2. The coverage probabilities are stored in **Coverage probabilities glioma r6 d1 d3 perturbed separately.txt** and **Coverage probabilities glioma r6 d1 d3 perturbed separately non zero only**. At the end, the histogram plots of the coverage probabilities for each parameter are given in the .png files in the folder. 

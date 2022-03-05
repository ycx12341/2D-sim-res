## Average results ##
This folder contains the final results from the three different runs of error calibrated ABC application on the reference dataset derived from the T98G glioma invasion patterns. 

File **Glioma simulation full reference densities.rds** contains the reference dataset derived from the T98G glioma invasion patterns in binary form.

File **Posterior results summary.R** reads in the **.rds** files **paras r6 average result run 1.rds**, **paras r6 average result run 2.rds**, **paras r6 average result run 3.rds**, the **.txt** files **Round 6 parameters log transform run 1.txt**, **Round 6 parameters log transform run 2.txt** and **Round 6 parameters log transform run 3.txt** which contain the final simulation outputs and their corresponding final parameter estimates from three different runs. It calculates the multivariate Bhattacharyya distance between the final outputs and the reference dataset, then decides which run's final simulation output and its corresponding parameter estimates will be chosen as the final result. The final parameter estimates from the 3 different runs and the ones chosen as final results are stored in **Glioma final parameter estimates err calib.rds**.  

Besides the numerical final outputs in the **.rds** files, their corresponding pattern plots are also given in the **.png** files in the folder (except **BW adaptive ESS mat ec.png** and **ESS adaptive ESS mat ec.png**).   

File **ESS_BW_plots.m** reads in the **.txt** files **Run 1 ESS BW.txt**, **Run 2 ESS BW.txt** and **Run 3 ESS BW.txt** which record the adaptive bandwidth factors and the actual ESS being used in each round of the simulation for the 3 different runs. The variations of these values for each round and each run are plotted in **ESS adaptive ESS mat ec.png** and **BW adpaitve ESS mat ec.png** to check the consistency of our ABC scheme.  

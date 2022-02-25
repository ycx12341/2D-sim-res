## Average results ##
This folder contains the final results from the three different runs of non error calibrated ABC application on the reference dataset derived from the T98G glioma invasion patterns. 

File **Posterior results summary.R** reads in the .rds files **paras r7 average result run 1.rds**, **paras r7 average result run 2.rds**, **paras r7 average result run 3.rds**, the .txt files **Round 7 parameters run 1.txt**, **Round 7 parameters run 2.txt** and **Round 7 parameters run 3.txt** which contain the final simulation outputs and their corresponding final parameter estimates from three different runs. It calculates the multivariate Bhattacharyya distance between the final outputs and the reference dataset, then decides which run's final simulation output and its corresponding parameter estimates will be chosen as the final result. The final parameter estimates from the 3 different runs and the ones chosen as final results are stored in **Glioma final parameter estimates non err calib.rds**.  

Besides the numerical final outputs in the .rds files, their corresponding pattern plots are also given in the .png files in the folder (except **BW adaptive ESS mat.png** and **ESS adaptive ESS mat.png**).   

File **ESS_BW_plots.m** reads in the .txt files **Run 1 ESS BW.txt**, **Run 2 ESS BW.txt** and **Run 3 ESS BW.txt** which record the adaptive bandwidth factors and the actual ESS being used in each round of the simulation for the 3 different runs. The variations of these values for each round and each run are plotted in **ESS adaptive ESS mat.png** and **BW adpaitve ESS mat.png** to check the consistency of our ABC scheme.  


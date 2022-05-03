## Average results ##
This folder contains the final results from the three different runs of non error calibrated ABC application on the reference dataset derived from the SCC invasion patterns. 

File **Summary pic day 3 p1-7.rds** contains the reference dataset derived from the SCC invasion patterns in binary form.

File **paras r11 d3 pattern plot run 1-3 correct init.png** and **grey shades plot paras r11 day 3 run 1-3.png** contains the plot of the final simulated patterns from the 3 different runs, presented in the form of both individual cell plots and grey shades density plots.  

File **Posterior results summary.R** reads in the .rds files **paras r11 average result run 1.rds**, **paras r11 average result run 2.rds**, **paras r11 average result run 3.rds**, the .txt files **Round 11 parameters run 1.txt**, **Round 11 parameters run 2.txt** and **Round 11 parameters run 3.txt** which contain the final simulation outputs and their corresponding final parameter estimates from three different runs. It calculates the Bhattacharyya distance between the final outputs and the reference dataset, then decides which run's final simulation output and its corresponding parameter estimates will be chosen as the final result. The chosen run's final simulation output is also used as the initial condition for the simulation of the next pattern. The final parameter estimates from the 3 different runs and the ones chosen as final results are stored in **Day 3 parameter estimates.rds**.  

Besides the numerical final outputs in the .rds files, their corresponding pattern plots are also given in the .png files in the folder (except **BW adaptive ESS day 3 mat.png** and **ESS adaptive ESS day 3 mat.png**).   

File **ESS_BW_plots.m** reads in the .txt files **Run 1 ESS BW.txt**, **Run 2 ESS BW.txt** and **Run 3 ESS BW.txt** which record the adaptive bandwidth factors and the actual ESS being used in each round of the simulation for the 3 different runs. The variations of these values for each round and each run are plotted in **ESS adaptive ESS day 3 mat.png** and **BW adpaitve ESS day 3 mat.png** to check the consistency of our ABC scheme.  

## Posterior analysis (run 3) ##
This folder contains the comparison plots of the prior and posterior density for each parameter. 

File **Posterior_comparisons.m** reads in the **.txt** files **Round 1 initial time varying parameters.txt**, **Round 9 parameters.txt**, **Time-dependent parameters bounds.txt** which stores the prior and posterior samples for each parameter and the lower&upper bounds for certain parameters' prior distributions. It then generates the density plots and save them as **.png** files in the main folder. 

File **Overlap SCC time-dependent.R** reads in the **.txt** files **Round 1 initial time varying parameters.txt**, **Round 9 parameters.txt** which stores the prior and posterior samples for each parameter and calculates the prior-posterior overlap values.

## Day 3 simulations ## 
This folder contains all the necessary code and simulation output regarding the application of non-error calibrated ABC on the reference dataset derived from the SCC
invasion patterns. Due to the stochastic nature of cell movements in the individual based model (IBM), the same analysis was carried out three times with different seeds being used in the random number generator in R (**set.seed()**). 

Once the final parameter samples are obtained for each run, we take the their sample means as the final parameter estimates in each corresponding run. We then calculate the Bhattacharyya distance between the simulation outputs generated from the 3 sets of final parameter estimates and the reference dataset. The one that yields the minimum Bhattacharyya distance and its corresponding parameter estimates is chosen as the final results for the current pattern. It is also be used as the initial condition for the evaluation of the next pattern. We then carry out the posterior diagnostic - "coverage test" on the final parameter samples from the run chosen as the final results of the current pattern to see if they can be considered as valid posteriors for the current simulations.  

Folder **Run 1-3** contains the detailed code and the simulation output. 

Folder **CT on run 3 d3** contains the coverage test results on the final parameter samples from the run chosen as the final result. 

Folder **Average results** contains the post-processing of the final results for each run. 

Folder **Day 6 simulations**, **Day 9 simulations**, **Day 12 simulations** and **Day 14 simulations** shares the same structure with the folder.  

See the README.md file for each folder for more details. 

### Note

All simulation results were generated using R **4.0.3** "Bunny-Wunnies Freak out". 

For R version **greater than 4.0.3**, the if statement in line 883 returns an error message instead of a warning message. We recommend users to change it to `if(length(temp.den.table) == 1)` so valid results can be obtained.

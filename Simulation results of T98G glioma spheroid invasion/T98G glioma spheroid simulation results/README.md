## T98G glioma spheroid simulation results ## 
This folder contains all the necessary code and simulation output regarding the application of non-error calibrated ABC on the reference dataset derived from the T98G glioma 
invasion patterns. Due to the stochastic nature of cell movements in the individual based model (IBM), the same analysis was carried out three times with different seed being used in the random number generator in R (**set.seed()**). Once the final parameter samples are obtained for each run, we carried out the posterior diagnostic - "coverage test" on them to see if they can be considered as valid posteriors from the perspective of Bayesian inference.  

Folder **Run 1-3** contains the detailed code and the simulation output. 

Folder **CT run 1-3** contains the coverage test results on the final parameter samples from the corresponding run. 

Folder **Average results** contains the post-processing of the final results for each run. 

See the README.md file for each folder for more details. 

All simulation results were generated using R 4.0.3 "Bunny-Wunnies Freak out".  

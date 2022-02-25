## T98G glioma spheroid simulation results ## 
This folder contains all the necessary code and simulation output regarding the application of error calibrated ABC on the reference dataset derived from the T98G glioma 
invasion patterns. Due to the stochastic nature of cell movements in the individual based model (IBM), the same analysis was carried out three times with different seed being used in the random number generator in R (**set.seed()**). The posterior diagnostics - "Coverage test" is only conducted on the final parameter samples from the run which was chosen as the final results.

Folder **Run 1-3** contains the detailed code and the simulation output. 

Folder **CT on run 2** contains the coverage test results on the final parameter samples from the run which was chosen as the final results.. 

Folder **Average results** contains the post-processing of the final results for each run. 

Folder **Posterior comparisons** contains the comparison between the posteriors obtained using non-error calibrated ABC and the ones obtained using error calibrated ABC.  

See the README.md file for each folder for more details. 

All simulation results were generated using R 4.0.3 "Bunny-Wunnies Freak out".  

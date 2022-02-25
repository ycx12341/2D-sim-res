## T98G glioma spheroid simulation results ## 
This folder contains all the necessary code and simulation output regarding the application of error calibrated ABC on the reference dataset derived from the T98G glioma 
invasion patterns. Due to the stochastic nature of cell movements in the individual based model (IBM), the same analysis was carried out three times with different seed being used in the random number generator in R (**set.seed()**). The posterior diagnostics - "Coverage test" is only conducted on the final parameter samples from the run which was chosen as the final results.

Folder **Run 1-3** contain the detailed code and the simulation output. 

Folder **CT run 2** contain the coverage test results on the final parameter samples from the run which was chosen as the final results.. 

Folder **Average results** contain the post-processing of the final results for each run. 

See the README.md file for each folder for more details. 

All simulation results were generated using R 4.0.3 "Bunny-Wunnies Freak out".  

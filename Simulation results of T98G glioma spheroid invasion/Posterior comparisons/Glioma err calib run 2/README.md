## Glioma err calib run 2 ## 
This folder contains the plots of probability densities for the parameter samples obtained at the end of each round regarding the simulation study of T98G glioma invasion patterns, calibrated using error calibrated ABC. 

The .txt files contain the parameter samples obtained at the end of each round. 

File **Overlap EC.R** calculates the prior-posterior overlap between the prior and posterior samples for each parameter. Low overlap values suggest the corresponding posteriors obtained were more influenced by the model and data, while high overlap values suggest the posteriors obtained were more influenced by the priors. 

File **Glioma_posteriors_run_2** plots the probability densities for the parameter samples obtained at the end of each round of the simulation study. It also plots the comparison between prior and posterior for each parameter within the model. 

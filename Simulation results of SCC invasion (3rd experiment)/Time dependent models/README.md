## Time dependent models ##
This folder contains the regression models fitted to the parameter estimates of each SCC pattern obtained at the end of the second experiment regarding the simulation study of SCC invasion patterns. 

File **Regression models.R** reads in the .rds files that contain the parameter estimates from the 2nd experiment (**Day x parameter estimates.rds**), fits different types of regression models to them (mean, linear, quadratic) and select the best-fitting model. The detailed output of these selected models were stored in **x regression models.rds** and their fitted values were stored in **Estimated values x.txt**. 

File **Parameter_estimations.m** reads in the final parameter estimates of the 2nd experiment stored in **Parameter estimates.txt** and the estimated regression values stored in **Estimated values x.txt**, then generates the comparison plot between the parameter estimates and the fitted regression values. The plots were stored in the .png files in the current folder. 

## Run 3 ##
This folder contains the necessary code and simulation output regarding the third run of the ABC application on the reference dataset derived from the SCC invasion patterns. Time-dependent parameters were introduced to the individual-based model. 

File **grey shades plot day 3 ref SCC.png**, **grey shades plot day 6 ref SCC.png**, **grey shades plot day 9 ref SCC.png**, **grey shades plot day 12 ref SCC.png**, **grey shades plot day 14 ref SCC.png**, **tumour pic p1 binary.png**, **tumour pic p2 binary.png**, **tumour pic p3 binary.png**, **tumour pic p4 binary.png** and **tumour pic p5 binary.png** are the reference SCC patterns presented in binary form and grey shades density regions.  

Folder **Reference densities** and file **Reference densities SCC.rds** contain the reference SCC cell densities of the patterns at different periods. 

File **alpha quadratic regression model.rds**, **dm constant regression model.rds**, **dn quadratic regression model.rds**, **gamma quadratic regression model.rds**, **prob.death constant regression model.rds**, **prob.prof quadratic regression model.rds** and **rn quadratic regression model.rds** contain the chosen regression models fitted to the final parameter estimates of the 2nd experiment regarding the simulation study of SCC invasion patterns. 

File **Initial time varying paras.R** generates the parameter vectors to be evaluated in the first round and stores them in **Round 1 initial time varying parameters.txt**. 

The main files **paras sampl rx.R** and function file **PDE 2D ABC functions adjusted (time-dependent parameters).R** evaluates the parameters of the current round (**Round x parameters parameters.txt**), record the discrepancy measurements in **Round x Least Sqaure.txt**, then generates the parameters to be evaluated in the next round (**Round x+1 parameters.txt**). 

File **Round x information list.rds** contains the parameter weights, tempering factors which further rescale these weights into proper resampling probabilities and the actual effective sample size for the parameters being evaluated in the current round. 

After the final parameter samples are obtained, their averaged values were taken and substituted back to the pattern generation function in order to produce the final results. The numerical outputs are stored in **paras r9 average results.rds**, the pattern plots are shown in **paras r9 tv post d3 pattern plot.png**, **paras r9 tv post d6 pattern plot.png**, **paras r9 tv post d9 pattern plot.png**, **paras r9 tv post d12 pattern plot.png**, **paras r9 tv post d14 pattern plot.png**, **grey shades paras r9 post d3.png**, **grey shades paras r9 post d6.png**, **grey shades paras r9 post d9.png**, **grey shades paras r9 post d12.png** and **grey shades paras r9 post d14.png**. 

File **Full parameter estimates r9 run 3.R** reads in the parameter samples obtained at the end of every round, stores the averaged parameter estimates at the end of every round in **Full parameter estimates 9 rounds run 3.txt** and the averaged final parameter estimates in **Round 9 final full parameter estimates run 3.txt**. Furthermore, it fits the corresponding regression models to the final parameters estimates of the time-dependent parameters, these fitted regression values are stored in **Estimated final parameter values regression.txt**. 

File **Parameter_estimations.m** reads in the final parameter estimates and the fitted regression values from **Round 9 final full parameter estimates run 3.txt** and **Estimated final parameter values regression.txt** and plots them in **paras r9 dn ests run 3.png**, **paras r9 gamma ests run 3.png**, **paras r9 rn ests run 3.png**, **paras r9 eta ests run 3.png**, **paras r9 dm ests run 3.png**, **paras r9 alpha ests run 3.png**, **paras r9 p_ext ests run 3.png** and **paras r9 p_mit ests run 3.png**.

File **Run 3 posterior summary.R** reads in the **Round x information list log transform.rds** files, extract the actual effective sample sizes and bandwidth factors for each round and store them in **Run 3 ESS BW.txt**. 

Folder **Run 1** and **Run 2** share similar structures with this folder. 

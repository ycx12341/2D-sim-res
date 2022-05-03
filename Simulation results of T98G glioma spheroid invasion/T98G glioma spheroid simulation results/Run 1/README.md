## Run 1 ##
This folder contains the necessary code and simulation output regarding the first run of the non-error calibrated ABC application on the reference dataset derived from the T98G glioma invasion patterns. 

File **pattern glioma t1 -40um 400_400 binary.png**, **pattern glioma t3 -40um 400_400 binary.png**, **grey shades plot day 1 ref.png** and **grey shades plot day 3 ref.png** are the reference T98G invasion patterns presented in binary form and grey shades density regions.  

File **Glioma simulation full reference densities.rds** contains the reference dataset derived from the T98G glioma invasion patterns in binary form.  

The main files **paras sampl rx.R** and function file **PDE 2D ABC functions adjusted.R** evaluates the parameters of the current round (**Round x parameters.txt**), record the discrepancy measurements in **Round x Least Sqaure.txt** and generates the parameters to be evaluated in the next round (**Round x+1 parameters.txt**). 

File **Round x information list.rds** contains the parameter weights, bandwidth factor which translates them into resampling probabilities and the actual effective sample size for the parameters being evaluated in the current round. 

After the final parameter samples were obtained, their averaged values were taken and substituted back to the pattern generation function in order to produce the final results. The numerical outputs are stored in **paras r7 average results.rds**, the pattern plots are shown in **paras r7 d1 pattern plot run 1.png**, **paras r7 d3 pattern plot run 1.png**, **grey shades plot day 1 paras r7 run 1.png** and **grey shades plot day 3 paras r7 run 1.png**. 

File **Paras ests.R** calculates the sample mean for the parameters being evaluated in each round. 

File **Run 1 posterior summary.R** reads in the **Round x information list.rds** files, extracts the actual effective sample sizes and bandwidth factors for each round and store them in **Run 1 ESS BW.txt**. 

Folder **Run 2** and **Run 3** share similar structures with this folder. 

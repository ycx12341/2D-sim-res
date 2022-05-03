## Run 1 ##
This folder contains the necessary code and simulation output regarding the first run of the non-error calibrated ABC application on the reference dataset derived from the SCC invasion patterns. 

File **tumour pic p1 binary.png** and **grey shades plot day 3 ref SCC.png** are the reference SCC invasion patterns presented in binary form and grey shades density regions.  

File **Summary pic day 3 p1-7.csv** contains the reference dataset derived from the T98G glioma invasion patterns in binary form.  

The main files **paras sampl rx.R** and function file **PDE 2D ABC functions adjusted.R** evaluates the parameters of the current round (**Round x parameters.txt**), record the discrepancy measurements in **Round x Least Sqaure.txt** and generates the parameters to be evaluated in the next round (**Round x+1 parameters.txt**). 

File **Round x information list.rds** contains the parameter weights, bandwidth factor which translates them into resampling probabilities and the actual effective sample size for the parameters being evaluated in the current round. 

After the final parameter samples are obtained, their average values were taken and substituted back to the pattern generation function in order to produce the final results. The numerical outputs are stored in **paras r11 average results.rds**, the pattern plots are shown in **paras r11 d3 pattern plot run 1.png**, **grey shades plot day 3 paras r11 run 1.png**. 

File **Paras ests d3 run 1.R** calculates the sample mean for the parameters being evaluated in each round. 

File **Run 1 posterior summary.R** reads in the **Round x information list.rds** files, extract the actual effective sample sizes and bandwidth factors for each round and store them in **Run 1 ESS BW.txt**. 

Folder **Run 2** and **Run 3** share similar structures with this folder. 

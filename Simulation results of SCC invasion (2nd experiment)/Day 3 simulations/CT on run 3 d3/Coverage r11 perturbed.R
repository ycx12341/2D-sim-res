# Coverage r11 perturbed.R
# This .R file conducts the coverage test on the final parameter samples 
# obtained at the end of the simulation for post-day 3 SCC invasion pattern to
# check if the final parameter samples can be regarded as appropriate posterior
# samples from the perspective of Bayesian inference. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(readr)
library(doParallel)
library(latex2exp)

# Set the directory to store the results. 
# save.sims.dir <- Cov_results_r11
save.sims.dir <- "Cov_results_r11_pert_all"
save.sims <- TRUE
if (save.sims) {
  if (!dir.exists(save.sims.dir)) dir.create(save.sims.dir)
}

# Set the cluster. 
n.thread <- detectCores()/2
cl <- makeCluster(n.thread)
registerDoParallel(cl)

# Read in the final parameter values. 
paras.r11 <- read.table("Round 11 parameters.txt", sep = "", header = TRUE)

# Directory which stores the simulation output. 
sims.dir <- "LS_results_r11"

# Check if there are still NaNs in the final output.
ls.r11 <- read.table("Round 11 Least Square.txt", 
                     sep = "", header = TRUE)
length(which(is.na(ls.r11[,2]))) # 22 singular values.  

# Remove the singular simulation output and their corresponding parameter 
# vectors. 
ls.r11.valid <- ls.r11[-which(is.na(ls.r11[,2])), ]
paras.r11.valid <- paras.r11[-which(is.na(ls.r11[,2])), ]

# Read in the file which contains the source functions and the observed 
# densities.  
source("PDE 2D ABC functions adjusted SCC.r")

# Read in the unperturbed model output and store them in a list. 
paras.r11.output <- vector(mode = "list")
for (i in 1:length(ls.r11.valid[,2])) {
  res.temp <- read_rds(paste0("./", sims.dir, 
                              "/Round_11_paras",
                              ls.r11.valid[i,1],"_res.rds"))
  ls.temp <- list(den.mat.d3 = res.temp$den.mat.d3)
  paras.r11.output[[i]] <- ls.temp
  print(i)
  names(paras.r11.output)[i] <- paste0("ls_r11_",ls.r11.valid[i,1])
}

# Write the unperturbed model output into a .rds file. 
write_rds(paras.r11.output, "Full output paras r11 run 3 valid.rds")

# Read in the output in round 11 that gave the minimum result of least square
# difference, and calculate the standard deviation (SD) based on it.
ls.r11.min <- read_rds("Round_11_paras9900_res.rds")
sd.r11.d3.min <- sd(ls.r11.min$den.mat.d3 - t3.ref.den)

# Add perturbations (gaussian) to the final unperturbed output based on the
# SD calculated in the previous step. 
paras.r11.output.pert <- vector(mode = "list")
set.seed(874512)
RNGkind(sample.kind = "Rejection")
for (i in 1:length(paras.r11.output)) {
  ls.temp <- paras.r11.output[[i]]
  # Read the unperturbed data
  den.mat.d3.temp.unpert <- ls.temp$den.mat.d3
  # Create empty matrices to store perturbed data
  den.mat.d3.temp.pert <- matrix(0, nrow = nrow(den.mat.d3.temp.unpert), 
                                 ncol = ncol(den.mat.d3.temp.unpert))
  # Add perturbations to the unperturbed output: (two different ways)
  # 1. Add perturbations to non-zero data only.
  # 2. Add perturbations to the full dataset. 
  # We believed the first method of adding perturbation is more reasonable, as
  # the SCC starts invading the surroundings from the left to the right. In 
  # the reference data (especially for the first few patterns), many regions 
  # close to the right boundary in the domain are still unoccupied. Hence, 
  # adding perturbation to all data points would corrupt the nature of the 
  # original dataset. 
  for (j in 1:length(den.mat.d3.temp.unpert)) {
    #if(den.mat.d3.temp.unpert[j] != 0) {
     den.mat.d3.temp.pert[j] <- rnorm(1, den.mat.d3.temp.unpert[j], sd = sd.r11.d3.min)
     while (den.mat.d3.temp.pert[j] < 0 || den.mat.d3.temp.pert[j] > 1) {
       den.mat.d3.temp.pert[j] <- rnorm(1, den.mat.d3.temp.unpert[j], sd = sd.r11.d3.min)
     }
    #}
  }
  
  # Sum of squared differences between the perturbed densities and the 
  # observed densities.
  sse.d3.temp <- sum((t3.ref.den - den.mat.d3.temp.pert)^2)
  # Store each set of perturbed data into the list
  paras.r11.output.pert[[i]] <- list(den.mat.d3.pert = den.mat.d3.temp.pert, 
                                    diff = sse.d3.temp)
  # Name the elements in the list with correct indices
  names(paras.r11.output.pert)[i] <- paste0("ls_r11_",ls.r11.valid[i,1],"_pert")
  print(i)
}

# Sort the sum of squared differences
ls.r11.valid.pert <- matrix(0, nrow = length(ls.r11.valid[,2])
                     ,ncol = 2)
for (i in 1:length(ls.r11.valid.pert[,2])) {
  ls.r11.valid.pert[i,] <- c(ls.r11.valid[i, 1], 
                            as.double(paras.r11.output.pert[[i]]$diff))
}
ls.r11.valid.pert.sort <- order(ls.r11.valid.pert[,2])
ls.r11.valid.pert.sort.min <- ls.r11.valid.pert.sort[1:200]

# Index in the full dataset (from 1 - 10000)
ls.r11.pert.sort.min.index <- ls.r11.valid.pert[ls.r11.valid.pert.sort[1:200], 
                                                1]

# Locate the corresponding 200 perturbed model outputs that have the 
# minimum discrepancy with the reference data.
paras.r11.output.pert.full <- paras.r11.output.pert
paras.r11.output.pert <- paras.r11.output.pert[ls.r11.valid.pert.sort.min]

# Discrepancies between unperturbed model output with each of these 200 
# perturbed model output which were treated as new reference data. 
ls.r11.diff.mat <- vector(mode = "list")
for (i in 1:length(paras.r11.output.pert)) {
  # For each perturbed dataset, record its discrepancy with other unperturbed 
  # model output, corresponding weights and resampling probabilities. 
  diff.mat <- matrix(0, nrow = length(paras.r11.valid[,2]), ncol = 2)
  paras.r11.output.pert.temp <- paras.r11.output.pert[[i]]
  for (j in 1:length(paras.r11.valid[,2])) {
    # Read in the unperturbed model output.
    res.temp <- paras.r11.output[[j]]
    # Calculate the sum of squared differences between the unperturbed model 
    # output and the selected set of perturbed data. 
    sse.d3.temp <- sum((res.temp$den.mat.d3 - 
                          paras.r11.output.pert.temp$den.mat.d3.pert)^2)
    # Append index, discrepancy into each row of the matrix.
    diff.mat[j, 1:2] <- c(ls.r11.valid.pert[j, 1], sse.d3.temp)
    # Optional, can be used to track the progress. 
    # print(j)
  }
  # Remove the parameter set that corresponds to the current perturbed 
  # model output
  diff.mat <- diff.mat[-ls.r11.valid.pert.sort.min[i],]
  
  # Store the info matrix in the list and name it with the corresponding index
  # in the full dataset. 
  ls.r11.diff.mat[[i]] <- diff.mat
  names(ls.r11.diff.mat)[i] <- paste0("diff_mat_r11_",
                                      ls.r11.pert.sort.min.index[i])
  # Progress tracking, optional
  print(i)
}

# Carry out another round of ABC but using these perturbed output as reference
# data. Calculate the resampling probabilities of each parameter vector based
# on the least square difference between their simulation output and the "new
# reference data". (Run in parallel). 
est <- foreach (i = 1:length(ls.r11.diff.mat), .combine = rbind) %dopar% {
  # Search for the bandwidth factor which will give the resampling weights and
  # yields the correct ESS. 
  diff.mat.temp <- ls.r11.diff.mat[[i]]
  # Initial interval for the bandwidth factors. 
  lb.bw.temp <- 6.12 # 2.53 
  ub.bw.temp <- 6.24 # 2.65 
  info.list.temp <- calculate.bw(ss.mat = diff.mat.temp, lb.bw = lb.bw.temp,
                                 ub.bw = ub.bw.temp, ess.target = 1640.25, 
                                 step.size = 0.01)
  # Change the lower and upper bounds of the interval if the current interval
  # did not yield the correct ESS. 
  while (info.list.temp$bw.obj == lb.bw.temp || info.list.temp$bw.obj == 
         ub.bw.temp) {
    if (info.list.temp$bw.obj == lb.bw.temp) {
      lb.bw.temp <- lb.bw.temp - 0.1
      ub.bw.temp <- lb.bw.temp + 0.01
      info.list.temp <- calculate.bw(ss.mat = diff.mat.temp, 
                                     lb.bw = lb.bw.temp,
                                     ub.bw = ub.bw.temp, 
                                     ess.target = 1640.25,
                                     step.size = 0.01)
    } else if (info.list.temp$bw.obj == ub.bw.temp) {
      lb.bw.temp <- ub.bw.temp - 0.01
      ub.bw.temp <- ub.bw.temp + 0.1
      info.list.temp <- calculate.bw(ss.mat = diff.mat.temp, 
                                     lb.bw = lb.bw.temp,
                                     ub.bw = ub.bw.temp, 
                                     ess.target = 1640.25,
                                     step.size = 0.01)
    }
  }
  # Save the results for each "pseudo-reference dataset" into .rds files in 
  # the correct directory. 
  readr::write_rds(info.list.temp, path = paste0("./", save.sims.dir, 
                                                 "/info_list_r11_", 
                                                 ls.r11.pert.sort.min.index[i], 
                                                 ".rds"))
  c(i, info.list.temp$bw.obj, info.list.temp$ess.obj)
}
stopCluster(cl)

# Read in the results obtained in the previous step, store them into one 
# single list. 
ls.r11.crage <- vector(mode = "list")
for (i in 1:length(ls.r11.diff.mat)) {
  info.list.temp <- read_rds(paste0("./", save.sims.dir, 
                                    "/info_list_r11_", 
                                    ls.r11.pert.sort.min.index[i], ".rds"))
  ls.r11.crage[[i]] <- info.list.temp
  names(ls.r11.crage)[i] <- paste0("info_list_r11_",
                                   ls.r11.pert.sort.min.index[i])
}

# Save the results into a .rds file
#write_rds(ls.r11.crage, "Coverage_test_disp_wt_prob_r11_pert_non_zero.rds")
write_rds(ls.r11.crage, "Coverage_test_disp_wt_prob_r11_pert.rds")

# ls.r11.crage <- read_rds("Coverage_test_disp_wt_prob_r11_pert_non_zero.rds")

# Calculation of coverage probabilities and check if they show deviations 
# to uniform distributions.

# Add indices to the parameter vectors. 
paras.r11.valid <- cbind(ls.r11.valid[,1], paras.r11.valid)

# An empty matrix used to store all the coverage probabilities
# regarding each perturbed dataset. 
cov.mat <- matrix(0, nrow = 200, ncol = (length(paras.r11.valid[1,]) - 1))

# Coverage probabilities calculations
for (i in 1:200) {
  # Read in the information matrix regarding each perturbed dataset
  info.mat.temp <- ls.r11.crage[[i]]$info.mat
  # Remove the parameter set which generated the perturbed dataset, 
  # in order to keep the indices consistent.
  paras.r11.valid.corres <- paras.r11.valid[-ls.r11.valid.pert.sort.min[i],]
  # For each parameter, calculate its coverage probability
  prob.vec <- vector()
  for (j in 1:(length(paras.r11.valid[1,]) - 1)) {
    cov.ind.temp <- which(paras.r11.valid.corres[, (j + 1)] <= 
                            paras.r11.valid[ls.r11.valid.pert.sort.min[i], (j + 1)])
    prob.temp <- sum(info.mat.temp[cov.ind.temp, 
                                   length(info.mat.temp[1,])])/
      sum(info.mat.temp[,length(info.mat.temp[1,])])
    prob.vec <- c(prob.vec, prob.temp)
  }
  # Store the coverage probability corresponding to each perturbed dataset into 
  # each row of the coverage probability matrix. 
  cov.mat[i, ] <- prob.vec
}

# Store the coverage probabilities
# write.table(cov.mat, "Coverage probabilities r11 non zero only.txt")
write.table(cov.mat, "Coverage probabilities r11.txt")

cov.mat <- read.table("Coverage probabilities r11 non zero only.txt", sep = "", header = TRUE)


# Histograms plots & uniformity check.

# Reference uniform distribution. 
set.seed(874512)
unif.sample <- runif(200, 0, 1)

hist(cov.mat[,1], main = TeX("Coverage check of $d_{n}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,1], unif.sample)
# No evidence against H0 (non-zero)
# No evidence against H0 (all)

hist(cov.mat[,2], main = TeX("Coverage check of $\\gamma$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,2], unif.sample) 
# No evidence against H0 (non-zero). 
# No evidence against H0. (all)

hist(cov.mat[,3], main = TeX("Coverage check of $r_{n}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,3], unif.sample) 
# No evidence against H0 (non-zero). 
# No evidence against H0. (all)

hist(cov.mat[,4], main = TeX("Coverage check of $\\eta$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,4], unif.sample) 
# No evidence against H0 (non-zero). 
# No evidence against H0. (all)

hist(cov.mat[,5], main = TeX("Coverage check of $d_{m}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,5], unif.sample) 
# No evidence against H0 (non-zero). 
# No evidence against H0. (all).

hist(cov.mat[,6], main = TeX("Coverage check of $\\alpha$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,6], unif.sample) 
# No evidence against H0 (non-zero).
# No evidence against H0. (all). 

hist(cov.mat[,7], main = TeX("Coverage check of $R_{init.}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,7], unif.sample) 
# No evidence against H0 (non-zero).
# Weak evidence against H0 (all).

hist(cov.mat[,8], main = TeX("Coverage check of $P_{ext.}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,8], unif.sample) 
# No evidence against H0 (non-zero).
# Weak evidence against H0 (all).

hist(cov.mat[,9], main = TeX("Coverage check of $P_{mit.}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,9], unif.sample) 
# No evidence against H0. 

# Non zero data perturbed only:
# 9/9 parameters passed the coverage test! 

# All data perturbed: 
# 8/9 parameters passed the coverage test!
# 1/9 parameters shown weak evidence against H0 of the coverage test.
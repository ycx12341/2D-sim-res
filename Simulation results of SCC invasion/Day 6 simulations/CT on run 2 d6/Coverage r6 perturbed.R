# Coverage r6 perturbed.R
# This .R file conducts the coverage test on the final parameter samples 
# obtained at the end of the simulation for post-day 6 SCC invasion pattern to
# check if the final parameter samples can be regarded as appropriate posterior
# samples from the perspective of Bayesian inference. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(readr)
library(doParallel)
library(latex2exp)

# Set the directory to store the results. 
# save.sims.dir <- "Cov_results_r6"
save.sims.dir <- "Cov_results_r6_pert_all"
save.sims <- TRUE
if (save.sims) {
  if (!dir.exists(save.sims.dir)) dir.create(save.sims.dir)
}

# Set the cluster. 
n.thread <- detectCores()/2
cl <- makeCluster(n.thread)
registerDoParallel(cl)

# Read in the final parameter values. 
paras.r6 <- read.table("Round 6 parameters.txt", sep = "", header = TRUE)

# Directory which stores the simulation output. 
sims.dir <- "LS_results_r6"

# Check if there are still NaNs in the final output.
ls.r6 <- read.table("Round 6 Least Square.txt", sep = "", header = TRUE)
length(which(is.na(ls.r6[,2]))) # 25 singular values. 

# Remove the singular simulation output and their corresponding parameter 
# vectors.
ls.r6.valid <- ls.r6[-which(is.na(ls.r6[,2])), ]
paras.r6.valid <- paras.r6[-which(is.na(ls.r6[,2])), ]

# Read in the file which contains the source functions and the observed 
# densities.  
source("PDE 2D ABC functions adjusted SCC.r")

# Read in the unperturbed model output and store them in a list. 
paras.r6.output <- vector(mode = "list")
for (i in 1:length(ls.r6.valid[,2])) {
  res.temp <- read_rds(paste0("./", sims.dir, 
                              "/Round_6_paras",
                              ls.r6.valid[i,1],"_res.rds"))
  ls.temp <- list(den.mat.d6 = res.temp$den.mat.d6)
  paras.r6.output[[i]] <- ls.temp
  print(i)
  names(paras.r6.output)[i] <- paste0("ls_r6_",ls.r6.valid[i,1])
}

# Write the unperturbed model output into a .rds file. 
write_rds(paras.r6.output, "Full output paras r6 run 2 valid.rds")

# Read in the output in round 11 that gave the minimum result of least square
# difference, and calculate the standard deviation (SD) based on it.
ls.r6.min <- read_rds("Round_6_paras5275_res.rds")
sd.r6.d6.min <- sd(ls.r6.min$den.mat.d6 - t6.ref.den)

# Add perturbations (gaussian) to the final unperturbed output based on the
# SD calculated in the previous step. 
paras.r6.output.pert <- vector(mode = "list")
set.seed(874514)
RNGkind(sample.kind = "Rejection")
for (i in 1:length(paras.r6.output)) {
  ls.temp <- paras.r6.output[[i]]
  # Read the unperturbed data
  den.mat.d6.temp.unpert <- ls.temp$den.mat.d6
  # Create empty matrices to store perturbed data
  den.mat.d6.temp.pert <- matrix(0, nrow = nrow(den.mat.d6.temp.unpert), 
                                 ncol = ncol(den.mat.d6.temp.unpert))
  # Add perturbations to the unperturbed output: (two different ways)
  # 1. Add perturbations to non-zero data only.
  # 2. Add perturbations to the full dataset. 
  # We believed the first method of adding perturbation is more reasonable, as
  # the SCC starts invading the surroundings from the left to the right. In 
  # the reference data (especially for the first few patterns), many regions 
  # close to the right boundary in the domain are still unoccupied. Hence, 
  # adding perturbation to all data points would corrupt the nature of the 
  # original dataset. 
  for (j in 1:length(den.mat.d6.temp.unpert)) {
    #if(den.mat.d6.temp.unpert[j] != 0) {
     den.mat.d6.temp.pert[j] <- rnorm(1, den.mat.d6.temp.unpert[j], sd = sd.r6.d6.min)
     while (den.mat.d6.temp.pert[j] < 0 || den.mat.d6.temp.pert[j] > 1) {
       den.mat.d6.temp.pert[j] <- rnorm(1, den.mat.d6.temp.unpert[j], sd = sd.r6.d6.min)
     }
    #}
  }
  
  # Sum of squared differences between the perturbed densities and the 
  # observed densities.
  sse.d6.temp <- sum((t6.ref.den - den.mat.d6.temp.pert)^2)
  # Store each set of perturbed data into the list
  paras.r6.output.pert[[i]] <- list(den.mat.d6.pert = den.mat.d6.temp.pert, 
                                    diff = sse.d6.temp)
  # Name the elements in the list with correct indices
  names(paras.r6.output.pert)[i] <- paste0("ls_r6_",ls.r6.valid[i,1],"_pert")
  print(i)
}

# Sort the sum of squared differences
ls.r6.valid.pert <- matrix(0, nrow = length(ls.r6.valid[,2])
                     ,ncol = 2)
for (i in 1:length(ls.r6.valid.pert[,2])) {
  ls.r6.valid.pert[i,] <- c(ls.r6.valid[i, 1], 
                            as.double(paras.r6.output.pert[[i]]$diff))
}
ls.r6.valid.pert.sort <- order(ls.r6.valid.pert[,2])
ls.r6.valid.pert.sort.min <- ls.r6.valid.pert.sort[1:200]

# Index in the full dataset (from 1 - 10000)
ls.r6.pert.sort.min.index <- ls.r6.valid.pert[ls.r6.valid.pert.sort[1:200], 1]

# Locate the corresponding 200 perturbed model outputs that have the 
# minimum discrepancy with the reference data.
paras.r6.output.pert.full <- paras.r6.output.pert
paras.r6.output.pert <- paras.r6.output.pert[ls.r6.valid.pert.sort.min]

# Discrepancies between unperturbed model output with each of these 200 
# perturbed model output which were treated as new reference data. 
ls.r6.diff.mat <- vector(mode = "list")
for (i in 1:length(paras.r6.output.pert)) {
  # For each perturbed dataset, record its discrepancy with other unperturbed 
  # model output, corresponding weights and resampling probabilities. 
  diff.mat <- matrix(0, nrow = length(paras.r6.valid[,2]), ncol = 2)
  paras.r6.output.pert.temp <- paras.r6.output.pert[[i]]
  for (j in 1:length(paras.r6.valid[,2])) {
    # Read in the unperturbed model output.
    res.temp <- paras.r6.output[[j]]
    # Calculate the sum of squared differences between the unperturbed model 
    # output and the selected set of perturbed data. 
    sse.d6.temp <- sum((res.temp$den.mat.d6 - 
                          paras.r6.output.pert.temp$den.mat.d6.pert)^2)
    # Append index, discrepancy into each row of the matrix
    diff.mat[j, 1:2] <- c(ls.r6.valid.pert[j, 1], sse.d6.temp)
    # Optional, can be used to track the progress. 
    # print(j)
  }
  # Remove the parameter set that corresponds to the current perturbed 
  # model output
  diff.mat <- diff.mat[-ls.r6.valid.pert.sort.min[i],]
  
  # Store the info matrix in the list and name it with the corresponding index
  # in the full dataset. 
  ls.r6.diff.mat[[i]] <- diff.mat
  names(ls.r6.diff.mat)[i] <- paste0("diff_mat_r6_",
                                     ls.r6.pert.sort.min.index[i])
  # Progress tracking, optional
  print(i)
}

# Carry out another round of ABC but using these perturbed output as reference
# data. Calculate the resampling probabilities of each parameter vector based
# on the least square difference between their simulation output and the "new
# reference data". (Run in parallel). 
est <- foreach (i = 1:length(ls.r6.diff.mat), .combine = rbind) %dopar% {
  # Search for the bandwidth factor which will give the resampling weights and
  # yields the correct ESS. 
  diff.mat.temp <- ls.r6.diff.mat[[i]]
  # Initial interval for the bandwidth factors. 
  lb.bw.temp <- 6.42 #3.57
  ub.bw.temp <- 6.54 #3.69 
  info.list.temp <- calculate.bw(ss.mat = diff.mat.temp, lb.bw = lb.bw.temp,
                                 ub.bw = ub.bw.temp, ess.target = 2500, 
                                 step.size = 0.01)
  # Change the lower and upper bounds of the interval if the current interval
  # did not yield the correct ESS. 
  while (info.list.temp$bw.obj == lb.bw.temp || info.list.temp$bw.obj == 
         ub.bw.temp) {
    if (info.list.temp$bw.obj == lb.bw.temp) {
      lb.bw.temp <- lb.bw.temp - 0.1
      ub.bw.temp <- lb.bw.temp + 0.01
      info.list.temp <- calculate.bw(ss.mat = diff.mat.temp, lb.bw = lb.bw.temp,
                                     ub.bw = ub.bw.temp, ess.target = 2500,
                                     step.size = 0.01)
    } else if (info.list.temp$bw.obj == ub.bw.temp) {
      lb.bw.temp <- ub.bw.temp - 0.01
      ub.bw.temp <- ub.bw.temp + 0.1
      info.list.temp <- calculate.bw(ss.mat = diff.mat.temp, lb.bw = lb.bw.temp,
                                     ub.bw = ub.bw.temp, ess.target = 2500,
                                     step.size = 0.01)
    }
  }
  # Save the results for each "pseudo-reference dataset" into .rds files in 
  # the correct directory. 
  readr::write_rds(info.list.temp, path = paste0("./", save.sims.dir, 
                                                 "/info_list_r6_", 
                                                 ls.r6.pert.sort.min.index[i], 
                                                 ".rds"))
  c(i, info.list.temp$bw.obj, info.list.temp$ess.obj)
}
stopCluster(cl)

# Read in the results obtained in the previous step, store them into one 
# single list. 
ls.r6.crage <- vector(mode = "list")
for (i in 1:length(ls.r6.diff.mat)) {
  info.list.temp <- read_rds(paste0("./", save.sims.dir, 
                                    "/info_list_r6_", 
                                    ls.r6.pert.sort.min.index[i], ".rds"))
  ls.r6.crage[[i]] <- info.list.temp
  names(ls.r6.crage)[i] <- paste0("info_list_r6_",ls.r6.pert.sort.min.index[i])
}

# Save the results into a .rds file
# write_rds(ls.r6.crage, "Coverage_test_disp_wt_prob_r6_pert_non_zero.rds")
write_rds(ls.r6.crage, "Coverage_test_disp_wt_prob_r6_pert.rds")

# ls.r6.crage <- read_rds("Coverage_test_disp_wt_prob_r6_pert_
# sep_non_zero.rds")

# Calculation of coverage probabilities and check if they show deviations 
# to uniform distributions.

# Add indices to the parameter vectors. 
paras.r6.valid <- cbind(ls.r6.valid[,1], paras.r6.valid)

# An empty matrix used to store all the coverage probabilities
# regarding each perturbed dataset. 
cov.mat <- matrix(0, nrow = 200, ncol = (length(paras.r6.valid[1,]) - 1))

# Coverage probabilities calculations
for (i in 1:200) {
  # Read in the information matrix regarding each perturbed dataset
  info.mat.temp <- ls.r6.crage[[i]]$info.mat
  # Remove the parameter set which generated the perturbed dataset, 
  # in order to keep the indices consistent.
  paras.r6.valid.corres <- paras.r6.valid[-ls.r6.valid.pert.sort.min[i],]
  # For each parameter, calculate its coverage probability
  prob.vec <- vector()
  for (j in 1:(length(paras.r6.valid[1,]) - 1)) {
    cov.ind.temp <- which(paras.r6.valid.corres[, (j + 1)] <= 
                            paras.r6.valid[ls.r6.valid.pert.sort.min[i], (j + 1)])
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
# write.table(cov.mat, "Coverage probabilities r6 non zero only.txt")
write.table(cov.mat, "Coverage probabilities r6.txt")

cov.mat <- read.table("Coverage probabilities r6 non zero only.txt", sep = "", header = TRUE)


# Histograms plots & uniformity check.

# Reference uniform distribution. 
set.seed(874514)
unif.sample <- runif(200, 0, 1)

hist(cov.mat[,1], main = TeX("Coverage check of $d_{n}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,1], unif.sample)
# No evidence against H0 (non-zero)
# Weak evidence against H0 (all)

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

hist(cov.mat[,4], main = TeX("Coverage check of $\\eta'$s final sample",
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

hist(cov.mat[,7], main = TeX("Coverage check of $P_{ext.}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,7], unif.sample) 
# No evidence against H0 (non-zero).
# Weak evidence against H0 (all).

hist(cov.mat[,8], main = TeX("Coverage check of $P_{mit.}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,8], unif.sample) 
# No evidence against H0 (non-zero).
# No evidence against H0 (all).

# Non zero data perturbed only:
# 8/8 parameters passed the coverage test! 

# All data perturbed: 
# 7/8 parameters passed the coverage test!
# 1/8 parameters shown weak evidence against H0 of the coverage test.
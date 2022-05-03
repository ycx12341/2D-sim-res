# Coverage r9 perturbed.R
# This .R file conducts the coverage test on the final parameter samples 
# obtained at the end of the time-dependent simulation for SCC invasion patterns 
# to check if the final parameter samples can be regarded as appropriate 
# posterior samples from the perspective of Bayesian inference. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(readr)
library(doParallel)
library(latex2exp)

# Set the directory to store the results. 
save.sims.dir <- "Cov_results_r9"
# save.sims.dir <- "Cov_results_r9_pert_all"
save.sims <- TRUE
if (save.sims) {
  if (!dir.exists(save.sims.dir)) dir.create(save.sims.dir)
}

# Set the cluster. 
n.thread <- detectCores()/2
cl <- makeCluster(n.thread)
registerDoParallel(cl)

# Read in the final parameter values. 
paras.r9 <- read.table("Round 9 parameters.txt", sep = "", header = TRUE)

# Directory which stores the simulation output. 
sims.dir <- "LS_results_r9"

# Check if there are still NaNs in the final output.
ls.r9 <- read.table("Round 9 Least Square.txt", 
                     sep = "", header = TRUE)
length(which(is.na(ls.r9[,2]))) # 740 singular values.  

# Remove the singular simulation output and their corresponding parameter 
# vectors. 
ls.r9.valid <- ls.r9[-which(is.na(ls.r9[,2])), ]
paras.r9.valid <- paras.r9[-which(is.na(ls.r9[,2])), ]

# Read in the file which contains the source functions and the observed 
# densities.  
source("PDE 2D ABC functions adjusted SCC (time-dependent parameters).r")

# Read in the unperturbed model output and store them in a list. 
paras.r9.output <- vector(mode = "list")
for (i in 1:length(ls.r9.valid[,2])) {
  res.temp <- read_rds(paste0("./", sims.dir, 
                              "/Round_9_paras",
                              ls.r9.valid[i,1],"_res.rds"))
  ls.temp <- list(den.mat.d3 = res.temp$den.mat.d3,
                  den.mat.d6 = res.temp$den.mat.d6,
                  den.mat.d9 = res.temp$den.mat.d9,
                  den.mat.d12 = res.temp$den.mat.d12,
                  den.mat.d14 = res.temp$den.mat.d14)
  paras.r9.output[[i]] <- ls.temp
  print(i)
  names(paras.r9.output)[i] <- paste0("ls_r9_",ls.r9.valid[i,1])
}

# Write the unperturbed model output into a .rds file. 
write_rds(paras.r9.output, "Full output paras r9 run 3 valid.rds")

# Read in the output in round 8 that gave the minimum result of least square
# difference, and calculate the standard deviation (SD) based on it.
ls.r9.min <- read_rds("Round_9_paras15251_res.rds")
sd.r9.d3.min <- sd(ls.r9.min$den.mat.d3 - t3.ref.den)
sd.r9.d6.min <- sd(ls.r9.min$den.mat.d6 - t6.ref.den)
sd.r9.d9.min <- sd(ls.r9.min$den.mat.d9 - t9.ref.den)
sd.r9.d12.min <- sd(ls.r9.min$den.mat.d12 - t12.ref.den)
sd.r9.d14.min <- sd(ls.r9.min$den.mat.d14 - t14.ref.den)

# Add perturbations (gaussian) to the final unperturbed output based on the
# SD calculated in the previous step. 
paras.r9.output.pert <- vector(mode = "list")
set.seed(874514)
RNGkind(sample.kind = "Rejection")
for (i in 1:length(paras.r9.output)) {
  ls.temp <- paras.r9.output[[i]]
  # Read the unperturbed data
  den.mat.d3.temp.unpert <- ls.temp$den.mat.d3
  den.mat.d6.temp.unpert <- ls.temp$den.mat.d6
  den.mat.d9.temp.unpert <- ls.temp$den.mat.d9
  den.mat.d12.temp.unpert <- ls.temp$den.mat.d12
  den.mat.d14.temp.unpert <- ls.temp$den.mat.d14
  
  # Create empty matrices to store perturbed data
  den.mat.d3.temp.pert <- matrix(0, nrow = nrow(den.mat.d3.temp.unpert), 
                                 ncol = ncol(den.mat.d3.temp.unpert))
  den.mat.d6.temp.pert <- matrix(0, nrow = nrow(den.mat.d6.temp.unpert), 
                                 ncol = ncol(den.mat.d6.temp.unpert))
  den.mat.d9.temp.pert <- matrix(0, nrow = nrow(den.mat.d9.temp.unpert), 
                                 ncol = ncol(den.mat.d9.temp.unpert))
  den.mat.d12.temp.pert <- matrix(0, nrow = nrow(den.mat.d12.temp.unpert), 
                                 ncol = ncol(den.mat.d12.temp.unpert))
  den.mat.d14.temp.pert <- matrix(0, nrow = nrow(den.mat.d14.temp.unpert), 
                                 ncol = ncol(den.mat.d14.temp.unpert))
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
   if(den.mat.d3.temp.unpert[j] != 0) {
    den.mat.d3.temp.pert[j] <- rnorm(1, den.mat.d3.temp.unpert[j], sd = sd.r9.d3.min)
    while (den.mat.d3.temp.pert[j] < 0 || den.mat.d3.temp.pert[j] > 1) {
      den.mat.d3.temp.pert[j] <- rnorm(1, den.mat.d3.temp.unpert[j], sd = sd.r9.d3.min)
    }
   }
    
   if(den.mat.d6.temp.unpert[j] != 0) {
    den.mat.d6.temp.pert[j] <- rnorm(1, den.mat.d6.temp.unpert[j], sd = sd.r9.d6.min)
    while (den.mat.d6.temp.pert[j] < 0 || den.mat.d6.temp.pert[j] > 1) {
      den.mat.d6.temp.pert[j] <- rnorm(1, den.mat.d6.temp.unpert[j], sd = sd.r9.d6.min)
    }
   }
    
   if(den.mat.d9.temp.unpert[j] != 0) {
    den.mat.d9.temp.pert[j] <- rnorm(1, den.mat.d9.temp.unpert[j], sd = sd.r9.d9.min)
    while (den.mat.d9.temp.pert[j] < 0 || den.mat.d9.temp.pert[j] > 1) {
      den.mat.d9.temp.pert[j] <- rnorm(1, den.mat.d9.temp.unpert[j], sd = sd.r9.d9.min)
    }
   }
    
   if(den.mat.d12.temp.unpert[j] != 0) {
    den.mat.d12.temp.pert[j] <- rnorm(1, den.mat.d12.temp.unpert[j], sd = sd.r9.d12.min)
    while (den.mat.d12.temp.pert[j] < 0 || den.mat.d12.temp.pert[j] > 1) {
      den.mat.d12.temp.pert[j] <- rnorm(1, den.mat.d12.temp.unpert[j], sd = sd.r9.d12.min)
    }
   }
    
   if(den.mat.d14.temp.unpert[j] != 0) {
    den.mat.d14.temp.pert[j] <- rnorm(1, den.mat.d14.temp.unpert[j], sd = sd.r9.d14.min)
    while (den.mat.d14.temp.pert[j] < 0 || den.mat.d14.temp.pert[j] > 1) {
      den.mat.d14.temp.pert[j] <- rnorm(1, den.mat.d14.temp.unpert[j], sd = sd.r9.d14.min)
    }
   }
  }
  
  # Sum of squared differences between the perturbed densities and the 
  # observed densities.
  sse.d3.temp <- sum((t3.ref.den - den.mat.d3.temp.pert)^2)
  sse.d6.temp <- sum((t6.ref.den - den.mat.d6.temp.pert)^2)
  sse.d9.temp <- sum((t9.ref.den - den.mat.d9.temp.pert)^2)
  sse.d12.temp <- sum((t12.ref.den - den.mat.d12.temp.pert)^2)
  sse.d14.temp <- sum((t14.ref.den - den.mat.d14.temp.pert)^2)
  diff.temp.pert <- sum(sse.d3.temp, sse.d6.temp, sse.d9.temp,
                        sse.d12.temp, sse.d14.temp)
  # Store each set of perturbed data into the list
  paras.r9.output.pert[[i]] <- list(den.mat.d3.pert = den.mat.d3.temp.pert, 
                                    den.mat.d6.pert = den.mat.d6.temp.pert,
                                    den.mat.d9.pert = den.mat.d9.temp.pert,
                                    den.mat.d12.pert = den.mat.d12.temp.pert,
                                    den.mat.d14.pert = den.mat.d14.temp.pert,
                                    diff = diff.temp.pert)
  # Name the elements in the list with correct indices
  names(paras.r9.output.pert)[i] <- paste0("ls_r9_",ls.r9.valid[i,1],"_pert")
  print(i)
}

# Sort the sum of squared differences
ls.r9.valid.pert <- matrix(0, nrow = length(ls.r9.valid[,2])
                            ,ncol = 2)
for (i in 1:length(ls.r9.valid.pert[,2])) {
  ls.r9.valid.pert[i,] <- c(ls.r9.valid[i, 1], 
                             as.double(paras.r9.output.pert[[i]]$diff))
}
ls.r9.valid.pert.sort <- order(ls.r9.valid.pert[,2])
ls.r9.valid.pert.sort.min <- ls.r9.valid.pert.sort[1:200]

# Index in the full dataset (from 1 - 10000)
ls.r9.pert.sort.min.index <- ls.r9.valid.pert[ls.r9.valid.pert.sort[1:200], 
                                                1]

# Locate the corresponding 200 perturbed model outputs that have the 
# minimum discrepancy with the reference data.
paras.r9.output.pert.full <- paras.r9.output.pert
paras.r9.output.pert <- paras.r9.output.pert[ls.r9.valid.pert.sort.min]

# Discrepancies between unperturbed model output with each of these 200 
# perturbed model output which were treated as new reference data. 
ls.r9.diff.mat <- vector(mode = "list")
for (i in 1:length(paras.r9.output.pert)) {
  # For each perturbed dataset, record its discrepancy with other unperturbed 
  # model output, corresponding weights and resampling probabilities. 
  diff.mat <- matrix(0, nrow = length(paras.r9.valid[,2]), ncol = 2)
  paras.r9.output.pert.temp <- paras.r9.output.pert[[i]]
  for (j in 1:length(paras.r9.valid[,2])) {
    # Read in the unperturbed model output.
    res.temp <- paras.r9.output[[j]]
    # Calculate the sum of squared differences between the unperturbed model 
    # output and the selected set of perturbed data. 
    sse.d3.temp <- sum((res.temp$den.mat.d3 - 
                          paras.r9.output.pert.temp$den.mat.d3.pert)^2)
    sse.d6.temp <- sum((res.temp$den.mat.d6 - 
                          paras.r9.output.pert.temp$den.mat.d6.pert)^2)
    sse.d9.temp <- sum((res.temp$den.mat.d9 - 
                          paras.r9.output.pert.temp$den.mat.d9.pert)^2)
    sse.d12.temp <- sum((res.temp$den.mat.d12 - 
                          paras.r9.output.pert.temp$den.mat.d12.pert)^2)
    sse.d14.temp <- sum((res.temp$den.mat.d14 - 
                          paras.r9.output.pert.temp$den.mat.d14.pert)^2)
    diff.temp <- sum(sse.d3.temp, sse.d6.temp, sse.d9.temp, sse.d12.temp,
                     sse.d14.temp)
    # Append index, discrepancy into each row of the matrix.
    diff.mat[j, 1:2] <- c(ls.r9.valid.pert[j, 1], diff.temp)
    # Optional, can be used to track the progress. 
    # print(j)
  }
  # Remove the parameter set that corresponds to the current perturbed 
  # model output
  diff.mat <- diff.mat[-ls.r9.valid.pert.sort.min[i],]
  
  # Store the info matrix in the list and name it with the corresponding index
  # in the full dataset. 
  ls.r9.diff.mat[[i]] <- diff.mat
  names(ls.r9.diff.mat)[i] <- paste0("diff_mat_r9_",
                                      ls.r9.pert.sort.min.index[i])
  # Progress tracking, optional
  print(i)
}

# Carry out another round of ABC but using these perturbed output as reference
# data. Calculate the resampling probabilities of each parameter vector based
# on the least square difference between their simulation output and the "new
# reference data". (Run in parallel). 
est <- foreach (i = 1:length(ls.r9.diff.mat), .combine = rbind) %dopar% {
  # Search for the bandwidth factor which will give the resampling weights and
  # yields the correct ESS. 
  diff.mat.temp <- ls.r9.diff.mat[[i]]
  # Initial interval for the bandwidth factors. 
  lb.bw.temp <- 15.85 # 10.05  
  ub.bw.temp <- 15.95 # 10.15 
  info.list.temp <- calculate.bw(ss.mat = diff.mat.temp, lb.bw = lb.bw.temp,
                                 ub.bw = ub.bw.temp, ess.target = 10125, 
                                 step.size = 0.01)
  bw.obj.curr <- unname(info.list.temp$bw.obj)
  # Change the lower and upper bounds of the interval if the current interval
  # did not yield the correct ESS. 
  while (isTRUE(all.equal(bw.obj.curr, lb.bw.temp)) || isTRUE(all.equal(bw.obj.curr, ub.bw.temp))) {
    if (isTRUE(all.equal(bw.obj.curr, lb.bw.temp))) {
      ub.bw.temp <- lb.bw.temp + 0.01
      lb.bw.temp <- lb.bw.temp - 0.1
      info.list.temp <- calculate.bw(ss.mat = diff.mat.temp, 
                                     lb.bw = lb.bw.temp,
                                     ub.bw = ub.bw.temp, 
                                     ess.target = 10125,
                                     step.size = 0.01)
      bw.obj.curr <- unname(info.list.temp$bw.obj)
    } else if (isTRUE(all.equal(bw.obj.curr, ub.bw.temp))) {
      lb.bw.temp <- ub.bw.temp - 0.01
      ub.bw.temp <- ub.bw.temp + 0.1
      info.list.temp <- calculate.bw(ss.mat = diff.mat.temp, 
                                     lb.bw = lb.bw.temp,
                                     ub.bw = ub.bw.temp, 
                                     ess.target = 10125,
                                     step.size = 0.01)
      bw.obj.curr <- unname(info.list.temp$bw.obj)
    }
  }
  # Save the results for each "pseudo-reference dataset" into .rds files in 
  # the correct directory. 
  readr::write_rds(info.list.temp, path = paste0("./", save.sims.dir, 
                                                 "/info_list_r9_", 
                                                 ls.r9.pert.sort.min.index[i], 
                                                 ".rds"))
  c(i, info.list.temp$bw.obj, info.list.temp$ess.obj)
}
stopCluster(cl)

# Read in the results obtained in the previous step, store them into one 
# single list. 
ls.r9.crage <- vector(mode = "list")
for (i in 1:length(ls.r9.diff.mat)) {
  info.list.temp <- read_rds(paste0("./", save.sims.dir, 
                                    "/info_list_r9_", 
                                    ls.r9.pert.sort.min.index[i], ".rds"))
  ls.r9.crage[[i]] <- info.list.temp
  names(ls.r9.crage)[i] <- paste0("info_list_r9_",
                                   ls.r9.pert.sort.min.index[i])
}

# Save the results into a .rds file
write_rds(ls.r9.crage, "Coverage_test_disp_wt_prob_r9_pert_non_zero.rds")
#write_rds(ls.r9.crage, "Coverage_test_disp_wt_prob_r9_pert_all.rds")

# Calculation of coverage probabilities and check if they show deviations 
# to uniform distributions.

# Add indices to the parameter vectors. 
paras.r9.valid <- cbind(ls.r9.valid[,1], paras.r9.valid)

# Generate the values of the time-dependent parameters at later periods based
# on the posterior first-period values.
dn.varying.paras.full.mat <- dn.varying.paras(dn.init.par = paras.r9.valid[,2],
                                              dn.quad.temp = paras.r9.valid[,11],
                                              dn.lin.temp = paras.r9.valid[,12])
dn.varying.paras.25 <- dn.varying.paras.full.mat[,2:5]

gamma.varying.paras.full.mat <- gamma.varying.paras(gamma.init.par = paras.r9.valid[,3],
                                              gamma.quad.temp = paras.r9.valid[,13],
                                              gamma.lin.temp = paras.r9.valid[,14])
gamma.varying.paras.25 <- gamma.varying.paras.full.mat[,2:5]

rn.varying.paras.full.mat <- rn.varying.paras(rn.init.par = paras.r9.valid[,4],
                                                    rn.quad.temp = paras.r9.valid[,15],
                                                    rn.lin.temp = paras.r9.valid[,16])
rn.varying.paras.25 <- rn.varying.paras.full.mat[,2:5]

eta.varying.paras.full.mat <- eta.varying.paras(eta.init.par = paras.r9.valid[,5],
                                              eta.quad.temp = paras.r9.valid[,17],
                                              eta.lin.temp = paras.r9.valid[,18])
eta.varying.paras.25 <- eta.varying.paras.full.mat[,2:5]

alpha.varying.paras.full.mat <- alpha.varying.paras(alpha.init.par = paras.r9.valid[,7],
                                                alpha.quad.temp = paras.r9.valid[,19],
                                                alpha.lin.temp = paras.r9.valid[,20])
alpha.varying.paras.25 <- alpha.varying.paras.full.mat[,2:5]

prob.prof.varying.paras.full.mat <- prob.prof.varying.paras(prob.prof.init.par = paras.r9.valid[,10],
                                                prob.prof.quad.temp = paras.r9.valid[,21],
                                                prob.prof.lin.temp = paras.r9.valid[,22])
prob.prof.varying.paras.25 <- prob.prof.varying.paras.full.mat[,2:5]

# Combine these time-dependent parameter values with the original parameter
# matrix to form the full parameter matrix. 
paras.r9.valid.full <- cbind(paras.r9.valid, dn.varying.paras.25,
                             gamma.varying.paras.25, rn.varying.paras.25,
                             eta.varying.paras.25, alpha.varying.paras.25,
                             prob.prof.varying.paras.25)

colnames(paras.r9.valid.full) <- c("Index", "dn", "gamma", "rn", "eta", "dm",
                                   "alpha", "R.init","P.ext", "P.mit", # 1 - 10
                                   "dn.quad", "dn.lin", "gamma.quad", "gamma.lin",
                                   "rn.quad", "rn.lin", "eta.quad", "eta.lin",
                                   "alpha.quad", "alpha.lin", 
                                   "P.mit.quad", "P.mit.lin", # 11 - 22
                                   "dn.p2", "dn.p3", "dn.p4", "dn.p5",
                                   "gamma.p2", "gamma.p3", "gamma.p4", "gamma.p5",
                                   "rn.p2", "rn.p3", "rn.p4", "rn.p5",
                                   "eta.p2", "eta.p3", "eta.p4", "eta.p5",
                                   "alpha.p2", "alpha.p3", "alpha.p4", "alpha.p5",
                                   "P.mit.p2", "P.mit.p3", "P.mit.p4", "P.mit.p5") # 23 - 46


# An empty matrix used to store all the coverage probabilities
# regarding each perturbed dataset. 
cov.mat <- matrix(0, nrow = 200, ncol = (length(paras.r9.valid.full[1,]) - 1))

# Coverage probabilities calculations
for (i in 1:200) {
  # Read in the information matrix regarding each perturbed dataset
  info.mat.temp <- ls.r9.crage[[i]]$info.mat
  # Remove the parameter set which generated the perturbed dataset, 
  # in order to keep the indices consistent.
  paras.r9.valid.full.corres <- paras.r9.valid.full[-ls.r9.valid.pert.sort.min[i],]
  # For each parameter, calculate its coverage probability
  prob.vec <- vector()
  for (j in 1:(length(paras.r9.valid.full[1,]) - 1)) {
    cov.ind.temp <- which(paras.r9.valid.full.corres[, (j + 1)] <= 
                            paras.r9.valid.full[ls.r9.valid.pert.sort.min[i], (j + 1)])
    prob.temp <- sum(info.mat.temp[cov.ind.temp, 
                                    length(info.mat.temp[1,])])/
      sum(info.mat.temp[,length(info.mat.temp[1,])])
    prob.vec <- c(prob.vec, prob.temp)
  }
  # Store the coverage probability corresponding to each perturbed dataset into 
  # each row of the coverage probability matrix. 
  cov.mat[i, ] <- prob.vec
  print(i)
}

# Store the coverage probabilities
write.table(cov.mat, "Coverage probabilities r9 non zero only.txt")
# write.table(cov.mat, "Coverage probabilities r9 pert all.txt")

# Uniformity check
# Reference uniform distribution. 
set.seed(874514)
unif.sample <- runif(200, 0, 1)

colnames(cov.mat) <- c("dn", "gamma", "rn", "eta", "dm",
                     "alpha", "R.init","P.ext", "P.mit", # 1 - 9
                     "dn.quad", "dn.lin", "gamma.quad", "gamma.lin",
                     "rn.quad", "rn.lin", "eta.quad", "eta.lin",
                     "alpha.quad", "alpha.lin", 
                     "P.mit.quad", "P.mit.lin", # 10 - 21
                     "dn.p2", "dn.p3", "dn.p4", "dn.p5",
                     "gamma.p2", "gamma.p3", "gamma.p4", "gamma.p5",
                     "rn.p2", "rn.p3", "rn.p4", "rn.p5",
                     "eta.p2", "eta.p3", "eta.p4", "eta.p5",
                     "alpha.p2", "alpha.p3", "alpha.p4", "alpha.p5",
                     "P.mit.p2", "P.mit.p3", "P.mit.p4", "P.mit.p5") # 22 - 45

col.name.cov.mat <- c("dn", "gamma", "rn", "eta", "dm",
                      "alpha", "R.init","P.ext", "P.mit", # 1 - 10
                      "dn.quad", "dn.lin", "gamma.quad", "gamma.lin",
                      "rn.quad", "rn.lin", "eta.quad", "eta.lin",
                      "alpha.quad", "alpha.lin", 
                      "P.mit.quad", "P.mit.lin", # 11 - 22
                      "dn.p2", "dn.p3", "dn.p4", "dn.p5",
                      "gamma.p2", "gamma.p3", "gamma.p4", "gamma.p5",
                      "rn.p2", "rn.p3", "rn.p4", "rn.p5",
                      "eta.p2", "eta.p3", "eta.p4", "eta.p5",
                      "alpha.p2", "alpha.p3", "alpha.p4", "alpha.p5",
                      "P.mit.p2", "P.mit.p3", "P.mit.p4", "P.mit.p5")

ks.test.results.mat <- matrix(0, nrow = length(cov.mat[1,]), ncol = 2)
for (i in 1:length(cov.mat[1,])) {
  ks.test.results.mat[i,] <- c(col.name.cov.mat[i], 
                               as.double(ks.test(cov.mat[, i], 
                                                  unif.sample)$p.value))
}
p.values <- as.double(ks.test.results.mat[,2])

# Final parameter samples showing no evidence against the H0 of the coverage 
# test: 
which(p.values > 0.1) 
# 37 final parameter samples show no evidence against the H0 of the 
# coverage test. (non-zero)
# 36 final parameter samples show no evidence against the H0 of the 
# coverage test. (pert all)

which(p.values < 0.1 & p.values > 0.05) 
# 2 of the final parameter samples show weak evidence against H0. (non-zero)
# None of the final parameter samples show weak evidence against H0. (pert all)

which(p.values < 0.05 & p.values > 0.01) 
# None of the final parameter samples show moderate evidence against H0. (non-zero)
# 2 final parameter samples show moderate evidence against H0. (pert all)

which(p.values < 0.01) 
# 6 final parameter samples show strong evidence against H0. (non-zero)
# 7 final parameter samples show strong evidence against H0. (pert all)

write.table(ks.test.results.mat, "Coverage test p values non zero only.txt")
# write.table(ks.test.results.mat, "Coverage test p values pert all.txt")

# Histogram plots
hist(cov.mat[,1], main = TeX("Coverage check $d_{n}$'s final sample",
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
# No evidence against H0 (all).

hist(cov.mat[,8], main = TeX("Coverage check of $P_{ext.}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,8], unif.sample) 
# No evidence against H0 (non-zero).
# No evidence against H0 (all).

hist(cov.mat[,9], main = TeX("Coverage check of $P_{mit.}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,9], unif.sample) 
# Strong evidence against H0. (non-zero)
# Moderate evidence against H0 (all)

hist(cov.mat[,10], main = TeX("Coverage check of $d_{n, quad}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,10], unif.sample)
# No evidence against H0. (non-zero)
# No evidence against H0. (all)

hist(cov.mat[,11], main = TeX("Coverage check of $d_{n, lin}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,11], unif.sample)
# Weak evidence against H0. (non-zero)
# No evidence against H0. (all)

hist(cov.mat[,12], main = TeX("Coverage check of $\\gamma_{quad}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,12], unif.sample)
# No evidence against H0. (non-zero)
# No evidence against H0. (all)

hist(cov.mat[,13], main = TeX("Coverage check of $\\gamma_{lin}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,13], unif.sample)
# Weak evidence against H0. (non-zero)
# No evidence against H0. (all)

hist(cov.mat[,14], main = TeX("Coverage check of $\\r_{n, quad}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,14], unif.sample)
# No evidence against H0. (non-zero)
# No evidence against H0. (all)

hist(cov.mat[,15], main = TeX("Coverage check of $\\r_{n, lin}$'s final sample",
                              bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,15], unif.sample)
# No evidence against H0. (non-zero)
# No evidence against H0. (all)

hist(cov.mat[,16], main = TeX("Coverage check of $\\eta_{, quad}$'s final sample",
                              bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,16], unif.sample)
# No evidence against H0. (non-zero)
# No evidence against H0. (all)

hist(cov.mat[,17], main = TeX("Coverage check of $\\eta_{, lin}$'s final sample",
                              bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,17], unif.sample)
# No evidence against H0. (non-zero)
# No evidence against H0. (all)

hist(cov.mat[,18], main = TeX("Coverage check of $\\alpha_{ quad.}$'s final sample",
                              bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,18], unif.sample)
# No evidence against H0. (non-zero) 
# No evidence against H0. (all)

hist(cov.mat[,19], main = TeX("Coverage check of $\\alpha_{ lin.}$'s final sample",
                              bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,19], unif.sample)
# No evidence against H0. (non-zero)
# No evidence against H0. (all)

hist(cov.mat[,20], main = TeX("Coverage check of $\\P_{mit, quad.}$'s final sample",
                              bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,20], unif.sample)
# No evidence against H0. (non-zero)
# No evidence against H0. (all)

hist(cov.mat[,21], main = TeX("Coverage check of $\\P_{mit, lin.}$'s final sample",
                              bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,21], unif.sample)
# No evidence against H0. (non-zero)
# No evidence against H0. (all)

# Non zero data perturbed only:
# 18/21 parameters passed the coverage test!
# 2/21 parameters' final samples showed weak evidence against the H0 of the coverage test!
# 1/21 parameter's final sample showed strong evidence against the H0 of the coverage test! 

# All data perturbed: 
# 20/21 parameters passed the coverage test!
# 1/21 parameters showed moderate evidence against H0 of the coverage test.


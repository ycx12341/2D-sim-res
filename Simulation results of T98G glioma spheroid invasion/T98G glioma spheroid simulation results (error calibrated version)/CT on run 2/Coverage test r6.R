# Coverage test r6.R
# This .R file conducts the coverage test (an ABC posterior diagnostic that 
# checks if the final parameter samples can be regarded as appropriate samples
# from the perspective of Bayesian inference) on the final parameter samples 
# obtained in run 2 (error-calibrated version).  

# Set the workspace, then load the necessary packages. 
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

# Load the file which contains the source function.
source("PDE 2D ABC functions adjusted.R")

# Read in the final parameter values. 
paras.r6 <- read.table("Round 6 parameters log transform.txt", 
                       sep = "", header = TRUE)

# Set the directory to the correct one. 
sims.dir <- "LS_results_r6"

# Check if there are still NaNs in the output from final sample
ls.r6 <- read.table("Round 6 Least Square.txt", sep = "", header = TRUE)
length(which(is.na(ls.r6[,2]))) 
# 0, all output from the final sample are valid. 

# Observed data (glioma cell density)
ref.den <- read_rds("Glioma simulation full reference densities.rds")

# Read in the unperturbed model output and store them in a list. 
paras.r6.output <- vector(mode = "list")
for (i in 1:length(ls.r6[,2])) {
  res.temp <- read_rds(paste0("./", sims.dir, 
                              "/Round_6_paras",i,"_res.rds"))
  ls.temp.d1 <- res.temp$den.mat.d1
  ls.temp.d3 <- res.temp$den.mat.d3
  ls.temp <- list(den.mat.d1 = ls.temp.d1,
                  den.mat.d3 = ls.temp.d3)
  paras.r6.output[[i]] <- ls.temp
  print(i)
  names(paras.r6.output)[i] <- paste0("ls_r6_",i)
}

# Write the unperturbed model output into a .rds file. 
write_rds(paras.r6.output, "Full output paras r6 run 2.rds")

# Read in the output in round 7 that gave the minimum result of least square
# difference, and calculate the standard deviation (SD) based on it.
ls.r6.min <- read_rds("Round_6_paras2518_res.rds")
sd.r6.d1.min <- sd(ls.r6.min$den.mat.d1 - ref.den$t1.ref.den)
sd.r6.d3.min <- sd(ls.r6.min$den.mat.d3 - ref.den$t3.ref.den)

# Add perturbations (gaussian) to the final unperturbed output based on the
# SD calculated in the previous step. 
paras.r6.output.pert <- vector(mode = "list")

set.seed(874512)
RNGkind(sample.kind = "Rejection")
for (i in 1:length(paras.r6.output)) {
  ls.temp <- paras.r6.output[[i]]
  # Read the unperturbed data
  den.mat.d1.temp.unpert <- ls.temp$den.mat.d1
  den.mat.d3.temp.unpert <- ls.temp$den.mat.d3
  # Create empty matrices to store perturbed data
  den.mat.d1.temp.pert <- matrix(0, nrow = nrow(den.mat.d1.temp.unpert), 
                                 ncol = ncol(den.mat.d1.temp.unpert))
  den.mat.d3.temp.pert <- matrix(0, nrow = nrow(den.mat.d3.temp.unpert),
                                 ncol = ncol(den.mat.d3.temp.unpert))
  # Add perturbations to the unperturbed output: (two different ways)
  # 1. Add perturbations to non-zero data only.
  # 2. Add perturbations to the full dataset. 
  # We believed the first method of adding perturbation is more reasonable, as
  # the glioma spheroid starts invading the surroundings from the centre. In 
  # the reference data, many regions close to the boundaries in the domain 
  # are still unoccupied. Hence, adding perturbation to all data points would 
  # corrupt the nature of the dataset.   
  for (j in 1:length(den.mat.d1.temp.unpert)) {
    # if(den.mat.d1.temp.unpert[j] != 0) {
     den.mat.d1.temp.pert[j] <- rnorm(1, den.mat.d1.temp.unpert[j], 
                                      sd = sd.r6.d1.min)
     while (den.mat.d1.temp.pert[j] < 0 || den.mat.d1.temp.pert[j] > 1) {
       den.mat.d1.temp.pert[j] <- rnorm(1, den.mat.d1.temp.unpert[j], sd = sd.r6.d1.min)
     }
    # }
    # if (den.mat.d3.temp.unpert[j] != 0) {
      den.mat.d3.temp.pert[j] <- rnorm(1, den.mat.d3.temp.unpert[j], sd = sd.r6.d3.min)
      while (den.mat.d3.temp.pert[j] < 0 || den.mat.d3.temp.pert[j] > 1) {
        den.mat.d3.temp.pert[j] <- rnorm(1, den.mat.d3.temp.unpert[j], sd = sd.r6.d3.min)
      }
    # }
  }
  
  # Sum of squared differences between the perturbed densities and the
  # observed densities
  sse.d1.temp <- sum((ref.den$t1.ref.den - den.mat.d1.temp.pert)^2)
  sse.d3.temp <- sum((ref.den$t3.ref.den - den.mat.d3.temp.pert)^2)
  diff.temp.pert <- sum(sse.d1.temp, sse.d3.temp) 
  
  # Store each set of perturbed data into the list
  paras.r6.output.pert[[i]] <- list(den.mat.d1.pert = den.mat.d1.temp.pert, 
                                    den.mat.d3.pert = den.mat.d3.temp.pert,
                                    diff = diff.temp.pert)
  # Name the elements in the list with correct indices
  names(paras.r6.output.pert)[i] <- paste0("ls_r6_",i,"_pert")
  print(i)
}

# Sort the sum of squared differences
ls.r6.pert <- matrix(0, nrow = length(ls.r6[,2])
                     ,ncol = 2)
for (i in 1:length(ls.r6[,2])) {
  ls.r6.pert[i,] <- c(i, as.double(paras.r6.output.pert[[i]]$diff))
}
ls.r6.pert.sort <- order(ls.r6.pert[,2])

# Locate the corresponding 200 perturbed model outputs that have the 
# minimum discrepancy with the reference data.
ls.r6.pert.sort.min <- ls.r6.pert.sort[1:200]
paras.r6.output.pert.full <- paras.r6.output.pert
paras.r6.output.pert <- paras.r6.output.pert[ls.r6.pert.sort.min]

# For each perturbed dataset, record its discrepancy with other unperturbed 
# model output, corresponding weights and resampling probabilities. 
ests <- foreach (i = 1:length(paras.r6.output.pert), .combine = rbind) %dopar% {
  # Create an empty matrix that stores the least square differences between 
  # the perturbed output and each of the unperturbed output. 
  diff.mat <- matrix(0, nrow = length(paras.r6[,2]), ncol = 4)
  paras.r6.output.pert.temp <- paras.r6.output.pert[[i]]
  for (j in 1:length(paras.r6[,2])) {
    # Read in the unperturbed model output.
    res.temp <- paras.r6.output[[j]]
    # Calculate the sum of squared differences between the unperturbed model output and the
    # selected set of perturbed data. 
    sse.d1.temp <- sum((res.temp$den.mat.d1 - 
                          paras.r6.output.pert.temp$den.mat.d1.pert)^2)
    sse.d3.temp <- sum((res.temp$den.mat.d3 - 
                          paras.r6.output.pert.temp$den.mat.d3.pert)^2)
    diff.temp <- sum(sse.d1.temp, sse.d3.temp)
    # Append index, discrepancy into each row of the matrix
    diff.mat[j, 1:4] <- c(j, diff.temp, sse.d1.temp, sse.d3.temp)
    # print(j)
  }
  # Remove the parameter set that corresponds to the current perturbed model output
  diff.mat <- diff.mat[-ls.r6.pert.sort.min[i],]
  
  # Obtain the info.list based on the diff.mat
  info.list <- calculate.bw(ss.mat = diff.mat, min.sd.d1 = sd.r6.d1.min,
                                min.sd.d3 = sd.r6.d3.min, 
                                ess.target = 2500, step.size = 0.0001)
  # Save the results for each "pseudo-reference dataset" into .rds files in 
  # the correct directory. 
  readr::write_rds(info.list, path = paste0("./", save.sims.dir, 
                                                 "/info_list_r6_", 
                                                 ls.r6.pert.sort.min[i], 
                                                 ".rds"))
  c(i, info.list$bw.obj, info.list$ess.obj)
  # Progress tracking, optional
  print(i)
}
stopCluster(cl)

# Read in the results obtained in the previous step, store them into one 
# single list. 
ls.r6.crage <- vector(mode = "list")
for (i in 1:length(ls.r6.pert.sort.min)) {
  info.list.temp <- read_rds(paste0("./", save.sims.dir, 
                                    "/info_list_r6_", 
                                    ls.r6.pert.sort.min[i], ".rds"))
  ls.r6.crage[[i]] <- info.list.temp
  names(ls.r6.crage)[i] <- paste0("info_list_r6_",
                                  ls.r6.pert.sort.min[i])
}


# Save the results into a .rds file
# write_rds(ls.r6.crage, "Coverage_test_r6_d1_d3_pert_sep_non_zero.rds")
write_rds(ls.r6.crage, "Coverage_test_r6_d1_d3_pert_sep.rds")

# Calculation of coverage probabilities and check if they show deviations 
# to uniform distributions. 

# Add indices to the parameter vectors.
paras.r6 <- cbind(seq(1,10000,by = 1), paras.r6)

# An empty matrix used to store all the coverage probabilities
# regarding each perturbed dataset. 
cov.mat <- matrix(0, nrow = 200, ncol = (length(paras.r6[1,]) - 1))

# Coverage probability calculations
for (i in 1:200) {
  # Read in the information matrix regarding each perturbed dataset.
  info.mat.temp <- ls.r6.crage[[i]]$info.mat
  # Remove the parameter set which generated the perturbed dataset, 
  # in order to keep the indices consistent
  paras.r6.corres <- paras.r6[-ls.r6.pert.sort.min[i],]
  # For each parameter, calculate its coverage probability
  prob.vec <- vector()
  for (j in 1:(length(paras.r6[1,]) - 1)) {
    cov.ind.temp <- which(paras.r6.corres[, (j + 1)] <= 
                            paras.r6[ls.r6.pert.sort.min[i], (j + 1)])
    prob.temp <- sum(info.mat.temp[cov.ind.temp, length(info.mat.temp[1,])])/
      sum(info.mat.temp[,length(info.mat.temp[1,])])
    prob.vec <- c(prob.vec, prob.temp)
  }
  # Store the coverage probability corresponding to each perturbed dataset into 
  # each row of the coverage probability matrix. 
  cov.mat[i, ] <- prob.vec
}

# Store the coverage probabilities
# write.table(cov.mat, 
# "Coverage probabilities glioma r6 d1 d3 perturbed separately non zero only.txt")
write.table(cov.mat, 
            "Coverage probabilities glioma r6 d1 d3 perturbed separately.txt")

cov.mat <- read.table("Coverage probabilities glioma r6 d1 d3 perturbed separately non zero only.txt", sep = "",
header = TRUE)

# Histograms plots & uniformity check.

# Reference uniform distribution. 
set.seed(874512)
unif.sample <- runif(200, 0, 1)

hist(cov.mat[,1], main = TeX("Coverage check of $d_{n}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,1], unif.sample)
# Weak evidence against H0 (non-zero)
# Moderate evidence against H0 (all)

hist(cov.mat[,2], main = TeX("Coverage check of $r_{n}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,2], unif.sample) 
# No evidence against H0 (non-zero). 
# Moderate evidence against H0. (all)

hist(cov.mat[,3], main = TeX("Coverage check of $R_{init.}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,3], unif.sample) 
# Weak evidence against H0 (non-zero). 
# Strong evidence against H0. (all)

hist(cov.mat[,4], main = TeX("Coverage check of $P_{ext.}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,4], unif.sample) 
# No evidence against H0 (non-zero). 
# No evidence against H0. (all)

hist(cov.mat[,5], main = TeX("Coverage check of $P_{mit.}$'s final sample",
                             bold = TRUE),
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,5], unif.sample) 
# Strong evidence against H0 (non-zero). 
# Strong evidence against H0. (all).

# Non zero data perturbed only:
# 2/5 parameters passed the coverage test! 
# 1/5 parameters shown strong evidence against H0.
# 2/5 parameter shown weak evidence against H0.

# All data perturbed: 
# 1/5 parameters passed the coverage test!
# 2/5 parameters shown strong evidence against H0.
# 2/5 parameters shown moderate evidence against H0. 
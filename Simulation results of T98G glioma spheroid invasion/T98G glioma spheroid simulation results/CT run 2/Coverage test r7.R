# Coverage test r7.R
# This .R file conducts the coverage test (an ABC posterior diagnostic that 
# checks if the final parameter samples can be regarded as appropriate samples
# from the perspective of Bayesian inference) on the final parameter samples 
# obtained in run 2. 

# Set the workspace, then load the necessary packages. 
rm(list = ls())
library(readr)
library(doParallel)
library(latex2exp)

# Set the directory to store the results. 
# save.sims.dir <- "Cov_results_r7"
save.sims.dir <- "Cov_results_r7_pert_all"
save.sims <- TRUE
if (save.sims) {
  if (!dir.exists(save.sims.dir)) dir.create(save.sims.dir)
}

# Set the cluster. 
n.thread <- detectCores()/2
cl <- makeCluster(n.thread)
registerDoParallel(cl)

# Calculate the correct bw that would yield the desirable ESS.
calculate.bw <- function(ss.mat, lb.bw, ub.bw, ess.target, step.size) {
  # Purpose: calculate the correct bandwidth factor which can yield the
  # desirable effective sample size (ESS).
  
  # Arguments: 
  # ss.mat: matrix of least square differences. 
  # lb.bw: lower bound of the decimal number sequence to be searched for
  #        the correct bandwidth factor. 
  # ub.bw: upper bound of the decimal number sequence to be searched for 
  #        the correct bandwidth factor. 
  # ess.target: desirable ESS. 
  # step.size: step of increment for the search of correct bandwidth factor
  #            in the sequence. 
  
  # Take away the singular outputs in the least square matrix. 
  ss.mat <- as.matrix(ss.mat)
  ind.nan <- which(is.na(ss.mat[,2]))
  if (length(ind.nan) == 0) {
    ss.mat.valid <- ss.mat
  } else {
    ss.mat.valid <- ss.mat[-ind.nan, ]
  }
  
  # Construct the sequence of bandwidth factor to be searched. 
  power <- seq(lb.bw, ub.bw, by = step.size)
  
  # An empty vector to store the ESS for each bandwidth factor. 
  ess.vec <- rep(0, length = length(power))
  
  # For every bandwidth factor in the sequence:
  for (i in 1:length(power)) {
    # Calculate the weights. 
    wt.temp <- ss.mat.valid[,2]^(-power[i])
    # Rescale the weights to resampling probabilities. 
    resamp.prob.temp <- rep(0, length = length(ss.mat.valid[,1]))
    
    for (j in 1:length(wt.temp)) {
      if (wt.temp[j] == min(wt.temp)) {
        resamp.prob.temp[j] <- 0
      } else if (wt.temp[j] == max(wt.temp)) {
        resamp.prob.temp[j] <- 1
      } else {
        resamp.prob.temp[j] <- (wt.temp[j]-min(wt.temp))/(max(wt.temp)-
                                                            min(wt.temp))
      }
    }
    
    # ESS calculation based on the rescaled weights obtained using the current
    # bandwidth factor. 
    ess.temp <- ((sum(resamp.prob.temp))^2)/sum(resamp.prob.temp^2)
    ess.vec[i] <- ess.temp
    # Optional, can be used to track the progress of the search. 
    # print(i)
  }
  
  # Matrix of bandwidth factors and their corresponding ESSs'.
  ess.mat <- cbind(power, ess.vec)
  # Take away the singular result. 
  ind.nan.ess <- which(is.na(ess.mat[,2]))
  if (length(ind.nan.ess) == 0) {
    ess.mat.valid <- ess.mat
  } else {
    ess.mat.valid <- ess.mat[-ind.nan.ess,]
  }
  
  # Locate the bandwidth factor which gives an ESS closest to the 
  # desirable one.
  ess.obj.ind <- which(abs(ess.mat.valid[,2] - ess.target) == 
                         min(abs(ess.mat.valid[,2] - ess.target)))
  ess.obj <- ess.mat.valid[ess.obj.ind,2]
  bw.obj <- ess.mat.valid[ess.obj.ind,1]
  
  # Calculate the corresponding weights and resampling probabilities using 
  # the bandwidth factor located previously. 
  wt.obj <- ss.mat.valid[,2]^(-bw.obj)
  resamp.prob.obj <- rep(0, length = length(ss.mat.valid[,1]))
  
  for (k in 1:length(wt.obj)) {
    if (wt.obj[k] == min(wt.obj)) {
      resamp.prob.obj[k] <- 0
    } else if (wt.obj[k] == max(wt.obj)) {
      resamp.prob.obj[k] <- 1
    } else {
      resamp.prob.obj[k] <- (wt.obj[k]-min(wt.obj))/(max(wt.obj)-min(wt.obj))
    }
  }
  
  # Output. 
  info.list <- list(info.mat = cbind(ss.mat.valid, wt.obj, resamp.prob.obj),
                    ess.mat = ess.mat,
                    ess.obj = ess.obj,
                    bw.obj = bw.obj)
  return(info.list)
}

# Read in the final parameter values. 
paras.r7 <- read.table("Round 7 parameters.txt", sep = "", header = TRUE)

# Set the directory to the correct one. 
sims.dir <- "LS_results_r7"

# Check if there are still NaNs in the output from posterior sample
ls.r7 <- read.table("Round 7 Least Square.txt", sep = "", header = TRUE)
length(which(is.na(ls.r7[,2]))) 
# 0, all output from the posterior sample are valid. 

# Observed data (glioma cell density)
ref.den <- read_rds("Glioma simulation full reference densities.rds")

# Read in the unperturbed model output and store them in a list. 
paras.r7.output <- vector(mode = "list")
for (i in 1:length(ls.r7[,2])) {
  res.temp <- read_rds(paste0("./", sims.dir, 
                              "/Round_7_paras",i,"_res.rds"))
  ls.temp.d1 <- res.temp$den.mat.d1
  ls.temp.d3 <- res.temp$den.mat.d3
  ls.temp <- list(den.mat.d1 = ls.temp.d1,
                  den.mat.d3 = ls.temp.d3)
  paras.r7.output[[i]] <- ls.temp
  print(i)
  names(paras.r7.output)[i] <- paste0("ls_r7_",i)
}

# Write the unperturbed model output into a .rds file. 
write_rds(paras.r7.output, "Full output paras r7 run 2.rds")


# Read in the output in round 7 that gave the minimum result of least square
# difference, and calculate the standard deviation (SD) based on it.
ls.r7.min <- read_rds("Round_7_paras425_res.rds")
sd.r7.d1.min <- sd(ls.r7.min$den.mat.d1 - ref.den$t1.ref.den)
sd.r7.d3.min <- sd(ls.r7.min$den.mat.d3 - ref.den$t3.ref.den)

# Add perturbations (gaussian) to the final unperturbed output based on the
# SD calculated in the previous step. 
paras.r7.output.pert <- vector(mode = "list")

set.seed(874512)
RNGkind(sample.kind = "Rejection")
for (i in 1:length(paras.r7.output)) {
  ls.temp <- paras.r7.output[[i]]
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
                                      sd = sd.r7.d1.min)
     while (den.mat.d1.temp.pert[j] < 0 || den.mat.d1.temp.pert[j] > 1) {
       den.mat.d1.temp.pert[j] <- rnorm(1, den.mat.d1.temp.unpert[j], 
                                        sd = sd.r7.d1.min)
     }
    # }
    # if (den.mat.d3.temp.unpert[j] != 0) {
      den.mat.d3.temp.pert[j] <- rnorm(1, den.mat.d3.temp.unpert[j], 
                                       sd = sd.r7.d3.min)
      while (den.mat.d3.temp.pert[j] < 0 || den.mat.d3.temp.pert[j] > 1) {
        den.mat.d3.temp.pert[j] <- rnorm(1, den.mat.d3.temp.unpert[j], 
                                         sd = sd.r7.d3.min)
      }
    # }
  }
  
  # Sum of squared differences between the perturbed densities and the
  # observed densities
  sse.d1.temp <- sum((ref.den$t1.ref.den - den.mat.d1.temp.pert)^2)
  sse.d3.temp <- sum((ref.den$t3.ref.den - den.mat.d3.temp.pert)^2)
  diff.temp.pert <- sum(sse.d1.temp, sse.d3.temp) 
  
  # Store each set of perturbed data into the list
  paras.r7.output.pert[[i]] <- list(den.mat.d1.pert = den.mat.d1.temp.pert, 
                                    den.mat.d3.pert = den.mat.d3.temp.pert,
                                    diff = diff.temp.pert)
  # Name the elements in the list with correct indices
  names(paras.r7.output.pert)[i] <- paste0("ls_r7_",i,"_pert")
  print(i)
}

# Sort the sum of squared differences
ls.r7.pert <- matrix(0, nrow = length(ls.r7[,2])
                     ,ncol = 2)
for (i in 1:length(ls.r7[,2])) {
  ls.r7.pert[i,] <- c(i, as.double(paras.r7.output.pert[[i]]$diff))
}
ls.r7.pert.sort <- order(ls.r7.pert[,2])

# Locate the corresponding 200 perturbed model outputs that have the 
# minimum discrepancy with the reference data.
ls.r7.pert.sort.min <- ls.r7.pert.sort[1:200]
paras.r7.output.pert.full <- paras.r7.output.pert
paras.r7.output.pert <- paras.r7.output.pert[ls.r7.pert.sort.min]

# Discrepancies between unperturbed model output with each of these 200 
# perturbed model output which were treated as new reference data. 
ls.r7.diff.mat <- vector(mode = "list")
# For each perturbed dataset, record its discrepancy with other unperturbed 
# model output.
for (i in 1:length(paras.r7.output.pert)) {
  # Create an empty matrix that stores the least square differences between 
  # the perturbed output and each of the unperturbed output. 
  diff.mat <- matrix(0, nrow = length(paras.r7[,2]), ncol = 2)
  paras.r7.output.pert.temp <- paras.r7.output.pert[[i]]
  for (j in 1:length(paras.r7[,2])) {
    # Read in the unperturbed model output.
    res.temp <- paras.r7.output[[j]]
    # Calculate the sum of squared differences between the unperturbed model 
    # output and the selected set of perturbed data. 
    sse.d1.temp <- sum((res.temp$den.mat.d1 - 
                          paras.r7.output.pert.temp$den.mat.d1.pert)^2)
    sse.d3.temp <- sum((res.temp$den.mat.d3 - 
                          paras.r7.output.pert.temp$den.mat.d3.pert)^2)
    diff.temp <- sum(sse.d1.temp, sse.d3.temp)
    # Append index, discrepancy into each row of the matrix
    diff.mat[j, 1:2] <- c(j, diff.temp)
    # Optional, can be used to track the progress.
    # print(j)
  }
  # Remove the parameter set that generated the current perturbed model output
  diff.mat <- diff.mat[-ls.r7.pert.sort.min[i],]
  
  # Store the least square difference matrix in the list and name it with 
  # the corresponding index.   
  ls.r7.diff.mat[[i]] <- diff.mat
  names(ls.r7.diff.mat)[i] <- paste0("diff_mat_r7_",ls.r7.pert.sort.min[i])
  # Progress tracking
  print(i)
}

# Carry out another round of ABC but using these perturbed output as reference
# data. Calculate the resampling probabilities of each parameter vector based
# on the least square difference between their simulation output and the "new
# reference data". 

ests <- foreach (i = 1:length(ls.r7.diff.mat), .combine = rbind) %dopar% {
  diff.mat.temp <- ls.r7.diff.mat[[i]]
  # Initial interval for the bandwidth factors. 
  lb.bw.temp <- 11.9 #4.1
  ub.bw.temp <- 12.0 #4.2
  info.list.temp <- calculate.bw(ss.mat = diff.mat.temp, lb.bw = lb.bw.temp,
                                 ub.bw = ub.bw.temp, ess.target = 2250, 
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
                                     ub.bw = ub.bw.temp, ess.target = 2250,
                                     step.size = 0.01)
    } else if (info.list.temp$bw.obj == ub.bw.temp) {
      lb.bw.temp <- ub.bw.temp - 0.01
      ub.bw.temp <- ub.bw.temp + 0.1
      info.list.temp <- calculate.bw(ss.mat = diff.mat.temp, 
                                     lb.bw = lb.bw.temp,
                                     ub.bw = ub.bw.temp, ess.target = 2250,
                                     step.size = 0.01)
    }
  }
  # Save the results for each "pseudo-reference dataset" into .rds files in 
  # the correct directory. 
  readr::write_rds(info.list.temp, path = paste0("./", save.sims.dir, 
                                                 "/info_list_r7_", 
                                                 ls.r7.pert.sort.min[i], 
                                                 ".rds"))
  c(i, info.list.temp$bw.obj, info.list.temp$ess.obj)
}
stopCluster(cl)

# Read in the results obtained in the previous step, store them into one 
# single list. 
ls.r7.crage <- vector(mode = "list")
for (i in 1:length(ls.r7.diff.mat)) {
  info.list.temp <- read_rds(paste0("./", save.sims.dir, 
                                    "/info_list_r7_", 
                                    ls.r7.pert.sort.min[i], ".rds"))
  ls.r7.crage[[i]] <- info.list.temp
  names(ls.r7.crage)[i] <- paste0("info_list_r7_",
                                  ls.r7.pert.sort.min[i])
}


# Save the results into a .rds file
# write_rds(ls.r7.crage, "Coverage_test_d1_d3_pert_non_zero.rds")
write_rds(ls.r7.crage, "Coverage_test_d1_d3_pert_all.rds")
# ls.r7.crage <- 
# read_rds("Coverage_test_r7_d1_d3_pert_sep_non_zero.rds")

# Calculation of coverage probabilities and check if they show deviations 
# to uniform distributions. 

# Add indices to the parameter vectors.
paras.r7 <- cbind(seq(1,10000,by = 1), paras.r7)

# An empty matrix used to store all the coverage probabilities
# regarding each perturbed dataset. 
cov.mat <- matrix(0, nrow = 200, ncol = (length(paras.r7[1,]) - 1))

# Coverage probabilities calculations
for (i in 1:200) {
  # Read in the information matrix regarding each perturbed dataset.
  info.mat.temp <- ls.r7.crage[[i]]$info.mat
  # Remove the parameter set which generated the perturbed dataset.
  # in order to keep the indices consistent.
  paras.r7.corres <- paras.r7[-ls.r7.pert.sort.min[i],]
  # For each parameter, calculate its coverage probability.
  prob.vec <- vector()
  for (j in 1:(length(paras.r7[1,]) - 1)) {
    cov.ind.temp <- which(paras.r7.corres[, (j + 1)] <= 
                            paras.r7[ls.r7.pert.sort.min[i], (j + 1)])
    prob.temp <- sum(info.mat.temp[cov.ind.temp, length(info.mat.temp[1,])])/
      sum(info.mat.temp[,length(info.mat.temp[1,])])
    prob.vec <- c(prob.vec, prob.temp)
  }
  # Store the coverage probability corresponding to each perturbed dataset into 
  # each row of the coverage probability matrix. 
  cov.mat[i, ] <- prob.vec
}

# Store the coverage probabilities
# write.table(cov.mat, "Coverage probabilities glioma r7 d1 d3 perturbed non zero only.txt")
write.table(cov.mat, "Coverage probabilities glioma r7 d1 d3 perturbed separately.txt")

cov.mat <- read.table("Coverage probabilities glioma r7 d1 d3 perturbed non zero only.txt", sep = "", header = TRUE)

# Histograms plots & uniformity check.

# Reference uniform distribution. 
set.seed(874512)
unif.sample <- runif(200, 0, 1)

hist(cov.mat[,1], main = TeX("Coverage check of $d_{n}$'s final sample",
                             bold = TRUE), 
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,1], unif.sample)
# No evidence against H0 (non-zero)
# Weak evidence against H0 (all)

hist(cov.mat[,2], main = TeX("Coverage check of $r_{n}$'s final sample",
                             bold = TRUE), 
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,2], unif.sample) 
# No evidence against H0 (non-zero). 
# No evidence against H0. (all)

hist(cov.mat[,3], main = TeX("Coverage check of $R_{init.}$'s final sample",
                             bold = TRUE), 
     xlim = c(0,1), xlab = "Probabilities", freq = FALSE)
ks.test(cov.mat[,3], unif.sample) 
# No evidence against H0 (non-zero). 
# Moderate evidence against H0. (all)

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
# No evidence against H0 (non-zero). 
# Strong evidence against H0. (all).

# Non zero data perturbed only:
# 5/5 parameters passed the coverage test! 


# All data perturbed: 
# 2/5 parameters passed the coverage test!
# 1/5 parameters shown strong evidence against H0.
# 1/5 parameters shown weak evidence against H0.
# 1/5 parameters shown moderate evidence against H0.
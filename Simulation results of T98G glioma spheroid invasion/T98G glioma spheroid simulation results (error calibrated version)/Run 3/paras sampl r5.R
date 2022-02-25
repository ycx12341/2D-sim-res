# paras sampl r5.R
# Author: Yunchen Xiao
# This .R file reads in the parameters used in round 5 of the simulations T98G
# glioma invasion patterns and returns the parameters to be evaluated in 
# round 5 (error-calibrated version). 

# Set the workspace, then load the necessary packages. 
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file which contains the source functions.
source("PDE 2D ABC functions adjusted.r")

# Set the directory to store the simulation results for each parameter vector
# generated. 
save.sims.dir <- "LS_results_r5"
save.sims <- TRUE
if(save.sims) {
  if(!dir.exists(save.sims.dir)) dir.create(save.sims.dir)
}

# Set the cluster.
n.thread <- detectCores()/2
n.sims <- 10000
cl <- makeCluster(n.thread)
registerDoParallel(cl)

# Read in the parameters to be evaluated in the current round. 
paras.r5 <- as.matrix(read.table("Round 5 parameters log transform.txt", 
                                 sep = "", header = TRUE))

# Simulation running in parallel. 
registerDoRNG(874514, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r5[i, 1:2], 
                                init.cells.rad = paras.r5[i, 3],
                                prob.death = paras.r5[i, 4], 
                                prob.prof = paras.r5[i, 5])
  # Write the results into .rds files, stored in the directory specified 
  # earlier. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_5_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 5835.081 sec elapsed. 

# Read in the least square differences of the parameter vectors. 
ls.r5 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_5_paras",i,"_res.rds"))
  ls.temp <- c(res.temp$diff, res.temp$sse.t1, res.temp$sse.t3)
  ls.r5 <- rbind(ls.r5, ls.temp)
  print(i)
}
ls.r5 <- unname(ls.r5)

# Combine the least square differences with their corresponding indices by 
# column. 
ls.r5 <- cbind(seq(1, n.sims, by = 1), ls.r5)

# Check if there are singular values in the output and calculate the average 
# least square difference. 
ind.nan <- which(is.na(ls.r5[, 2]))
mean(ls.r5[,2]) # 1.603679

# Locate the minimum one. 
# minimum : 4758th 0.9519926

# Write the least square matrix into a .txt file. 
write.table(ls.r5, "Round 5 Least Square.txt")

# Minimum result
ls.r5.min <- read_rds("Round_5_paras4758_res.rds")

# Calculate the SDs of the discrepancy between the best-fitted output
# and the reference data. 
den.mat.d1.r5.min <- ls.r5.min$den.mat.d1
den.mat.d3.r5.min <- ls.r5.min$den.mat.d3
sd.den.mat.d1.r5.min <- sd(den.mat.d1.r5.min - t1.ref.den)
sd.den.mat.d3.r5.min <- sd(den.mat.d3.r5.min - t3.ref.den)

# Compute the bandwidth factor that is used to calculate the resampling weights,
# based on the desirable effective sample size (ESS).
info.list.r5 <- calculate.bw(ss.mat = ls.r5, min.sd.d1 = sd.den.mat.d1.r5.min,
                             min.sd.d3 = sd.den.mat.d3.r5.min, 
                             ess.target = 1640.25, step.size = 0.0001)
write_rds(info.list.r5, "Round 5 information list log transform.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r6 <- abc_bcd(info.mat = info.list.r5$info.mat, paras.r5)
write.table(paras.r6, "Round 6 parameters log transform.txt")

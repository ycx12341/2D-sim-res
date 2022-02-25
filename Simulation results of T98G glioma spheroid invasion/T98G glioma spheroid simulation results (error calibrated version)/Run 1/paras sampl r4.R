# paras sampl r4.R
# Author: Yunchen Xiao
# This .R file reads in the parameters used in round 4 of the simulations T98G
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
save.sims.dir <- "LS_results_r4"
save.sims <- TRUE
if(save.sims) {
  if(!dir.exists(save.sims.dir)) dir.create(save.sims.dir)
}

# Set the cluster.
n.thread <- detectCores() - 1
n.sims <- 10000
cl <- makeCluster(n.thread)
registerDoParallel(cl)

# Read in the parameters to be evaluated in the current round. 
paras.r4 <- as.matrix(read.table("Round 4 parameters log transform.txt", 
                                 sep = "", header = TRUE))

# Simulation running in parallel. 
registerDoRNG(874513, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r4[i, 1:2], 
                                init.cells.rad = paras.r4[i, 3],
                                prob.death = paras.r4[i, 4], 
                                prob.prof = paras.r4[i, 5])
  # Write the results into .rds files, stored in the directory specified 
  # earlier. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_4_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 11926.1 sec elapsed. 

# Read in the least square differences of the parameter vectors. 
ls.r4 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_4_paras",i,"_res.rds"))
  ls.temp <- c(res.temp$diff, res.temp$sse.t1, res.temp$sse.t3)
  ls.r4 <- rbind(ls.r4, ls.temp)
  print(i)
}
ls.r4 <- unname(ls.r4)

# Combine the least square differences with their corresponding indices by 
# column. 
ls.r4 <- cbind(seq(1, n.sims, by = 1), ls.r4)

# Take away the singular results (temporarily) and calculate the average 
# least square difference.
ind.nan <- which(is.na(ls.r4[, 2]))
ls.r4.no.nan <- ls.r4[-ind.nan, ]
mean(ls.r4[,2]) # 1.954337

# Locate the minimum one. 
# minimum : 7421th 0.9706663

# Write the least square matrix into a .txt file. 
write.table(ls.r4, "Round 4 Least Square.txt")

# Minimum result
ls.r4.min <- read_rds("Round_4_paras7421_res.rds")

# Calculate the SDs of the discrepancy between the best-fitted output
# and the reference data. 
den.mat.d1.r4.min <- ls.r4.min$den.mat.d1
den.mat.d3.r4.min <- ls.r4.min$den.mat.d3
sd.den.mat.d1.r4.min <- sd(den.mat.d1.r4.min - t1.ref.den)
sd.den.mat.d3.r4.min <- sd(den.mat.d3.r4.min - t3.ref.den)

# Compute the bandwidth factor that is used to calculate the resampling weights,
# based on the desirable effective sample size (ESS).
info.list.r4 <- calculate.bw(ss.mat = ls.r4, min.sd.d1 = sd.den.mat.d1.r4.min,
                             min.sd.d3 = sd.den.mat.d3.r4.min, 
                             ess.target = 1822.5, step.size = 0.0001)
write_rds(info.list.r4, "Round 4 information list log transform.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r5 <- abc_bcd(info.mat = info.list.r4$info.mat, paras.r4)
write.table(paras.r5, "Round 5 parameters log transform.txt")

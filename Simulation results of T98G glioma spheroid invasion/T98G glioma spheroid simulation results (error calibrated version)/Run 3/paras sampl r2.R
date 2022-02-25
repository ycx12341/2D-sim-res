# paras sampl r2.R
# Author: Yunchen Xiao
# This .R file reads in the parameters used in round 2 of the simulations T98G
# glioma invasion patterns and returns the parameters to be evaluated in 
# round 3 (error-calibrated version). 

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
save.sims.dir <- "LS_results_r2"
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
paras.r2 <- as.matrix(read.table("Round 2 parameters log transform.txt", 
                                 sep = "", header = TRUE))

# Simulation running in parallel. 
registerDoRNG(874514, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r2[i, 1:2], 
                                init.cells.rad = paras.r2[i, 3],
                                prob.death = paras.r2[i, 4], 
                                prob.prof = paras.r2[i, 5])
  # Write the results into .rds files, stored in the directory specified 
  # earlier. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_2_paras", i, "_res.rds"))
  
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 10231.03 sec elapsed. 

# Read in the least square differences of the parameter vectors. 
ls.r2 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_2_paras",i,"_res.rds"))
  ls.temp <- c(res.temp$diff, res.temp$sse.t1, res.temp$sse.t3)
  ls.r2 <- rbind(ls.r2, ls.temp)
  print(i)
}
ls.r2 <- unname(ls.r2)

# Combine the least square differences with their corresponding indices by 
# column. 
ls.r2 <- cbind(seq(1, n.sims, by = 1), ls.r2)

# Take away the singular results (temporarily) and calculate the average 
# least square difference. 
ind.nan <- which(is.na(ls.r2[, 2]))
ls.r2.no.nan <- ls.r2[-ind.nan, ]
mean(ls.r2.no.nan[,2]) # 4.709282

# Locate the minimum one. 
# minimum : 6804th 1.127705

# Write the least square matrix into a .txt file. 
write.table(ls.r2, "Round 2 Least Square.txt")

# Minimum result
ls.r2.min <- read_rds("Round_2_paras6804_res.rds")

# Calculate the SDs of the discrepancy between the best-fitted output
# and the reference data. 
den.mat.d1.r2.min <- ls.r2.min$den.mat.d1
den.mat.d3.r2.min <- ls.r2.min$den.mat.d3
sd.den.mat.d1.r2.min <- sd(den.mat.d1.r2.min - t1.ref.den)
sd.den.mat.d3.r2.min <- sd(den.mat.d3.r2.min - t3.ref.den)

# Compute the bandwidth factor that is used to calculate the resampling weights,
# based on the desirable effective sample size (ESS).
info.list.r2 <- calculate.bw(ss.mat = ls.r2, 
                             min.sd.d1 = sd.den.mat.d1.r2.min,
                             min.sd.d3 = sd.den.mat.d3.r2.min, 
                             ess.target = 2250, step.size = 0.0001)
write_rds(info.list.r2, "Round 2 information list log transform.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r3 <- abc_bcd(info.mat = info.list.r2$info.mat, paras.r2)
write.table(paras.r3, "Round 3 parameters log transform.txt")

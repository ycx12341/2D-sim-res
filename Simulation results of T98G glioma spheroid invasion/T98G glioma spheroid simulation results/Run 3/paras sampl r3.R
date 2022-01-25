# paras sampl r3.R
# Author: Yunchen Xiao
# This .R file reads in the parameters used in round 3 of the simulations T98G
# glioma invasion patterns and returns the parameters to be evaluated in 
# round 4.

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
save.sims.dir <- "LS_results_r3"
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
paras.r3 <- as.matrix(read.table("Round 3 parameters.txt", sep = "", 
                                 header = TRUE))

# Simulation running in parallel. 
registerDoRNG(874514, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r3[i, 1:2], 
                                init.cells.rad = paras.r3[i, 3],
                                prob.death = paras.r3[i, 4], 
                                prob.prof = paras.r3[i, 5])
  # Write the results into .rds files, stored in the directory specified 
  # earlier. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_3_paras", i, "_res.rds"))
  
  c(i, sse.res.temp$diff)
}
toc()

stopCluster(cl)

# 11839.06 sec elapsed. 

# Read in the least square differences of the parameter vectors. 
ls.r3 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_3_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r3 <- rbind(ls.r3, ls.temp)
  print(i)
}
ls.r3 <- unname(ls.r3)

# Combine the least square differences with their corresponding indices by 
# column.
ls.r3 <- cbind(seq(1, n.sims, by = 1), ls.r3)

# Take away the singular results (temporarily) and calculate the average 
# least square difference.
ind.nan <- which(is.na(ls.r3[, 2]))
ls.r3.no.nan <- ls.r3[-ind.nan, ]
mean(ls.r3.no.nan[,2]) # 3.64384

# Locate the minimum one. 
# minimum : 6813th 0.9996551

# Write the least square matrix into a .txt file. 
write.table(ls.r3, "Round 3 Least Square.txt")

# Compute the bandwidth factor that is used to calculate the resampling weights,
# based on the desirable effective sample size (ESS).
info.list.r3 <- calculate.bw(ss.mat = ls.r3, lb.bw = 3.9, ub.bw = 4.0, 
                             ess.target = 2025, step.size = 0.01)
write_rds(info.list.r3, "Round 3 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r4 <- abc_bcd(info.mat = info.list.r3$info.mat, paras.r3)
write.table(paras.r4, "Round 4 parameters.txt")

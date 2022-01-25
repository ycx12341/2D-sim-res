# paras sampl r6.R
# Author: Yunchen Xiao
# This .R file reads in the parameters used in round 5 of the simulations T98G
# glioma invasion patterns and returns the parameters to be evaluated in 
# round 7.

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
save.sims.dir <- "LS_results_r6"
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
paras.r6 <- as.matrix(read.table("Round 6 parameters.txt", sep = "", 
                                 header = TRUE))

# Simulation running in parallel. 
registerDoRNG(874514, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r6[i, 1:2], 
                                init.cells.rad = paras.r6[i, 3],
                                prob.death = paras.r6[i, 4], 
                                prob.prof = paras.r6[i, 5])
  # Write the results into .rds files, stored in the directory specified 
  # earlier. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_6_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 4685.375 sec elapsed. 

# Read in the least square differences of the parameter vectors. 
ls.r6 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_6_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r6 <- rbind(ls.r6, ls.temp)
  print(i)
}
ls.r6 <- unname(ls.r6)

# Combine the least square differences with their corresponding indices by 
# column.
ls.r6 <- cbind(seq(1, n.sims, by = 1), ls.r6)

# Check if the results contain singular values and calculate the average least
# square difference. 
ind.nan <- which(is.na(ls.r6[, 2]))
mean(ls.r6[,2]) # 1.549375

# Locate the minimum one. 
# minimum : 6056th 0.8525401

# Write the least square matrix into a .txt file. 
write.table(ls.r6, "Round 6 Least Square.txt")

# Compute the bandwidth factor that is used to calculate the resampling weights,
# based on the desirable effective sample size (ESS).
info.list.r6 <- calculate.bw(ss.mat = ls.r6, lb.bw = 8.1, ub.bw = 8.4, 
                             ess.target = 2500, step.size = 0.01)
write_rds(info.list.r6, "Round 6 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r7 <- abc_bcd(info.mat = info.list.r6$info.mat, paras.r6)
write.table(paras.r7, "Round 7 parameters.txt")






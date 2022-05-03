# paras sampl r6.R
# Author: Yunchen Xiao
# This .R file reads in the parameter values to be evaluated for the 
# simulation of SCC pattern in the current round and returns the parameters 
# to be evaluated in round 7. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC.R")

# Set the directory to store the simulation results. 
save.sims.dir <- "LS_results_r6"
save.sims <- TRUE
if(save.sims) {
  if(!dir.exists(save.sims.dir)) dir.create(save.sims.dir)
}

# Set the cluster. 
n.thread <- detectCores() - 1
n.sims <- 10000
cl <- makeCluster(n.thread)
registerDoParallel(cl)

# Read in the parameter values to be evaluated in the current round. 
paras.r6 <- as.matrix(read.table("Round 6 parameters.txt", sep = "",
                                 header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874512, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r6[i, 1:6], 
                                init.cells.cols = paras.r6[i, 7],
                                prob.death = paras.r6[i, 8], 
                                prob.prof = paras.r6[i, 9])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_6_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 4922.24 sec elapsed. 

# Read in the least square differences.
ls.r6 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_6_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r6 <- rbind(ls.r6, ls.temp)
  print(i)
}
ls.r6 <- unname(ls.r6)

# Combine the least square differences with the indices by column. 
ls.r6 <- cbind(seq(1, n.sims, by = 1), ls.r6)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r6[, 2]))
ls.r6.no.nan <- ls.r6[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r6.no.nan[,2]) # 0.4858792
# minimum :5492th 0.09910419 

# Write the least square differences into a .txt file.  
write.table(ls.r6, "Round 6 Least Square.txt")
# ls.r6 <- as.matrix(read.table("Round 6 Least Square.txt", sep = "",
#                               header = TRUE))

# Calculate the bandwidth factor which gives the desirable ESS.
info.list.r6 <- calculate.bw(ss.mat = ls.r6, lb.bw = 3, ub.bw = 3.5, 
                             ess.target = 1640.25, step.size = 0.01)
write_rds(info.list.r6, "Round 6 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r7 <- abc_bcd(info.mat = info.list.r6$info.mat, paras.r6)
write.table(paras.r7, "Round 7 parameters.txt")

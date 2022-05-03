# paras sampl r3.R
# Author: Yunchen Xiao
# This .R file reads in the parameter values to be evaluated for the 
# simulation of SCC pattern in the current round and returns the parameters 
# to be evaluated in round 4. 

# Clear the workspace and load the necessary packages.
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC.R")

# Set the directory to store the simulation results. 
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

# Read in the parameter values to be evaluated in the current round. 
paras.r3 <- as.matrix(read.table("Round 3 parameters.txt", sep = "",
                                 header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874512, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r3[i, 1:6], 
                                prob.death = paras.r3[i, 7], 
                                prob.prof = paras.r3[i, 8])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_3_paras", i, "_res.rds"))
  
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 12907.859 sec elapsed. 

# Read in the least square differences.
ls.r3 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_3_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r3 <- rbind(ls.r3, ls.temp)
  print(i)
}
ls.r3 <- unname(ls.r3)

# Combine the least square differences with the indices by column. 
ls.r3 <- cbind(seq(1, n.sims, by = 1), ls.r3)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r3[, 2]))
ls.r3.no.nan <- ls.r3[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r3.no.nan[,2]) # 0.6278714
# minimum : 4347th 0.1928125

# Write the least square differences into a .txt file.  
write.table(ls.r3, "Round 3 Least Square.txt")
# ls.r3 <- as.matrix(read.table("Round 3 Least Square.txt", sep = "",
#                               header = TRUE))

# Calculate the bandwidth factor which gives the desirable ESS.
info.list.r3 <- calculate.bw(ss.mat = ls.r3, lb.bw = 5, ub.bw = 5.2, 
                             ess.target = 2025, step.size = 0.01)
write_rds(info.list.r3, "Round 3 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r4 <- abc_bcd(info.mat = info.list.r3$info.mat, paras.r3)
write.table(paras.r4, "Round 4 parameters.txt")

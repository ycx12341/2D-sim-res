# paras sampl r10.R
# Author: Yunchen Xiao
# This .R file reads in the parameter values to be evaluated for the 
# simulation of SCC pattern in the current round and returns the parameters 
# to be evaluated in round 11. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC.R")

# Set the directory to store the simulation results. 
save.sims.dir <- "LS_results_r10"
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
paras.r10 <- as.matrix(read.table("Round 10 parameters.txt", sep = "",
                                 header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874513, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r10[i, 1:6], 
                                init.cells.cols = paras.r10[i, 7],
                                prob.death = paras.r10[i, 8], 
                                prob.prof = paras.r10[i, 9])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_10_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 4947.91 sec elapsed. 

# Read in the least square differences.
ls.r10 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_10_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r10 <- rbind(ls.r10, ls.temp)
  print(i)
}
ls.r10 <- unname(ls.r10)

# Combine the least square differences with the indices by column. 
ls.r10 <- cbind(seq(1, n.sims, by = 1), ls.r10)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r10[, 2]))
ls.r10.no.nan <- ls.r10[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r10.no.nan[,2]) # 0.2211397
# minimum : 1209th 0.06635699 

# ls.r10 <- as.matrix(read.table("Round 10 Least Square.txt", sep = "",
#                               header = TRUE))

# Write the least square differences into a .txt file.  
write.table(ls.r10, "Round 10 Least Square.txt")

# Calculate the bandwidth factor which gives the desirable ESS.
info.list.r10 <- calculate.bw(ss.mat = ls.r10, lb.bw = 4, ub.bw = 4.5, 
                             ess.target = 1822.5, step.size = 0.01)
write_rds(info.list.r10, "Round 10 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r11 <- abc_bcd(info.mat = info.list.r10$info.mat, paras.r10)
write.table(paras.r11, "Round 11 parameters.txt")

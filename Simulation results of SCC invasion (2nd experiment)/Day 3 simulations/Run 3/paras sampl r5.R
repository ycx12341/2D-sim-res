# paras sampl r5.R
# Author: Yunchen Xiao
# This .R file reads in the parameter values to be evaluated for the 
# simulation of SCC pattern in the current round and returns the parameters 
# to be evaluated in round 6. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC.R")

# Set the directory to store the simulation results. 
save.sims.dir <- "LS_results_r5"
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
paras.r5 <- as.matrix(read.table("Round 5 parameters.txt", sep = "",
                                 header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874512, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r5[i, 1:6], 
                                init.cells.cols = paras.r5[i, 7],
                                prob.death = paras.r5[i, 8], 
                                prob.prof = paras.r5[i, 9])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_5_paras", i, "_res.rds"))
  
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 4828.87 sec elapsed. 

# Read in the least square differences.
ls.r5 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_5_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r5 <- rbind(ls.r5, ls.temp)
  print(i)
}
ls.r5 <- unname(ls.r5)

# Combine the least square differences with the indices by column. 
ls.r5 <- cbind(seq(1, n.sims, by = 1), ls.r5)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r5[, 2]))
ls.r5.no.nan <- ls.r5[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r5.no.nan[,2]) # 0.7785499
# minimum : 2702nd 0.1092530 

# Write the least square differences into a .txt file.  
write.table(ls.r5, "Round 5 Least Square.txt")
# ls.r5 <- as.matrix(read.table("Round 5 Least Square.txt", sep = "",
#                               header = TRUE))

# Calculate the bandwidth factor which gives the desirable ESS.
info.list.r5 <- calculate.bw(ss.mat = ls.r5, lb.bw = 2.5, ub.bw = 3, 
                             ess.target = 1822.5, step.size = 0.01)
write_rds(info.list.r5, "Round 5 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r6 <- abc_bcd(info.mat = info.list.r5$info.mat, paras.r5)
write.table(paras.r6, "Round 6 parameters.txt")

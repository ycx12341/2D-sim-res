# paras sampl r4.R
# Author: Yunchen Xiao
# This .R file reads in the parameter values to be evaluated for the 
# simulation of SCC pattern in the current round and returns the parameters 
# to be evaluated in round 5. 

# Clear the workspace and load the necessary packages.
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC.R")

# Set the directory to store the simulation results. 
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

# Read in the parameter values to be evaluated in the current round. 
paras.r4 <- as.matrix(read.table("Round 4 parameters.txt", sep = "",
                                 header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874514, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r4[i, 1:6],
                                prob.death = paras.r4[i, 7], 
                                prob.prof = paras.r4[i, 8])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_4_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 21213.86 sec elapsed. 

# Read in the least square differences.
ls.r4 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_4_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r4 <- rbind(ls.r4, ls.temp)
  print(i)
}
ls.r4 <- unname(ls.r4)

# Combine the least square differences with the indices by column. 
ls.r4 <- cbind(seq(1, n.sims, by = 1), ls.r4)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r4[, 2]))
ls.r4.no.nan <- ls.r4[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r4.no.nan[,2]) # 0.6062692
# minimum : 88th 0.2461443 

# Write the least square differences into a .txt file.  
write.table(ls.r4, "Round 4 Least Square.txt")
# ls.r4 <- as.matrix(read.table("Round 4 Least Square.txt", sep = "",
#                               header = TRUE))

# Calculate the bandwidth factor which gives the desirable ESS.
info.list.r4 <- calculate.bw(ss.mat = ls.r4, lb.bw = 6, ub.bw = 6.5, 
                             ess.target = 1822.5, step.size = 0.01)
write_rds(info.list.r4, "Round 4 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r5 <- abc_bcd(info.mat = info.list.r4$info.mat, paras.r4)
write.table(paras.r5, "Round 5 parameters.txt")

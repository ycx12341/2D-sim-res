# paras sampl r8.R
# Author: Yunchen Xiao
# This .R file reads in the parameter values to be evaluated for the 
# simulation of SCC pattern in the current round and returns the parameters 
# to be evaluated in round 9. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC.R")

# Set the directory to store the simulation results. 
save.sims.dir <- "LS_results_r8"
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
paras.r8 <- as.matrix(read.table("Round 8 parameters.txt", sep = "",
                                 header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874514, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r8[i, 1:6], 
                                init.cells.cols = paras.r8[i, 7],
                                prob.death = paras.r8[i, 8], 
                                prob.prof = paras.r8[i, 9])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_8_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 5563.65 sec elapsed. 

# Read in the least square differences.
ls.r8 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_8_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r8 <- rbind(ls.r8, ls.temp)
  print(i)
}
ls.r8 <- unname(ls.r8)

# Combine the least square differences with the indices by column. 
ls.r8 <- cbind(seq(1, n.sims, by = 1), ls.r8)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r8[, 2]))
ls.r8.no.nan <- ls.r8[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r8.no.nan[,2]) # 0.2504458
# minimum : 4875th 0.07675459

# Write the least square differences into a .txt file.  
write.table(ls.r8, "Round 8 Least Square.txt")

# ls.r8 <- as.matrix(read.table("Round 8 Least Square.txt", sep = "",
#                               header = TRUE))

# Calculate the bandwidth factor which gives the desirable ESS.
info.list.r8 <- calculate.bw(ss.mat = ls.r8, lb.bw = 4, ub.bw = 4.5, 
                             ess.target = 2250, step.size = 0.01)
write_rds(info.list.r8, "Round 8 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r9 <- abc_bcd(info.mat = info.list.r8$info.mat, paras.r8)
write.table(paras.r9, "Round 9 parameters.txt")

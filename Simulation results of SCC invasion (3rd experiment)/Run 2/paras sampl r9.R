# paras sampl r9.R
# Author: Yunchen Xiao
# This .R file reads in the post-round 8 parameters to be evaluated for the 
# simulation of SCC pattern and returns the parameters to be evaluated in 
# round 10. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC (time-dependent parameters).R")

# Set the directory to store the simulation results. 
save.sims.dir <- "LS_results_r9"
save.sims <- TRUE
if(save.sims) {
  if(!dir.exists(save.sims.dir)) dir.create(save.sims.dir)
}

# Set the cluster. 
n.thread <- detectCores()/2
n.sims <- 50000
cl <- makeCluster(n.thread)
registerDoParallel(cl)

# Read in the initial parameters from the stored .txt file. 
paras.r9 <- as.matrix(read.table("Round 9 parameters.txt", 
                                 sep = "", header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874512, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- generate.pattern(par = paras.r9[i, 1:6], 
                                   init.cells.cols = paras.r9[i, 7],
                                   prob.death = paras.r9[i, 8], 
                                   prob.prof = paras.r9[i, 9],
                                   slopes = paras.r9[i, 10:21])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_9_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff.tot)
}
toc()
stopCluster(cl)

# 70111.131 sec elapsed. 

# Read in the least square differences.
ls.r9 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_9_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff.tot
  ls.r9 <- rbind(ls.r9, ls.temp)
  print(i)
}
ls.r9 <- unname(ls.r9)

# Combine the least square differences with the indices by column. 
ls.r9 <- cbind(seq(1, n.sims, by = 1), ls.r9)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r9[, 2]))
ls.r9.no.nan <- ls.r9[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r9.no.nan[,2]) # 3.977847
# minimum : 45476th 2.494227

# Write the least square differences into a .txt file.  
write.table(ls.r9, "Round 9 Least Square.txt")
# ls.r1 <- as.matrix(read.table("Round 1 Least Square.txt", sep = "",
# header = TRUE))

# Calculate the bandwidth factor which gives the desirable ESS.
info.list.r9 <- calculate.bw(ss.mat = ls.r9, lb.bw = 12.70, ub.bw = 12.80, 
                             ess.target = 10125, step.size = 0.01)

write_rds(info.list.r9, "Round 9 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r10 <- abc_bcd(info.mat = info.list.r9$info.mat, paras.r9)
write.table(paras.r10, "Round 10 parameters.txt")
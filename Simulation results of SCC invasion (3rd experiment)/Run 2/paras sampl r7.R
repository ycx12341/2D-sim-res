# paras sampl r7.R
# Author: Yunchen Xiao
# This .R file reads in the post-round 6 parameters to be evaluated for the 
# simulation of SCC pattern and returns the parameters to be evaluated in 
# round 8. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC (time-dependent parameters).R")

# Set the directory to store the simulation results. 
save.sims.dir <- "LS_results_r7"
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
paras.r7 <- as.matrix(read.table("Round 7 parameters.txt", 
                                 sep = "", header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874512, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- generate.pattern(par = paras.r7[i, 1:6], 
                                   init.cells.cols = paras.r7[i, 7],
                                   prob.death = paras.r7[i, 8], 
                                   prob.prof = paras.r7[i, 9],
                                   slopes = paras.r7[i, 10:21])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_7_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff.tot)
}
toc()
stopCluster(cl)

# 63807.842 sec elapsed. 

# Read in the least square differences.
ls.r7 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_7_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff.tot
  ls.r7 <- rbind(ls.r7, ls.temp)
  print(i)
}
ls.r7 <- unname(ls.r7)

# Combine the least square differences with the indices by column. 
ls.r7 <- cbind(seq(1, n.sims, by = 1), ls.r7)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r7[, 2]))
ls.r7.no.nan <- ls.r7[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r7.no.nan[,2]) # 4.65573
# minimum : 40361th 2.637506

# Write the least square differences into a .txt file.  
write.table(ls.r7, "Round 7 Least Square.txt")
# ls.r1 <- as.matrix(read.table("Round 1 Least Square.txt", sep = "",
# header = TRUE))

# Calculate the bandwidth factor which gives the desirable ESS.
info.list.r7 <- calculate.bw(ss.mat = ls.r7, lb.bw = 9.40, ub.bw = 9.45, 
                             ess.target = 12500, step.size = 0.01)

write_rds(info.list.r7, "Round 7 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r8 <- abc_bcd(info.mat = info.list.r7$info.mat, paras.r7)
write.table(paras.r8, "Round 8 parameters.txt")
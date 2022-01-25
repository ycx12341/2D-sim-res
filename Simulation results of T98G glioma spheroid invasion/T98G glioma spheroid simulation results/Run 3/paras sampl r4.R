# paras sampl r4.R
# Author: Yunchen Xiao
# This .R file reads in the parameters used in round 4 of the simulations T98G
# glioma invasion patterns and returns the parameters to be evaluated in 
# round 5.

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
save.sims.dir <- "LS_results_r4"
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
paras.r4 <- as.matrix(read.table("Round 4 parameters.txt", sep = "", 
                                 header = TRUE))

# Simulation running in parallel. 
registerDoRNG(874514, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r4[i, 1:2], 
                                init.cells.rad = paras.r4[i, 3],
                                prob.death = paras.r4[i, 4],
                                prob.prof = paras.r4[i, 5])
  # Write the results into .rds files, stored in the directory specified 
  # earlier. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_4_paras", i, "_res.rds"))
  
  c(i, sse.res.temp$diff)
}
toc()

stopCluster(cl)

# 5489.378 sec elapsed. 

# Read in the least square differences of the parameter vectors. 
ls.r4 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_4_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r4 <- rbind(ls.r4, ls.temp)
  print(i)
}
ls.r4 <- unname(ls.r4)

# Combine the least square differences with their corresponding indices by 
# column.
ls.r4 <- cbind(seq(1, n.sims, by = 1), ls.r4)

# Take away the singular results (temporarily) and calculate the average 
# least square difference.
ind.nan <- which(is.na(ls.r4[, 2]))
ls.r4.no.nan <- ls.r4[-ind.nan, ]
mean(ls.r4[,2]) # 2.076628

# Locate the minimum one. 
# minimum : 9208th 0.9587063

# Write the least square matrix into a .txt file. 
write.table(ls.r4, "Round 4 Least Square.txt")

# Compute the bandwidth factor that is used to calculate the resampling weights,
# based on the desirable effective sample size (ESS).
info.list.r4 <- calculate.bw(ss.mat = ls.r4, lb.bw = 7.5, ub.bw = 7.7, 
                             ess.target = 1822.5, step.size = 0.01)
write_rds(info.list.r4, "Round 4 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r5 <- abc_bcd(info.mat = info.list.r4$info.mat, paras.r4)
write.table(paras.r5, "Round 5 parameters.txt")

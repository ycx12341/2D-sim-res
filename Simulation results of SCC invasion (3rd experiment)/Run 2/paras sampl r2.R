# paras sampl r2.R
# Author: Yunchen Xiao
# This .R file reads in the post-round 1 parameters to be evaluated for the 
# simulation of SCC pattern and returns the parameters to be evaluated in 
# round 3. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC (time-dependent parameters).R")

# Set the directory to store the simulation results. 
save.sims.dir <- "LS_results_r2"
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
paras.r2 <- as.matrix(read.table("Round 2 parameters.txt", 
                                    sep = "", header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874512, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- generate.pattern(par = paras.r2[i, 1:6], 
                                   init.cells.cols = paras.r2[i, 7],
                                   prob.death = paras.r2[i, 8], 
                                   prob.prof = paras.r2[i, 9],
                                   slopes = paras.r2[i, 10:21])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_2_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 121879.218 sec elapsed. 

# Read in the least square differences.
ls.r2 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_2_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff.tot
  ls.r2 <- rbind(ls.r2, ls.temp)
  print(i)
}
ls.r2 <- unname(ls.r2)

# Combine the least square differences with the indices by column. 
ls.r2 <- cbind(seq(1, n.sims, by = 1), ls.r2)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r2[, 2]))
ls.r2.no.nan <- ls.r2[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r2.no.nan[,2]) # 89.87517
# minimum : 10673th 4.343284

# Write the least square differences into a .txt file.  
write.table(ls.r2, "Round 2 Least Square.txt")
# ls.r2 <- as.matrix(read.table("Round 2 Least Square.txt", sep = "",
# header = TRUE))

# Calculate the bandwidth factor which gives the desirable ESS.
info.list.r2 <- calculate.bw(ss.mat = ls.r2, lb.bw = 0.9, ub.bw = 1.0, 
                             ess.target = 12500, step.size = 0.01)

write_rds(info.list.r2, "Round 2 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r3 <- abc_bcd(info.mat = info.list.r2$info.mat, paras.r2)
write.table(paras.r3, "Round 3 parameters.txt")
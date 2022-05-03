# paras sampl.R
# Author: Yunchen Xiao
# This .R file generates the initial parameters to be evaluated for the 
# simulation of SCC pattern and returns the parameters to be evaluated in 
# round 2. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC (time-dependent parameters).R")

# Set the directory to store the simulation results. 
save.sims.dir <- "LS_results"
save.sims <- TRUE
if(save.sims) {
  if(!dir.exists(save.sims.dir)) dir.create(save.sims.dir)
}

# Set the cluster. 
n.thread <- detectCores() - 1
n.sims <- 50000
cl <- makeCluster(n.thread)
registerDoParallel(cl)


# Read in the initial parameters from the stored .txt file. 
paras.table <- as.matrix(read.table("Round 1 initial time varying parameters.txt", 
                          sep = "", header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874513, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- generate.pattern(par = paras.table[i, 1:6], 
                                init.cells.cols = paras.table[i, 7],
                                prob.death = paras.table[i, 8], 
                                prob.prof = paras.table[i, 9],
                                slopes = paras.table[i, 10:21])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_1_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 36046.264 sec elapsed. 

# Read in the least square differences.
ls.r1 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_1_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff.tot
  ls.r1 <- rbind(ls.r1, ls.temp)
  print(i)
}
ls.r1 <- unname(ls.r1)

# Combine the least square differences with the indices by column. 
ls.r1 <- cbind(seq(1, n.sims, by = 1), ls.r1)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r1[, 2]))
ls.r1.no.nan <- ls.r1[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r1.no.nan[,2]) # 138.2118
# minimum : 27614th 5.204999

# Write the least square differences into a .txt file.  
write.table(ls.r1, "Round 1 Least Square.txt")

ls.r1 <- as.matrix(read.table("Round 1 Least Square.txt", sep = "",
 header = TRUE))

# Calculate the bandwidth factor which gives the desirable ESS. Note that since 
# the parameter vectors which can return valid simulation output is less
# than 12500 in the current round, we chose the desirable ESS to be 
# approximately half of the number of valid output. 
info.list.r1 <- calculate.bw(ss.mat = ls.r1, lb.bw = 0, ub.bw = 0.2, 
                             ess.target = 4550, step.size = 0.01)

write_rds(info.list.r1, "Round 1 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r2 <- abc_bcd(info.mat = info.list.r1$info.mat, paras.table)
write.table(paras.r2, "Round 2 parameters.txt")

# paras sampl r1.R
# Author: Yunchen Xiao
# This .R file generates the initial parameters used in the simulation of T98G
# glioma invasion patterns and returns the parameters to be evaluated in 
# round 2. 

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
save.sims.dir <- "LS_results"
save.sims <- TRUE
if(save.sims) {
  if(!dir.exists(save.sims.dir)) dir.create(save.sims.dir)
}

# Set the cluster.
n.thread <- detectCores() - 1
n.sims <- 10000
cl <- makeCluster(n.thread)
registerDoParallel(cl)

# Generate the initial parameters.
set.seed(874513)
RNGkind(sample.kind = "Rejection")
dn <- runif(n.sims, 0.000069, 0.02)
rn <- runif(n.sims,0,0.0035)
init.cells.rads <- runif(n.sims, 0.01, 0.15)
prob.death <- runif(n.sims, 0.01, 0.1)
prob.prof <- runif(n.sims, 0.2, 1)
paras.table <- cbind(dn, rn, init.cells.rads,
                     prob.death, prob.prof)
write.table(paras.table, "Round 1 initial parameters.txt")

# Simulation running in parallel. 
registerDoRNG(874512, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.table[i, 1:2], 
                                init.cells.rad = paras.table[i, 3],
                                prob.death = paras.table[i, 4], 
                                prob.prof = paras.table[i, 5])
  # Write the results into .rds files, stored in the directory specified 
  # earlier. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_1_paras", i, "_res.rds"))
  
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 6768.979 sec elapsed.

# Read in the least square differences of the parameter vectors. 
ls.r1 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_1_paras",i,"_res.rds"))
  sse.d1.temp <- sum((res.temp$den.mat.d1 - t1.ref.den)^2)
  sse.d3.temp <- sum((res.temp$den.mat.d3 - t3.ref.den)^2)
  ls.temp <- c(res.temp$diff, sse.d1.temp, sse.d3.temp)
  ls.r1 <- rbind(ls.r1, ls.temp)
  print(i)
}
ls.r1 <- unname(ls.r1)

# Combine the least square differences with their corresponding indices by 
# column.
ls.r1 <- cbind(seq(1, n.sims, by = 1), ls.r1)

# Take away the singular results (temporarily) and calculate the average 
# least square difference.
ind.nan <- which(is.na(ls.r1[, 2]))
ls.r1.no.nan <- ls.r1[-ind.nan, ]
mean(ls.r1.no.nan[,2]) # 6.280812

# Locate the minimum one. 
# minimum : 9904th 1.302929

# Write the least square matrix into a .txt file. 
write.table(ls.r1, "Round 1 Least Square.txt")

# Minimum result
ls.r1.min <- read_rds("Round_1_paras9904_res.rds")

# Calculate the SDs of the discrepancy between the best-fitted output
# and the reference data. 
den.mat.d1.r1.min <- ls.r1.min$den.mat.d1
den.mat.d3.r1.min <- ls.r1.min$den.mat.d3
sd.den.mat.d1.r1.min <- sd(den.mat.d1.r1.min - t1.ref.den)
sd.den.mat.d3.r1.min <- sd(den.mat.d3.r1.min - t3.ref.den)

# Compute the bandwidth factor that is used to calculate the resampling weights,
# based on the desirable effective sample size (ESS).
info.list.r1 <- calculate.bw(ss.mat = ls.r1, min.sd.d1 = sd.den.mat.d1.r1.min,
                             min.sd.d3 = sd.den.mat.d3.r1.min, 
                             ess.target = 2500, step.size = 0.0001)
write_rds(info.list.r1, "Round 1 information list log transform.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r2 <- abc_bcd(info.mat = info.list.r1$info.mat, paras.table)
write.table(paras.r2, "Round 2 parameters log transform.txt")

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
source("PDE 2D ABC functions adjusted SCC.R")

# Set the directory to store the simulation results. 
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
gamma <- runif(n.sims, 0.005, 0.26)
rn <- runif(n.sims, 0.0008, 0.08)
eta <- runif(n.sims, 7,18)
dm <- runif(n.sims,0.0001,0.033)
alpha <- runif(n.sims,0.07,0.18)
prob.death <- runif(n.sims, 0.01, 0.1)
prob.prof <- runif(n.sims, 0.2, 1)
paras.table <- cbind(dn, gamma, rn, eta, dm, alpha,
                     prob.death, prob.prof)

# Write the initial parameters into a .txt file. 
write.table(paras.table, "Round 1 initial parameters.txt")

# Run the simulation in parallel. 
registerDoRNG(874514, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests_ver2 <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.table[i, 1:6],
                                prob.death = paras.table[i, 7], 
                                prob.prof = paras.table[i, 8])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_1_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()

stopCluster(cl)

# 20554.82 sec elapsed.

# Read in the least square differences.
ls.r1 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_1_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r1 <- rbind(ls.r1, ls.temp)
}
ls.r1 <- unname(ls.r1)

# Combine the least square differences with the indices by column. 
ls.r1 <- cbind(seq(1, n.sims, by = 1), ls.r1)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r1[, 2]))
ls.r1.no.nan <- ls.r1[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r1.no.nan[,2]) # 7.36273
# minimum : 2357th 0.4085451

# Write the least square differences into a .txt file.  
write.table(ls.r1, "Round 1 Least Square.txt")

# Calculate the bandwidth factor which gives the desirable ESS.
info.list.r1 <- calculate.bw(ss.mat = ls.r1, lb.bw = 1.2, ub.bw = 1.5, 
                             ess.target = 2500, step.size = 0.01)
write_rds(info.list.r1, "Round 1 information list.rds")

# Resample and record the parameters to be evaluated in the next round.
set.seed(123)
RNGkind(sample.kind = "Rejection")
paras.r2 <- abc_bcd(info.mat = info.list.r1$info.mat, paras.table)
write.table(paras.r2, "Round 2 parameters.txt")

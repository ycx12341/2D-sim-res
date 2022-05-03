# paras sampl r5.R
# Author: Yunchen Xiao
# This .R file reads in the parameter values to be evaluated for the 
# simulation of SCC pattern in the current round and returns final simulation 
# output for run 2 of post-day 14 SCC pattern. 

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
registerDoRNG(874514, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r5[i, 1:6],
                                prob.death = paras.r5[i,7], 
                                prob.prof = paras.r5[i,8])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_5_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 20904.03 sec elapsed. 

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
mean(ls.r5.no.nan[,2]) # 0.7772115
# minimum : 4566th 0.3956359

# Write the least square differences into a .txt file.  
write.table(ls.r5, "Round 5 Least Square.txt")
# ls.r5 <- as.matrix(read.table("Round 5 Least Square.txt", sep = "",
#                               header = TRUE))

# Reduction in average least square difference is less than 5%, simulation for 
# the current run can be terminated. 
paras.r5.mean <- apply(paras.r5, 2, mean)

# Generate the final simulation output. 
set.seed(874514)
RNGkind(sample.kind = "Rejection")
paras.r5.mean.res <- calculate.sse(pars = paras.r5.mean[1:6], 
                                   prob.death = paras.r5.mean[7], 
                                   prob.prof = paras.r5.mean[8])

# Write the final output into an .rds file. 
write_rds(paras.r5.mean.res, "paras r5 average result.rds")
# paras.r5.mean.res <- read_rds("paras r5 average result.rds")

# Plot of individual cell positions. 
sse.t14.r5.mean.grid <- paras.r5.mean.res$ind.position.d14
cell.coord <- which(sse.t14.r5.mean.grid == 1, arr.ind = TRUE)
cell.coord.x <- cell.coord[,1]
cell.coord.y <- cell.coord[,2]
space.length.x <- 35
space.length.y <- 60
x <- seq(1, space.length.x, by = 1)
y <- seq(1, space.length.y, by = 1)
plot(x[cell.coord.y[1]], y[space.length.y + 1 - cell.coord.x[1]], pch = 15,
     bg = "black", xlim = c(1, space.length.x), ylim = c(1, space.length.y), 
     cex = 1.2, xaxt = 'n', yaxt = 'n')

for (i in 2:length(cell.coord.x)) {
  points(x[cell.coord.y[i]], y[space.length.y + 1 - cell.coord.x[i]],
         pch = 15, cex = 1.2, bg = "black")
}

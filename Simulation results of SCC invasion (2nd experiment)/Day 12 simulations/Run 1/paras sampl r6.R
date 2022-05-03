# paras sampl r6.R
# Author: Yunchen Xiao
# This .R file reads in the parameter values to be evaluated for the 
# simulation of SCC pattern in the current round and returns final simulation 
# output for run 1 of post-day 12 SCC pattern. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC.R")

# Set the directory to store the simulation results. 
save.sims.dir <- "LS_results_r6"
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
paras.r6 <- as.matrix(read.table("Round 6 parameters.txt", sep = "",
                                 header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874513, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r6[i, 1:6],
                                prob.death = paras.r6[i, 7], 
                                prob.prof = paras.r6[i, 8])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_6_paras", i, "_res.rds"))
  
  c(i, sse.res.temp$diff)
}
toc()

stopCluster(cl)

# 28030.05 sec elapsed. 

# Read in the least square differences.
ls.r6 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_6_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r6 <- rbind(ls.r6, ls.temp)
  print(i)
}
ls.r6 <- unname(ls.r6)

# Combine the least square differences with the indices by column. 
ls.r6 <- cbind(seq(1, n.sims, by = 1), ls.r6)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r6[, 2]))
ls.r6.no.nan <- ls.r6[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r6.no.nan[,2]) # 0.4805667
# minimum : 9282th 0.2348165 

# Write the least square differences into a .txt file.  
write.table(ls.r6, "Round 6 Least Square.txt")
# ls.r6 <- as.matrix(read.table("Round 6 Least Square.txt", sep = "",
#                               header = TRUE))

# Reduction in average least square difference is less than 5%, simulation for 
# the current run can be terminated. 
paras.r6.mean <- apply(paras.r6, 2, mean)

# Generate the final simulation output. 
set.seed(874513)
RNGkind(sample.kind = "Rejection")
paras.r6.mean.res <- calculate.sse(pars = paras.r6.mean[1:6], 
                                   prob.death = paras.r6.mean[7], 
                                   prob.prof = paras.r6.mean[8])

# Write the final output into an .rds file. 
write_rds(paras.r6.mean.res, "paras r6 average result.rds")
# paras.r6.mean.res <- read_rds("paras r6 average result.rds")

# Plot of individual cell positions. 
sse.t12.r6.mean.grid <- paras.r6.mean.res$ind.position.d12
cell.coord <- which(sse.t12.r6.mean.grid == 1, arr.ind = TRUE)
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
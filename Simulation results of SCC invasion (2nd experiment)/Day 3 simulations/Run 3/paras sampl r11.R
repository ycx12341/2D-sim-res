# paras sampl r11.R
# Author: Yunchen Xiao
# This .R file reads in the parameter values to be evaluated for the 
# simulation of SCC pattern in the current round and returns final simulation 
# output for run 3 of post-day 3 SCC pattern. 

# Clear the workspace and load the necessary packages.
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC.R")

# Set the directory to store the simulation results. 
save.sims.dir <- "LS_results_r11"
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
paras.r11 <- as.matrix(read.table("Round 11 parameters.txt", sep = "",
                                  header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874512, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r11[i, 1:6], 
                                init.cells.cols = paras.r11[i, 7],
                                prob.death = paras.r11[i, 8], 
                                prob.prof = paras.r11[i, 9])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_11_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 5007.86 sec elapsed. 

# Read in the least square differences.
ls.r11 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_11_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r11 <- rbind(ls.r11, ls.temp)
  print(i)
}
ls.r11 <- unname(ls.r11)

# Combine the least square differences with the indices by column. 
ls.r11 <- cbind(seq(1, n.sims, by = 1), ls.r11)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r11[, 2]))
ls.r11.no.nan <- ls.r11[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r11.no.nan[,2]) # 0.2087382
# minimum : 9900th 0.08170499

# Write the least square differences into a .txt file.  
write.table(ls.r11, "Round 11 Least Square.txt")
# ls.r11 <- as.matrix(read.table("Round 11 Least Square.txt", sep = "",
#                               header = TRUE))

# Reduction in average least square difference is less than 5%, simulation for 
# the current run can be terminated. 
paras.r11.mean <- apply(paras.r11, 2, mean)

# Generate the final simulation output. 
set.seed(874512)
RNGkind(sample.kind = "Rejection")
paras.r11.mean.res <- calculate.sse(pars = paras.r11.mean[1:6], 
                                    init.cells.cols = paras.r11.mean[7],
                                    prob.death = paras.r11.mean[8], 
                                    prob.prof = paras.r11.mean[9])
# Write the final output into an .rds file. 
write_rds(paras.r11.mean.res, "paras r11 average result.rds")

# Plot of individual cell positions. 
sse.t3.r11.mean.grid <- paras.r11.mean.res$ind.position.d3
cell.coord <- which(sse.t3.r11.mean.grid == 1, arr.ind = TRUE)
cell.coord.x <- cell.coord[,1]
cell.coord.y <- cell.coord[,2]
space.length.x <- 35
space.length.y <- 60
x <- seq(1, space.length.x, by = 1)
y <- seq(1, space.length.y, by = 1)
plot(x[cell.coord.y[1]], y[space.length.y + 1 - cell.coord.x[1]], pch = 15,
     bg = "black", xlim = c(1, space.length.x), ylim = c(1, space.length.y), cex = 1.2, xaxt = 'n', yaxt = 'n')

for (i in 2:length(cell.coord.x)) {
  points(x[cell.coord.y[i]], y[space.length.y + 1 - cell.coord.x[i]],
         pch = 15, cex = 1.2, bg = "black")
}
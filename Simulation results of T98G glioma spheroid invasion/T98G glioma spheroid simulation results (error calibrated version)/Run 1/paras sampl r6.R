# paras sampl r6.R
# Author: Yunchen Xiao
# This .R file reads in the parameters used in round 6 of the simulations T98G
# glioma invasion patterns and returns the final simulated patterns based on 
# the average parameter estimates. 

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
save.sims.dir <- "LS_results_r6"
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
paras.r6 <- as.matrix(read.table("Round 6 parameters log transform.txt", 
                                 sep = "", header = TRUE))

# Simulation running in parallel. 
registerDoRNG(874513, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r6[i, 1:2], 
                                init.cells.rad = paras.r6[i, 3],
                                prob.death = paras.r6[i, 4], 
                                prob.prof = paras.r6[i, 5])
  # Write the results into .rds files, stored in the directory specified 
  # earlier. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_6_paras", i, "_res.rds"))
  
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 5899.08 sec elapsed. 

# Read in the least square differences of the parameter vectors. 
ls.r6 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_6_paras",i,"_res.rds"))
  ls.temp <- c(res.temp$diff, res.temp$sse.t1, res.temp$sse.t3)
  ls.r6 <- rbind(ls.r6, ls.temp)
  print(i)
}
ls.r6 <- unname(ls.r6)

# Combine the least square differences with their corresponding indices by 
# column. 
ls.r6 <- cbind(seq(1, n.sims, by = 1), ls.r6)

# Check if there are singular results in the simulation output and calculate 
# the average least square difference.
ind.nan <- which(is.na(ls.r6[, 2]))
mean(ls.r6[,2]) # 1.54459

# Minimum result
# minimum : 5694th 0.9290988

# Write the least square matrix into a .txt file. 
write.table(ls.r6, "Round 6 Least Square.txt")

# Reduction in average least square difference is lower than 5%, algorithm can 
# be terminated now! (4.87%)
(1.623746 - 1.54459)/1.623746*100  

# Use the averaged final parameter estimates to simulate the pattern.
paras.r6.mean <- apply(paras.r6, 2, mean)
set.seed(874513)
RNGkind(sample.kind = "Rejection")
paras.r6.mean.res <- calculate.sse(pars = paras.r6.mean[1:2],
                                   init.cells.rad = paras.r6.mean[3],
                                   prob.death = paras.r6.mean[4], 
                                   prob.prof = paras.r6.mean[5])

# Write the final simulated results into a .rds file. 
write_rds(paras.r6.mean.res, "paras r6 average result.rds")
# paras.r6.mean.res <- read_rds("paras r6 average result.rds")


# Read in the inidividual cell positions in the final simulated pattern. 
sse.t1.r6.mean.grid <- paras.r6.mean.res$ind.position.d1
sse.t3.r6.mean.grid <- paras.r6.mean.res$ind.position.d3

# Individual cells plot. 
# cell.coord <- which(sse.t1.r6.mean.grid == 1, arr.ind = TRUE)
cell.coord <- which(sse.t3.r6.mean.grid == 1, arr.ind = TRUE)
cell.coord.x <- cell.coord[,1]
cell.coord.y <- cell.coord[,2]
space.length.x <- 40
space.length.y <- 40
x <- seq(1, space.length.x, by = 1)
y <- seq(1, space.length.y, by = 1)
plot(x[cell.coord.y[1]], y[space.length.y + 1 - cell.coord.x[1]], pch = 15,
     bg = "black", xlim = c(1, space.length.x), ylim = c(1, space.length.y), 
     cex = 1.2, xaxt = 'n', yaxt = 'n')

for (i in 2:length(cell.coord.x)) {
  points(x[cell.coord.y[i]], y[space.length.y + 1 - cell.coord.x[i]],
         pch = 15, cex = 1.2, bg = "black")
}

# Grey shades plots
paras.r6.mean.res <- read_rds("paras r6 average result.rds")
rotate <- function(x) t(apply(x, 2, rev))
t1.ref.den.rotate <- rotate(t1.ref.den)
t3.ref.den.rotate <- rotate(t3.ref.den)
den.mat.d1 <- paras.r6.mean.res$den.mat.d1
den.mat.d3 <- paras.r6.mean.res$den.mat.d3
den.mat.d1.rotate <- rotate(den.mat.d1)
den.mat.d3.rotate <- rotate(den.mat.d3)
grey.spectrum.full <- colorRampPalette(c("grey100", "grey0"))
grey.spectrum.50 <- grey.spectrum.full(50)
image(t1.ref.den.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(den.mat.d1.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(t3.ref.den.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(den.mat.d3.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')

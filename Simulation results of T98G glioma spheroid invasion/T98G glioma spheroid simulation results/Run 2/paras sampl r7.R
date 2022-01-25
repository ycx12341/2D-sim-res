# paras sampl r7.R
# Author: Yunchen Xiao
# This .R file reads in the parameters used in round 7 of the simulations T98G
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
save.sims.dir <- "LS_results_r7"
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
paras.r7 <- as.matrix(read.table("Round 7 parameters.txt", sep = "", 
                                 header = TRUE))

# Simulation running in parallel. 
registerDoRNG(874512, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- calculate.sse(pars = paras.r7[i, 1:2], 
                                init.cells.rad = paras.r7[i, 3],
                                prob.death = paras.r7[i, 4], 
                                prob.prof = paras.r7[i, 5])
  # Write the results into .rds files, stored in the directory specified 
  # earlier. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_7_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff)
}
toc()
stopCluster(cl)

# 4654.247 sec elapsed. 

# Read in the least square differences of the parameter vectors. 
ls.r7 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_7_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff
  ls.r7 <- rbind(ls.r7, ls.temp)
  print(i)
}
ls.r7 <- unname(ls.r7)

# Combine the least square differences with their corresponding indices by 
# column.
ls.r7 <- cbind(seq(1, n.sims, by = 1), ls.r7)

# Check if the simulation results still contain singular values. 
ind.nan <- which(is.na(ls.r7[, 2]))
# Calculate the averaged least square difference. 
mean(ls.r7[,2]) # 1.506051

# Reduction in average least square difference is lower than 5%, algorithm can 
# be terminated now! (2.27%)
(1.5461068 - 1.506051)/1.5461068*100  

# Locate the minimum one. 
# minimum : 425th 0.8174601

# Write the least square matrix into a .txt file. 
write.table(ls.r7, "Round 7 Least Square.txt")

# Use the averaged final parameter estimates to simulate the pattern.
paras.r7.mean <- apply(paras.r7, 2, mean)
set.seed(874512)
RNGkind(sample.kind = "Rejection")
paras.r7.mean.res <- calculate.sse(pars = paras.r7.mean[1:2],
                                   init.cells.rad = paras.r7.mean[3],
                                   prob.death = paras.r7.mean[4], 
                                   prob.prof = paras.r7.mean[5])

# Write the final simulated results into a .rds file. 
write_rds(paras.r7.mean.res, "paras r7 average result.rds")
# paras.r7.mean.res <- read_rds("paras r7 average result.rds")

# Read in the inidividual cell positions in the final simulated pattern. 
sse.t1.r7.mean.grid <- paras.r7.mean.res$ind.position.d1
sse.t3.r7.mean.grid <- paras.r7.mean.res$ind.position.d3

# Individual cells plot. 
# cell.coord <- which(sse.t1.r7.mean.grid == 1, arr.ind = TRUE)
cell.coord <- which(sse.t3.r7.mean.grid == 1, arr.ind = TRUE)
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
paras.r7.mean.res <- read_rds("paras r7 average result.rds")
rotate <- function(x) t(apply(x, 2, rev))
t1.ref.den.rotate <- rotate(t1.ref.den)
t3.ref.den.rotate <- rotate(t3.ref.den)
den.mat.d1 <- paras.r7.mean.res$den.mat.d1
den.mat.d3 <- paras.r7.mean.res$den.mat.d3
den.mat.d1.rotate <- rotate(den.mat.d1)
den.mat.d3.rotate <- rotate(den.mat.d3)
grey.spectrum.full <- colorRampPalette(c("grey100", "grey0"))
grey.spectrum.50 <- grey.spectrum.full(50)
image(t1.ref.den.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(den.mat.d1.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(t3.ref.den.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(den.mat.d3.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')

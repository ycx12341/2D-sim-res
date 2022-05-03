# paras sampl r8.R
# Author: Yunchen Xiao
# This .R file reads in the post-round 7 parameters to be evaluated for the 
# simulation of SCC pattern and returns the parameters to be evaluated in 
# round 9. 

# Clear the workspace and load the necessary packages. 
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)

# Load the file that contains the functions. 
source("PDE 2D ABC functions adjusted SCC (time-dependent parameters).R")

# Set the directory to store the simulation results. 
save.sims.dir <- "LS_results_r8"
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
paras.r8 <- as.matrix(read.table("Round 8 parameters.txt", 
                                 sep = "", header = TRUE))

# Run the simulation in parallel. 
registerDoRNG(874513, once = TRUE)
RNGkind(sample.kind = "Rejection")
tic()
ests <- foreach (i = 1:n.sims, .combine = rbind) %dopar% {
  sse.res.temp <- generate.pattern(par = paras.r8[i, 1:6], 
                                   init.cells.cols = paras.r8[i, 7],
                                   prob.death = paras.r8[i, 8], 
                                   prob.prof = paras.r8[i, 9],
                                   slopes = paras.r8[i, 10:21])
  # Write the simulation results into .rds files. 
  readr::write_rds(sse.res.temp, 
                   path = paste0("./", save.sims.dir, 
                                 "/Round_8_paras", i, "_res.rds"))
  c(i, sse.res.temp$diff.tot)
}
toc()
stopCluster(cl)

# 70464.788 sec elapsed. 

# Read in the least square differences.
ls.r8 <- vector()
for (i in 1:n.sims) {
  res.temp <- read_rds(paste0("./", save.sims.dir, 
                              "/Round_8_paras",i,"_res.rds"))
  ls.temp <- res.temp$diff.tot
  ls.r8 <- rbind(ls.r8, ls.temp)
  print(i)
}
ls.r8 <- unname(ls.r8)

# Combine the least square differences with the indices by column. 
ls.r8 <- cbind(seq(1, n.sims, by = 1), ls.r8)

# Locate the singular values and take them away temporarily.
ind.nan <- which(is.na(ls.r8[, 2]))
ls.r8.no.nan <- ls.r8[-ind.nan, ]

# Calculate the average least square differences and locate the minimum one. 
mean(ls.r8.no.nan[,2]) # 4.459617
# minimum : 28248th 3.145709

# Write the least square differences into a .txt file.  
write.table(ls.r8, "Round 8 Least Square.txt")
# ls.r1 <- as.matrix(read.table("Round 1 Least Square.txt", sep = "",
# header = TRUE))

# Reduction in average least square difference is less than 5%, simulation for 
# the current run can be terminated. 
paras.r8.mean <- apply(paras.r8, 2, mean)

# Generate the final simulation output. 
set.seed(874513)
RNGkind(sample.kind = "Rejection")
paras.r8.mean.res <- generate.pattern(par = paras.r8.mean[1:6], 
                                      init.cells.cols = paras.r8.mean[7],
                                      prob.death = paras.r8.mean[8], 
                                      prob.prof = paras.r8.mean[9],
                                      slopes = paras.r8.mean[10:21])

# Write the final output into an .rds file. 
write_rds(paras.r8.mean.res, "paras r8 average result.rds")

# Final results plots
paras.r8.mean.res <- read_rds("paras r8 average result.rds")

# Plot of individual cell positions. 
sse.t3.r8.mean.grid <- paras.r8.mean.res$ind.position.d3
sse.t6.r8.mean.grid <- paras.r8.mean.res$ind.position.d6
sse.t9.r8.mean.grid <- paras.r8.mean.res$ind.position.d9
sse.t12.r8.mean.grid <- paras.r8.mean.res$ind.position.d12
sse.t14.r8.mean.grid <- paras.r8.mean.res$ind.position.d14

cell.coord <- which(sse.t14.r8.mean.grid == 1, arr.ind = TRUE)
cell.coord.x <- cell.coord[,1]
cell.coord.y <- cell.coord[,2]
space.length.x <- 28
space.length.y <- 48
x <- seq(1, space.length.x, by = 1)
y <- seq(1, space.length.y, by = 1)
plot(x[cell.coord.y[1]], y[space.length.y + 1 - cell.coord.x[1]], pch = 15,
     bg = "black", xlim = c(1, space.length.x), ylim = c(1, space.length.y), 
     cex = 1.2, xaxt = 'n', yaxt = 'n')

for (i in 2:length(cell.coord.x)) {
  points(x[cell.coord.y[i]], y[space.length.y + 1 - cell.coord.x[i]],
         pch = 15, cex = 1.2, bg = "black")
}

# Grey shades density plot. 
rotate <- function(x) t(apply(x, 2, rev))
grey.spectrum.full <- colorRampPalette(c("grey100", "grey0"))
grey.spectrum.50 <- grey.spectrum.full(50)
den.mat.d3 <- paras.r8.mean.res$den.mat.d3
den.mat.d6 <- paras.r8.mean.res$den.mat.d6
den.mat.d9 <- paras.r8.mean.res$den.mat.d9
den.mat.d12 <- paras.r8.mean.res$den.mat.d12
den.mat.d14 <- paras.r8.mean.res$den.mat.d14

t3.ref.den.rotate <- rotate(t3.ref.den)
den.mat.d3.rotate <- rotate(den.mat.d3)

t6.ref.den.rotate <- rotate(t6.ref.den)
den.mat.d6.rotate <- rotate(den.mat.d6)

t9.ref.den.rotate <- rotate(t9.ref.den)
den.mat.d9.rotate <- rotate(den.mat.d9)

t12.ref.den.rotate <- rotate(t12.ref.den)
den.mat.d12.rotate <- rotate(den.mat.d12)

t14.ref.den.rotate <- rotate(t14.ref.den)
den.mat.d14.rotate <- rotate(den.mat.d14)

image(t3.ref.den.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(den.mat.d3.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')

image(t6.ref.den.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(den.mat.d6.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')

image(t9.ref.den.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(den.mat.d9.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')

image(t12.ref.den.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(den.mat.d12.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')

image(t14.ref.den.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(den.mat.d14.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')

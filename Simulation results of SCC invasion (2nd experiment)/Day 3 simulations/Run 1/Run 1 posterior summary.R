# Run 1 posterior summary.R
# This .R file reads in the related information of the least square differences
# obtained at the end of each round in run 1, extracts the actual ESS and 
# corresponding bandwidth factors and combine them into a matrix. Additionally, 
# the grey shades density plot of the final output in run 1 is generated in 
# this file. 

# Clear the workspace and load the necessary package. 
rm(list = ls())
library(readr)

# Read in the information list regarding the least square differences obtained
# at the end of each round. 
info.list.r1 <- read_rds("Round 1 information list.rds")
info.list.r2 <- read_rds("Round 2 information list.rds")
info.list.r3 <- read_rds("Round 3 information list.rds")
info.list.r4 <- read_rds("Round 4 information list.rds")
info.list.r5 <- read_rds("Round 5 information list.rds")
info.list.r6 <- read_rds("Round 6 information list.rds")
info.list.r7 <- read_rds("Round 7 information list.rds")
info.list.r8 <- read_rds("Round 8 information list.rds")
info.list.r9 <- read_rds("Round 9 information list.rds")
info.list.r10 <- read_rds("Round 10 information list.rds")

# Extract the actual ESS and corresponding bandwidth factors, combine them
# into a matrix by column and record it. 
ess.vec <- c(info.list.r1$ess.obj, info.list.r2$ess.obj, info.list.r3$ess.obj,
             info.list.r4$ess.obj, info.list.r5$ess.obj, info.list.r6$ess.obj,
             info.list.r7$ess.obj, info.list.r8$ess.obj, info.list.r9$ess.obj,
             info.list.r10$ess.obj)
bw.vec <- c(info.list.r1$bw.obj, info.list.r2$bw.obj, info.list.r3$bw.obj,
             info.list.r4$bw.obj, info.list.r5$bw.obj, info.list.r6$bw.obj,
             info.list.r7$bw.obj, info.list.r8$bw.obj, info.list.r9$bw.obj,
             info.list.r10$bw.obj)
info.mat.run1 <- unname(cbind(ess.vec, bw.vec))
write.table(info.mat.run1, "Run 1 ESS BW.txt")

# Grey shades density plot. 
source("PDE 2D ABC functions adjusted SCC.R")
paras.r11.mean.res <- read_rds("paras r11 average result.rds")
rotate <- function(x) t(apply(x, 2, rev))
grey.spectrum.full <- colorRampPalette(c("grey100", "grey0"))
grey.spectrum.50 <- grey.spectrum.full(50)
den.mat.d3 <- paras.r11.mean.res$den.mat.d3
t3.ref.den.rotate <- rotate(t3.ref.den)
den.mat.d3.rotate <- rotate(den.mat.d3)
image(t3.ref.den.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(den.mat.d3.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')

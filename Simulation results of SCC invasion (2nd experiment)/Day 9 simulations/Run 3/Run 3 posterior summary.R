# Run 3 posterior summary.R
# This .R file reads in the related information of the least square differences
# obtained at the end of each round in run 3, extracts the actual ESS and 
# corresponding bandwidth factors and combine them into a matrix. Additionally, 
# the grey shades density plot of the final output in run 3 is generated in 
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

# Extract the actual ESS and corresponding bandwidth factors, combine them
# into a matrix by column and record it.
ess.vec <- c(info.list.r1$ess.obj, info.list.r2$ess.obj, info.list.r3$ess.obj,
             info.list.r4$ess.obj)
bw.vec <- c(info.list.r1$bw.obj, info.list.r2$bw.obj, info.list.r3$bw.obj,
            info.list.r4$bw.obj)
info.mat.run3 <- unname(cbind(ess.vec, bw.vec))
write.table(info.mat.run3, "Run 3 ESS BW.txt")

# Grey shades density plot. 
source("PDE 2D ABC functions adjusted SCC.R")
paras.r5.mean.res <- read_rds("paras r5 average result.rds")
rotate <- function(x) t(apply(x, 2, rev))
grey.spectrum.full <- colorRampPalette(c("grey100", "grey0"))
grey.spectrum.50 <- grey.spectrum.full(50)
den.mat.d9 <- paras.r5.mean.res$den.mat.d9
t9.ref.den.rotate <- rotate(t9.ref.den)
den.mat.d9.rotate <- rotate(den.mat.d9)
image(t9.ref.den.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
image(den.mat.d9.rotate, col = grey.spectrum.50, xaxt = 'n', yaxt = 'n')
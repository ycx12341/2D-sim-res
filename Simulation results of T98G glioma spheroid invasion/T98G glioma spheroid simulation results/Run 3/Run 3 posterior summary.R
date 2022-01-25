# Run 3 posterior summary.R
# Author: Yunchen Xiao
# This .R file reads in the information lists which contain the related 
# calculations on the least square differences obtained at the end of each 
# round in run 3 (e.g. actual ESS, bandwidth factor, weights of simulation 
# results, resampling probabilities), extracts the value of actual ESS, 
# bandwidth factors then form them into a matrix. 

# Set the workspace, then load the necessary packages. 
rm(list = ls())
library(readr)

# Read in the information lists. 
info.list.r1 <- read_rds("Round 1 information list.rds")
info.list.r2 <- read_rds("Round 2 information list.rds")
info.list.r3 <- read_rds("Round 3 information list.rds")
info.list.r4 <- read_rds("Round 4 information list.rds")
info.list.r5 <- read_rds("Round 5 information list.rds")
info.list.r6 <- read_rds("Round 6 information list.rds")

# Read in the actual ESS and the corresponding bandwidth factor. 
ess.vec <- c(info.list.r1$ess.obj, info.list.r2$ess.obj, 
             info.list.r3$ess.obj, info.list.r4$ess.obj, 
             info.list.r5$ess.obj, info.list.r6$ess.obj)

bw.vec <- c(info.list.r1$bw.obj, info.list.r2$bw.obj, 
            info.list.r3$bw.obj, info.list.r4$bw.obj, 
            info.list.r5$bw.obj, info.list.r6$bw.obj)

# Combines these values by column then write them into a .txt file. 
info.mat.run3 <- unname(cbind(ess.vec, bw.vec))
write.table(info.mat.run3, "Run 3 ESS BW.txt")

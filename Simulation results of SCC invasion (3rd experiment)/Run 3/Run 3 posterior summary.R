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
info.list.r5 <- read_rds("Round 5 information list.rds")
info.list.r6 <- read_rds("Round 6 information list.rds")
info.list.r7 <- read_rds("Round 7 information list.rds")
info.list.r8 <- read_rds("Round 8 information list.rds")

# Extract the actual ESS and corresponding bandwidth factors, combine them
# into a matrix by column and record it. 
ess.vec <- c(info.list.r1$ess.obj, info.list.r2$ess.obj, info.list.r3$ess.obj,
             info.list.r4$ess.obj, info.list.r5$ess.obj, info.list.r6$ess.obj,
             info.list.r7$ess.obj, info.list.r8$ess.obj)

bw.vec <- c(info.list.r1$bw.obj, info.list.r2$bw.obj, info.list.r3$bw.obj,
            info.list.r4$bw.obj, info.list.r5$bw.obj, info.list.r6$bw.obj,
            info.list.r7$bw.obj, info.list.r8$bw.obj)

info.mat.run3 <- unname(cbind(ess.vec, bw.vec))
write.table(info.mat.run3, "Run 3 ESS BW.txt")

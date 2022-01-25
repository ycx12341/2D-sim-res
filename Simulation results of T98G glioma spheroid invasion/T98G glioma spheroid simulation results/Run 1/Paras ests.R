# Paras ests.R
# Author: Yunchen Xiao
# This .R file reads in the parameter values being evaluated in each round
# and calculates each parameter's average estimate. 

# Clear the workspace. 
rm(list = ls())

# Reads in the parameter values being evaluated in each round. 
paras.init.run1 <- read.table("Round 1 initial parameters.txt", sep = "", 
                         header = TRUE)
paras.r2.run1 <- read.table("Round 2 parameters.txt", sep = "", 
                         header = TRUE)
paras.r3.run1 <- read.table("Round 3 parameters.txt", sep = "", 
                         header = TRUE)
paras.r4.run1 <- read.table("Round 4 parameters.txt", sep = "", 
                         header = TRUE)
paras.r5.run1 <- read.table("Round 5 parameters.txt", sep = "", 
                         header = TRUE)
paras.r6.run1 <- read.table("Round 6 parameters.txt", sep = "", 
                         header = TRUE)
paras.r7.run1 <- read.table("Round 7 parameters.txt", sep = "", 
                         header = TRUE)

# Computes the average parameter estimates. 
paras.init.mean <- apply(paras.init.run1, 2, mean)
paras.r2.mean <- apply(paras.r2.run1, 2, mean)
paras.r3.mean <- apply(paras.r3.run1, 2, mean)
paras.r4.mean <- apply(paras.r4.run1, 2, mean)
paras.r5.mean <- apply(paras.r5.run1, 2, mean)
paras.r6.mean <- apply(paras.r6.run1, 2, mean)
paras.r7.mean <- apply(paras.r7.run1, 2, mean)

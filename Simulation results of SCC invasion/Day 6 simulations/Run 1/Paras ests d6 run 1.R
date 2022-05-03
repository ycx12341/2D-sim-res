# Paras est d6 run 1. R
# Author: Yunchen Xiao
# This .R file reads in the parameter values being evaluated in each round
# and calculates each parameter's average estimate. 

# Clear the workspace. 
rm(list = ls())

# Reads in the parameter values being evaluated in each round. 
paras.r1.run1 <- as.matrix(read.table("Round 1 initial parameters.txt", 
                                      sep = "",
                                      header = TRUE))
paras.r2.run1 <- as.matrix(read.table("Round 2 parameters.txt", sep = "",
                                      header = TRUE))
paras.r3.run1 <- as.matrix(read.table("Round 3 parameters.txt", sep = "",
                                      header = TRUE))
paras.r4.run1 <- as.matrix(read.table("Round 4 parameters.txt", sep = "",
                                      header = TRUE))
paras.r5.run1 <- as.matrix(read.table("Round 5 parameters.txt", sep = "",
                                      header = TRUE))
paras.r6.run1 <- as.matrix(read.table("Round 6 parameters.txt", sep = "",
                                      header = TRUE))

# Computes the average parameter estimates. 
paras.r1.run1.mean <- apply(paras.r1.run1, 2, mean)
paras.r2.run1.mean <- apply(paras.r2.run1, 2, mean)
paras.r3.run1.mean <- apply(paras.r3.run1, 2, mean)
paras.r4.run1.mean <- apply(paras.r4.run1, 2, mean)
paras.r5.run1.mean <- apply(paras.r5.run1, 2, mean)
paras.r6.run1.mean <- apply(paras.r6.run1, 2, mean)

paras.run1.mean <- rbind(paras.r1.run1.mean, paras.r2.run1.mean, 
                         paras.r3.run1.mean, paras.r4.run1.mean, 
                         paras.r5.run1.mean, paras.r6.run1.mean)

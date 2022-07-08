# Overlap SCC d14.R
# Author: Yunchen Xiao
# This .R file calculates the prior-posterior overlap for the second experiment
# of the simulation study regarding the SCC invasion patterns. 
# Post day 14 pattern.

# Clear the current workspace and load the necessary packages.
rm(list = ls())
library("overlapping")

# Read the prior and posterior parameters.
paras.r1 <- as.matrix(read.table("Round 1 initial parameters.txt", sep = "",
                                 header = TRUE))
paras.r5 <- as.matrix(read.table("Round 5 parameters.txt", sep = "", 
                                 header = TRUE))

# Prior-posterior overlap calculations.
dn.pri.post <- list(dn.prior = paras.r1[,1], dn.post = paras.r5[,1])
dn.overlap <- overlap(dn.pri.post)
# 25.36%

gamma.pri.post <- list(gamma.prior = paras.r1[,2], gamma.post = paras.r5[,2])
gamma.overlap <- overlap(gamma.pri.post)
# 83.05%

rn.pri.post <- list(rn.prior = paras.r1[,3], rn.post = paras.r5[,3])
rn.overlap <- overlap(rn.pri.post)
# 87.30%

eta.pri.post <- list(eta.prior = paras.r1[,4], eta.post = paras.r5[,4])
eta.overlap <- overlap(eta.pri.post)
# 86.48%

dm.pri.post <- list(dm.prior = paras.r1[,5], dm.post = paras.r5[,5])
dm.overlap <- overlap(dm.pri.post)
# 76.84%

alpha.pri.post <- list(alpha.prior = paras.r1[,6], r.init.post = paras.r5[,6])
alpha.overlap <- overlap(alpha.pri.post)
# 86.87%

p.ext.pri.post <- list(p.ext.prior = paras.r1[,7], p.ext.post = paras.r5[,7])
p.ext.overlap <- overlap(p.ext.pri.post)
# 85.23%

p.mit.pri.post <- list(p.mit.prior = paras.r1[,8], p.mit.post = paras.r5[,8])
p.mit.overlap <- overlap(p.mit.pri.post)
# 19.49%
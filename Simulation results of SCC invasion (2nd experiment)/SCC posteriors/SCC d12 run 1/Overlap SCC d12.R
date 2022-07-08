# Overlap SCC d12.R
# Author: Yunchen Xiao
# This .R file calculates the prior-posterior overlap for the second experiment
# of the simulation study regarding the SCC invasion patterns. 
# Post day 12 pattern.

# Clear the current workspace and load the necessary packages.
rm(list = ls())
library("overlapping")

# Read the prior and posterior parameters.
paras.r1 <- as.matrix(read.table("Round 1 initial parameters.txt", sep = "",
                                 header = TRUE))
paras.r6 <- as.matrix(read.table("Round 6 parameters.txt", sep = "", 
                                 header = TRUE))

# Prior-posterior overlap calculations.
dn.pri.post <- list(dn.prior = paras.r1[,1], dn.post = paras.r6[,1])
dn.overlap <- overlap(dn.pri.post)
# 13.24%

gamma.pri.post <- list(gamma.prior = paras.r1[,2], gamma.post = paras.r6[,2])
gamma.overlap <- overlap(gamma.pri.post)
# 72.49%

rn.pri.post <- list(rn.prior = paras.r1[,3], rn.post = paras.r6[,3])
rn.overlap <- overlap(rn.pri.post)
# 76.56%

eta.pri.post <- list(eta.prior = paras.r1[,4], eta.post = paras.r6[,4])
eta.overlap <- overlap(eta.pri.post)
# 70.68%

dm.pri.post <- list(dm.prior = paras.r1[,5], dm.post = paras.r6[,5])
dm.overlap <- overlap(dm.pri.post)
# 70.36%

alpha.pri.post <- list(alpha.prior = paras.r1[,6], r.init.post = paras.r6[,6])
alpha.overlap <- overlap(alpha.pri.post)
# 82.76%

p.ext.pri.post <- list(p.ext.prior = paras.r1[,7], p.ext.post = paras.r6[,7])
p.ext.overlap <- overlap(p.ext.pri.post)
# 67.58%

p.mit.pri.post <- list(p.mit.prior = paras.r1[,8], p.mit.post = paras.r6[,8])
p.mit.overlap <- overlap(p.mit.pri.post)
# 15.68%
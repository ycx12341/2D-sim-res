# Overlap.R
# Author: Yunchen Xiao
# This .R file calculates the prior-posterior overlap for the simulation study
# of T98G glioma patterns. 

# Clear the current workspace and load the necessary packages.
rm(list = ls())
library("overlapping")

# Read the prior and posterior parameters.
paras.r1 <- as.matrix(read.table("Round 1 initial parameters.txt", sep = "",
                                 header = TRUE))
paras.r7 <- as.matrix(read.table("Round 7 parameters.txt", sep = "", 
                                 header = TRUE))

# Prior-posterior overlap calculations.
dn.pri.post <- list(dn.prior = paras.r1[,1], dn.post = paras.r7[,1])
dn.overlap <- overlap(dn.pri.post)
# 0.58%

rn.pri.post <- list(rn.prior = paras.r1[,2], rn.post = paras.r7[,2])
rn.overlap <- overlap(rn.pri.post)
# 39.79%

r.init.pri.post <- list(r.init.prior = paras.r1[,3], r.init.post = paras.r7[,3])
r.init.overlap <- overlap(r.init.pri.post)
# 23.51%

p.ext.pri.post <- list(p.ext.prior = paras.r1[,4], p.ext.post = paras.r7[,4])
p.ext.overlap <- overlap(p.ext.pri.post)
# 61.91%

p.mit.pri.post <- list(p.mit.prior = paras.r1[,5], p.mit.post = paras.r7[,5])
p.mit.overlap <- overlap(p.mit.pri.post)
# 28.98%
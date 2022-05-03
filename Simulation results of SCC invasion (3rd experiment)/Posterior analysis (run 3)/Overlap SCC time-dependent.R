# Overlap SCC time-dependent.R
# Author: Yunchen Xiao
# This .R file calculates the prior-posterior overlap for the third experiment
# of the simulation study regarding the SCC invasion patterns. 

# Clear the current workspace and load the necessary packages.
rm(list = ls())
library("overlapping")

# Read the prior and posterior parameters.
paras.r1 <- as.matrix(read.table("Round 1 initial time varying parameters.txt", sep = "",
                                 header = TRUE))
paras.r9 <- as.matrix(read.table("Round 9 parameters.txt", sep = "", 
                                 header = TRUE))

# Prior-posterior overlap calculations.
dn.pri.post <- list(dn.prior = paras.r1[,1], dn.post = paras.r9[,1])
dn.overlap <- overlap(dn.pri.post)
# 1.89%

gamma.pri.post <- list(gamma.prior = paras.r1[,2], gamma.post = paras.r9[,2])
gamma.overlap <- overlap(gamma.pri.post)
# 1.89%

rn.pri.post <- list(rn.prior = paras.r1[,3], rn.post = paras.r9[,3])
rn.overlap <- overlap(rn.pri.post)
# 8.83%

eta.pri.post <- list(eta.prior = paras.r1[,4], eta.post = paras.r9[,4])
eta.overlap <- overlap(eta.pri.post)
# 10.31%

dm.pri.post <- list(dm.prior = paras.r1[,5], dm.post = paras.r9[,5])
dm.overlap <- overlap(dm.pri.post)
# 8.65%

alpha.pri.post <- list(alpha.prior = paras.r1[,6], r.init.post = paras.r9[,6])
alpha.overlap <- overlap(alpha.pri.post)
# 12.13%

r.init.pri.post <- list(r.init.prior = paras.r1[,7], r.init.post = paras.r9[,7])
r.init.overlap <- overlap(r.init.pri.post)
# 4.89%

p.ext.pri.post <- list(p.ext.prior = paras.r1[,8], p.ext.post = paras.r9[,8])
p.ext.overlap <- overlap(p.ext.pri.post)
# 7.44%

p.mit.pri.post <- list(p.mit.prior = paras.r1[,9], p.mit.post = paras.r9[,9])
p.mit.overlap <- overlap(p.mit.pri.post)
# 7.44%

dn.quad.pri.post <- list(dn.quad.prior = paras.r1[,10], dn.quad.post = paras.r9[,10])
dn.quad.overlap <- overlap(dn.quad.pri.post)
# 7.50%

dn.lin.pri.post <- list(dn.lin.prior = paras.r1[,11], dn.lin.post = paras.r9[,11])
dn.lin.overlap <- overlap(dn.lin.pri.post)
# 8.53%

gamma.quad.pri.post <- list(gamma.quad.prior = paras.r1[,12], gamma.quad.post = paras.r9[,12])
gamma.quad.overlap <- overlap(gamma.quad.pri.post)
# 10.16%

gamma.lin.pri.post <- list(gamma.lin.prior = paras.r1[,13], gamma.lin.post = paras.r9[,13])
gamma.lin.overlap <- overlap(gamma.lin.pri.post)
# 10.36%

rn.quad.pri.post <- list(rn.quad.prior = paras.r1[,14], rn.quad.post = paras.r9[,14])
rn.quad.overlap <- overlap(rn.quad.pri.post)
# 8.42%

rn.lin.pri.post <- list(rn.lin.prior = paras.r1[,15], rn.lin.post = paras.r9[,15])
rn.lin.overlap <- overlap(rn.lin.pri.post)
# 10.03%

eta.quad.pri.post <- list(eta.quad.prior = paras.r1[,16], eta.quad.post = paras.r9[,16])
eta.quad.overlap <- overlap(eta.quad.pri.post)
# 7.48%

eta.lin.pri.post <- list(eta.lin.prior = paras.r1[,17], eta.lin.post = paras.r9[,17])
eta.lin.overlap <- overlap(eta.lin.pri.post)
# 9.47%

alpha.quad.pri.post <- list(alpha.quad.prior = paras.r1[,18], alpha.quad.post = paras.r9[,18])
alpha.quad.overlap <- overlap(alpha.quad.pri.post)
# 6.74%

alpha.lin.pri.post <- list(alpha.lin.prior = paras.r1[,19], alpha.lin.post = paras.r9[,19])
alpha.lin.overlap <- overlap(alpha.lin.pri.post)
# 8.01%

p.mit.quad.pri.post <- list(p.mit.quad.prior = paras.r1[,20], p.mit.quad.post = paras.r9[,20])
p.mit.quad.overlap <- overlap(p.mit.quad.pri.post)
# 7.56%

p.mit.lin.pri.post <- list(p.mit.lin.prior = paras.r1[,21], p.mit.lin.post = paras.r9[,21])
p.mit.lin.overlap <- overlap(p.mit.lin.pri.post)
# 7.18%
# Overlap SCC d3.R
# Author: Yunchen Xiao
# This .R file calculates the prior-posterior overlap for the second experiment
# of the simulation study regarding the SCC invasion patterns. 
# Post day 3 pattern.

# Clear the current workspace and load the necessary packages.
rm(list = ls())
library("overlapping")

# Read the prior and posterior parameters.
paras.r1 <- as.matrix(read.table("Round 1 initial parameters.txt", sep = "",
                                 header = TRUE))
paras.r11 <- as.matrix(read.table("Round 11 parameters.txt", sep = "", 
                                 header = TRUE))

# Prior-posterior overlap calculations.
dn.pri.post <- list(dn.prior = paras.r1[,1], dn.post = paras.r11[,1])
dn.overlap <- overlap(dn.pri.post)
# 0.40%

gamma.pri.post <- list(gamma.prior = paras.r1[,2], gamma.post = paras.r11[,2])
gamma.overlap <- overlap(gamma.pri.post)
# 0.18%

rn.pri.post <- list(rn.prior = paras.r1[,3], rn.post = paras.r11[,3])
rn.overlap <- overlap(rn.pri.post)
# 11.83%

eta.pri.post <- list(eta.prior = paras.r1[,4], eta.post = paras.r11[,4])
eta.overlap <- overlap(eta.pri.post)
# 5.29%

dm.pri.post <- list(dm.prior = paras.r1[,5], dm.post = paras.r11[,5])
dm.overlap <- overlap(dm.pri.post)
# 11.09%

alpha.pri.post <- list(alpha.prior = paras.r1[,6], r.init.post = paras.r11[,6])
alpha.overlap <- overlap(alpha.pri.post)
# 10.27%

r.init.pri.post <- list(r.init.prior = paras.r1[,7], r.init.post = paras.r11[,7])
r.init.overlap <- overlap(r.init.pri.post)
# 5.97%

p.ext.pri.post <- list(p.ext.prior = paras.r1[,8], p.ext.post = paras.r11[,8])
p.ext.overlap <- overlap(p.ext.pri.post)
# 10.42%

p.mit.pri.post <- list(p.mit.prior = paras.r1[,9], p.mit.post = paras.r11[,9])
p.mit.overlap <- overlap(p.mit.pri.post)
# 8.03%
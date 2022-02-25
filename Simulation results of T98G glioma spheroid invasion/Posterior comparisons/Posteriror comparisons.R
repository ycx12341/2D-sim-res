# Posterior comparisons.R
# Author: Yunchen Xiao
# This .R file compares the prior and posterior samples obtained using non-error
# calibrated and error-calibrated ABC scheme.

rm(list = ls())

# Prior parameters
paras.r1 <- as.matrix(read.table("Round 1 initial parameters.txt", sep = "",
                                 header = TRUE))


# Posteriror parameter samples
paras.final.ec <- as.matrix(read.table("Round 6 parameters log transform.txt",
                                       sep = "", header = TRUE))
paras.final.non.ec <- as.matrix(read.table("Round 7 parameters.txt",
                                           sep = "", header = TRUE))
# dn comparisons
dn.post.ec <- paras.final.ec[,1]
dn.post.non.ec <- paras.final.non.ec[,1]

plot(density(dn.post.ec, adjust = 2.5), 
     xlim = c(0.000069, 0.00025), lwd = 3, main = "Posterior comparison of dn")
lines(density(dn.post.non.ec, adjust = 2.5), col = "red", lwd = 3)

# rn comparisons
rn.post.ec <- paras.final.ec[,2]
rn.post.non.ec <- paras.final.non.ec[,2]
plot(density(rn.post.ec, adjust = 5), 
     xlim = c(0, 0.0035), ylim = c(0, 720), lwd = 3, 
     main = "Posterior comparison of rn")
lines(density(rn.post.non.ec, adjust = 5), col = "red", lwd = 3)

# R_init. comparisons
r.init.post.ec <- paras.final.ec[,3]
r.init.post.non.ec <- paras.final.non.ec[,3]
plot(density(r.init.post.ec, adjust = 6), 
     lwd = 3, ylim = c(0,45), xlim = c(0.04, 0.15),
     main = "Posterior comparison of R_{init.}")
lines(density(r.init.post.non.ec, adjust = 6), col = "red", lwd = 3)

# P_ext. comparisons
p.ext.post.ec <- paras.final.ec[,4]
p.ext.post.non.ec <- paras.final.non.ec[,4]
plot(density(p.ext.post.ec, adjust = 5),
     lwd = 3, ylim = c(0, 17), xlim = c(0.01, 0.1),
     main = "Posterior comparison of P_{ext.}")
lines(density(p.ext.post.non.ec, adjust = 5), col = "red", lwd = 3)

# P_mit. comparisons
p.mit.post.ec <- paras.final.ec[,5]
p.mit.post.non.ec <- paras.final.non.ec[,5]
plot(density(p.mit.post.ec, adjust = 6),
     lwd = 3, xlim = c(0.2, 1), ylim = c(0,4.5),
     main = "Posterior comparison of P_{mit.}")
lines(density(p.mit.post.non.ec, adjust = 6), col = "red", lwd = 3)

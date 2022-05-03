# Full parameter estimates r8 run1.R
# Author: Yunchen Xiao
# This .R file reads in the parameter estimates obtained at the end of every
# round. It then generates the corresponding values at later periods of invasion
# using the averaged estimated regression coefficients. Finally, the values are
# combined into table and write into .txt files.

# Clear the environment and read the source files. 
rm(list = ls())
source("PDE 2D ABC functions adjusted SCC (time-dependent parameters).R")

# Read in the corresponding parameter samples obtained at the end of every
# round. 
paras.r1.run3 <- as.matrix(read.table("Round 1 initial time varying parameters.txt", 
                                      sep = "", header = TRUE))
paras.r2.run3 <- as.matrix(read.table("Round 2 parameters.txt", sep = "", 
                                      header = TRUE))
paras.r3.run3 <- as.matrix(read.table("Round 3 parameters.txt", sep = "", 
                                      header = TRUE))
paras.r4.run3 <- as.matrix(read.table("Round 4 parameters.txt", sep = "", 
                                      header = TRUE))
paras.r5.run3 <- as.matrix(read.table("Round 5 parameters.txt", sep = "", 
                                      header = TRUE))
paras.r6.run3 <- as.matrix(read.table("Round 6 parameters.txt", sep = "", 
                                      header = TRUE))
paras.r7.run3 <- as.matrix(read.table("Round 7 parameters.txt", sep = "", 
                                      header = TRUE))
paras.r8.run3 <- as.matrix(read.table("Round 8 parameters.txt", sep = "", 
                                      header = TRUE))

# Calculate the samples means for the parameters. 
paras.r1.run3.mean <- apply(paras.r1.run3, 2 ,mean)
paras.r2.run3.mean <- apply(paras.r2.run3, 2 ,mean)
paras.r3.run3.mean <- apply(paras.r3.run3, 2 ,mean)
paras.r4.run3.mean <- apply(paras.r4.run3, 2 ,mean)
paras.r5.run3.mean <- apply(paras.r5.run3, 2 ,mean)
paras.r6.run3.mean <- apply(paras.r6.run3, 2 ,mean)
paras.r7.run3.mean <- apply(paras.r7.run3, 2 ,mean)
paras.r8.run3.mean <- apply(paras.r8.run3, 2 ,mean)

# Combine the averaged parameter estimates for each round into a single 
# matrix. 
paras.mean.mat <- rbind(paras.r1.run3.mean, paras.r2.run3.mean, 
                        paras.r3.run3.mean, paras.r4.run3.mean, 
                        paras.r5.run3.mean, paras.r6.run3.mean,
                        paras.r7.run3.mean, paras.r8.run3.mean)

# For the time-dependent parameters, calculate their values at later periods 
# based on the averaged estimate of the 1st period and the averaged estimates
# of regression coefficients in each round. 
dn.varying.paras.vec <- dn.varying.paras(dn.init.par = paras.mean.mat[,1],
                                         dn.quad.temp = paras.mean.mat[,10],
                                         dn.lin.temp = paras.mean.mat[,11])

gamma.varying.paras.vec <- gamma.varying.paras(gamma.init.par = paras.mean.mat[,2],
                                           gamma.quad.temp = paras.mean.mat[,12],
                                           gamma.lin.temp = paras.mean.mat[,13])

rn.varying.paras.vec <- rn.varying.paras(rn.init.par = paras.mean.mat[,3],
                                         rn.quad.temp = paras.mean.mat[,14],
                                         rn.lin.temp = paras.mean.mat[,15])

eta.varying.paras.vec <- eta.varying.paras(eta.init.par = paras.mean.mat[,4],
                                           eta.quad.temp = paras.mean.mat[,16],
                                           eta.lin.temp = paras.mean.mat[,17])

alpha.varying.paras.vec <- alpha.varying.paras(alpha.init.par = paras.mean.mat[,6],
                                               alpha.quad.temp = paras.mean.mat[,18],
                                               alpha.lin.temp = paras.mean.mat[,19])

prob.prof.varying.paras.vec <- prob.prof.varying.paras(prob.prof.init.par = paras.mean.mat[,9],
                                                       prob.prof.quad.temp = paras.mean.mat[,20],
                                                       prob.prof.lin.temp = paras.mean.mat[,21])

# Combine these parameter values into a single matrix. 
paras.mean.mat.full <- cbind(dn.varying.paras.vec, gamma.varying.paras.vec,
                             rn.varying.paras.vec, eta.varying.paras.vec,
                             paras.mean.mat[,5], 
                             alpha.varying.paras.vec,
                             paras.mean.mat[,7], 
                             paras.mean.mat[,8],
                             prob.prof.varying.paras.vec)

# Extract the parameter estimates of the final round (final parameter estimates)
dn.vals.r8 <- paras.mean.mat.full[8, 1:7]
gamma.vals.r8 <- paras.mean.mat.full[8, 8:14]
rn.vals.r8 <- paras.mean.mat.full[8, 15:21]
eta.vals.r8 <- paras.mean.mat.full[8, 22:28]
dm.vals.r8 <- paras.mean.mat.full[8, 29]
alpha.vals.r8 <- paras.mean.mat.full[8, 30:36]
r.init.vals.r8 <- paras.mean.mat.full[8, 37]
p.ext.vals.r8 <- paras.mean.mat.full[8, 38]
p.mit.vals.r8 <- paras.mean.mat.full[8, 39:45]

paras.r8.mean.mat.full <- rbind(dn.vals.r8, gamma.vals.r8, rn.vals.r8,
                                eta.vals.r8, c(rep(dm.vals.r8, 5),rep(NA, 2)), 
                                alpha.vals.r8, c(rep(r.init.vals.r8, 5),rep(NA, 2)),
                                c(rep(p.ext.vals.r8, 5), rep(NA, 2)), p.mit.vals.r8)

rownames(paras.r8.mean.mat.full) <- c("dn.vals.r8", "gamma.vals.r8", "rn.vals.r8",
                                      "eta.vals.r8", "dm.vals.r8", "alpha.vals.r8",
                                      "r.init.vals.r8", "p.ext.vals.r8", "p.mit.vals.r8")

colnames(paras.r8.mean.mat.full) <- c("init.par", "p2", "p3", "p4", 
                                      "p5", "quad.coef", "lin.coef")

# Write the matrix that contains the full parameter estimates and the one that
# contains the parameter estimates of the final round into .txt files.

write.table(paras.mean.mat.full, "Full parameter estimates 8 rounds run 1.txt")
write.table(paras.r8.mean.mat.full, "Round 8 final full parameter estimates run 1.txt")

# Save the regression models fitted to the final parameter estimates of 
# time-dependent parameters in .rds files

# Explanatory variables and the x-values to be predicted by the regression 
# model.
x <- c(1,2,3,4,5)
x2 <- x^2
periods.cont <- seq(0, 6, by = 0.01)

# dn
dn.quad.mod <- lm(dn.vals.r8[1:5] ~ x + x2)
fitted.vals.predict.dn.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.dn.quad[i] <- predict(dn.quad.mod,
                                 newdata = data.frame(x = periods.cont[i],
                                 x2 = periods.cont[i]^2))
}

# gamma 
gamma.quad.mod <- lm(gamma.vals.r8[1:5] ~ x + x2)
fitted.vals.predict.gamma.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.gamma.quad[i] <- predict(gamma.quad.mod,
                                    newdata = data.frame(x = periods.cont[i],
                                    x2 = periods.cont[i]^2))
}

# rn 
rn.quad.mod <- lm(rn.vals.r8[1:5] ~ x + x2)
fitted.vals.predict.rn.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.rn.quad[i] <- predict(rn.quad.mod,
                                    newdata = data.frame(x = periods.cont[i],
                                    x2 = periods.cont[i]^2))
}

# eta 
eta.quad.mod <- lm(eta.vals.r8[1:5] ~ x + x2)
fitted.vals.predict.eta.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.eta.quad[i] <- predict(eta.quad.mod,
                                     newdata = data.frame(x = periods.cont[i],
                                     x2 = periods.cont[i]^2))
}

# alpha
alpha.quad.mod <- lm(alpha.vals.r8[1:5] ~ x + x2)
fitted.vals.predict.alpha.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.alpha.quad[i] <- predict(alpha.quad.mod,
                                     newdata = data.frame(x = periods.cont[i],
                                     x2 = periods.cont[i]^2))
}

# p.mit 
p.mit.quad.mod <- lm(p.mit.vals.r8[1:5] ~ x + x2)
fitted.vals.predict.p.mit.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.p.mit.quad[i] <- predict(p.mit.quad.mod,
                                       newdata = data.frame(x = periods.cont[i],
                                       x2 = periods.cont[i]^2))
}

# Combine these fitted values into a matrix and write it into a .txt file.
reg.est.vals <- cbind(fitted.vals.predict.dn.quad, 
                      fitted.vals.predict.gamma.quad, 
                      fitted.vals.predict.rn.quad,
                      fitted.vals.predict.eta.quad, 
                      rep(dm.vals.r8, length(periods.cont)),
                      fitted.vals.predict.alpha.quad,
                      rep(r.init.vals.r8, length(periods.cont)),
                      rep(p.ext.vals.r8, length(periods.cont)),
                      fitted.vals.predict.p.mit.quad)

write.table(reg.est.vals, "Estimated final parameter values regression.txt")

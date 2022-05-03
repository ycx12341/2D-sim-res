# Initial time varying paras.R
# Author: Yunchen Xiao

# This .R file reads in the regression models (mostly quadratic) fitted to the
# parameter estimates for each SCC invasion pattern. Generates time-varying
# parameters at different periods of invasion which fall in the bounds of the
# prior distribution

# Clear the workspace and read in the source functions. 
rm(list = ls())
library(doParallel)
library(doRNG)
library(tictoc)
library(readr)
source("PDE 2D ABC functions adjusted SCC (time-dependent parameters).R")

# Regression coefficients
dn.regression.model <- read_rds("dn quadratic regression model.rds")
dn.lin.coef <- unname(dn.regression.model$coefficients[2])
dn.quad.coef <- unname(dn.regression.model$coefficients[3])

gamma.regression.model <- read_rds("gamma quadratic regression model.rds")
gamma.lin.coef <- unname(gamma.regression.model$coefficients[2])
gamma.quad.coef <- unname(gamma.regression.model$coefficients[3])

rn.regression.model <- read_rds("rn quadratic regression model.rds")
rn.lin.coef <- unname(rn.regression.model$coefficients[2])
rn.quad.coef <- unname(rn.regression.model$coefficients[3])

eta.regression.model <- read_rds("eta quadratic regression model.rds")
eta.lin.coef <- unname(eta.regression.model$coefficients[2])
eta.quad.coef <- unname(eta.regression.model$coefficients[3])

alpha.regression.model <- read_rds("alpha quadratic regression model.rds")
alpha.lin.coef <- unname(alpha.regression.model$coefficients[2])
alpha.quad.coef <- unname(alpha.regression.model$coefficients[3])

prob.prof.regression.model <- read_rds("prob.prof quadratic regression model.rds")
prob.prof.lin.coef <- unname(prob.prof.regression.model$coefficients[2])
prob.prof.quad.coef <- unname(prob.prof.regression.model$coefficients[3])

# dn parameters

# Set the counter to be 0 initially. 
valid.counter.dn <- 0

# An empty matrix used to store the dn parameters and their corresponding 
# regression coefficients that can generate valid time-varying parameters.
dn.paras.table <- vector()

set.seed(874514)
RNGkind(sample.kind = "Rejection")
# Objective: 50000 sets of initial dn parameter values and their corresponding
# regression coefficients that can generate valid time-varying parameters. 
while(valid.counter.dn < 50000) {
  # Sample initial dn values from its prior. 
  dn.temp <- runif(10000, 0.000069, 0.02)
  # Sample regression coefficients from uniform distributions which use bounds
  # associated to the corresponding regression coefficients estimates. 
  slope.dn.quad.temp <- runif(10000, dn.quad.coef/2, dn.quad.coef*2)
  slope.dn.linear.temp <- runif(10000, dn.lin.coef*2, dn.lin.coef/2)
  # Check how many sets out of 10000 can generate parameter values at 5 
  # different time points which are all supported by the prior distribution. 
  dn.paras.temp <- dn.varying.paras(dn.init.par = dn.temp, 
                                      dn.lin.temp = slope.dn.linear.temp,
                                      dn.quad.temp = slope.dn.quad.temp)
  # Update the parameter matrix and the counter. 
  dn.paras.table <- rbind(dn.paras.table, dn.paras.temp)
  valid.counter.dn <- valid.counter.dn + length(dn.paras.temp[,1])
}  

# Check if these parameters and their associated regression coefficients
# are valid again. 
dn.paras.check <- dn.varying.paras(dn.init.par = dn.paras.table[,1],
                                   dn.quad.temp = dn.paras.table[,6],
                                   dn.lin.temp = dn.paras.table[,7])

# Truncate the size of the parameter matrix to 50000 rows if the parameter 
# matrix generated in the end contains more than 50000 rows.  
dn.paras <- dn.paras.check[1:50000, c(1,6,7)]

# gamma parameters
valid.counter.gamma <- 0
gamma.paras.table <- vector()
set.seed(874514)
RNGkind(sample.kind = "Rejection")
while(valid.counter.gamma < 50000) {
  gamma.temp <- runif(10000, 0.005, 0.26)
  slope.gamma.quad.temp <- runif(10000, gamma.quad.coef*2, gamma.quad.coef/2)
  slope.gamma.linear.temp <- runif(10000, gamma.lin.coef/2, gamma.lin.coef*2)
  gamma.paras.temp <- gamma.varying.paras(gamma.init.par = gamma.temp, 
                                    gamma.lin.temp = slope.gamma.linear.temp,
                                    gamma.quad.temp = slope.gamma.quad.temp)
  gamma.paras.table <- rbind(gamma.paras.table, gamma.paras.temp)
  valid.counter.gamma <- valid.counter.gamma + length(gamma.paras.temp[,1])
}  

gamma.paras.check <- gamma.varying.paras(gamma.init.par = gamma.paras.table[,1],
                                   gamma.quad.temp = gamma.paras.table[,6],
                                   gamma.lin.temp = gamma.paras.table[,7])

gamma.paras <- gamma.paras.check[1:50000, c(1,6,7)]

# rn parameters
valid.counter.rn <- 0
rn.paras.table <- vector()
set.seed(874514)
RNGkind(sample.kind = "Rejection")
while(valid.counter.rn < 50000) {
  rn.temp <- runif(10000, 0.0008, 0.08)
  slope.rn.quad.temp <- runif(10000, rn.quad.coef*2, rn.quad.coef/2)
  slope.rn.linear.temp <- runif(10000, rn.lin.coef/2, rn.lin.coef*2)
  rn.paras.temp <- rn.varying.paras(rn.init.par = rn.temp, 
                                          rn.lin.temp = slope.rn.linear.temp,
                                          rn.quad.temp = slope.rn.quad.temp)
  rn.paras.table <- rbind(rn.paras.table, rn.paras.temp)
  valid.counter.rn <- valid.counter.rn + length(rn.paras.temp[,1])
}  

rn.paras.check <- rn.varying.paras(rn.init.par = rn.paras.table[,1],
                                         rn.quad.temp = rn.paras.table[,6],
                                         rn.lin.temp = rn.paras.table[,7])

rn.paras <- rn.paras.check[1:50000, c(1,6,7)]

# eta parameters 
valid.counter.eta <- 0
eta.paras.table <- vector()
set.seed(874514)
RNGkind(sample.kind = "Rejection")
while(valid.counter.eta < 50000) {
  eta.temp <- runif(10000, 7, 18)
  slope.eta.quad.temp <- runif(10000, eta.quad.coef/2, eta.quad.coef*2)
  slope.eta.linear.temp <- runif(10000, eta.lin.coef*2, eta.lin.coef/2)
  eta.paras.temp <- eta.varying.paras(eta.init.par = eta.temp, 
                                    eta.lin.temp = slope.eta.linear.temp,
                                    eta.quad.temp = slope.eta.quad.temp)
  eta.paras.table <- rbind(eta.paras.table, eta.paras.temp)
  valid.counter.eta <- valid.counter.eta + length(eta.paras.temp[,1])
}  
eta.paras.check <- eta.varying.paras(eta.init.par = eta.paras.table[,1],
                                   eta.quad.temp = eta.paras.table[,6],
                                   eta.lin.temp = eta.paras.table[,7])
eta.paras <- eta.paras.check[1:50000, c(1,6,7)]

# alpha parameters
valid.counter.alpha <- 0
alpha.paras.table <- vector()
set.seed(874514)
RNGkind(sample.kind = "Rejection")
while(valid.counter.alpha < 50000) {
  alpha.temp <- runif(10000, 0.07, 0.18)
  slope.alpha.quad.temp <- runif(10000, alpha.quad.coef/2, alpha.quad.coef*2)
  slope.alpha.linear.temp <- runif(10000, alpha.lin.coef*2, alpha.lin.coef/2)
  alpha.paras.temp <- alpha.varying.paras(alpha.init.par = alpha.temp, 
                                      alpha.lin.temp = slope.alpha.linear.temp,
                                      alpha.quad.temp = slope.alpha.quad.temp)
  alpha.paras.table <- rbind(alpha.paras.table, alpha.paras.temp)
  valid.counter.alpha <- valid.counter.alpha + length(alpha.paras.temp[,1])
}  
alpha.paras.check <- alpha.varying.paras(alpha.init.par = alpha.paras.table[,1],
                                     alpha.quad.temp = alpha.paras.table[,6],
                                     alpha.lin.temp = alpha.paras.table[,7])
alpha.paras <- alpha.paras.check[1:50000, c(1,6,7)]

# prob.prof parameters 
valid.counter.prob.prof <- 0
prob.prof.paras.table <- vector()

set.seed(874514)
RNGkind(sample.kind = "Rejection")
while(valid.counter.prob.prof < 50000) {
  prob.prof.temp <- runif(10000, 0.2, 1)
  slope.prob.prof.quad.temp <- runif(10000, prob.prof.quad.coef*2, 
                                     prob.prof.quad.coef/2)
  slope.prob.prof.linear.temp <- runif(10000, prob.prof.lin.coef/2, 
                                       prob.prof.lin.coef*2)
  prob.prof.paras.temp <- prob.prof.varying.paras(prob.prof.init.par = prob.prof.temp, 
                                                  prob.prof.lin.temp = slope.prob.prof.linear.temp,
                                                  prob.prof.quad.temp = slope.prob.prof.quad.temp)
  prob.prof.paras.table <- rbind(prob.prof.paras.table, prob.prof.paras.temp)
  valid.counter.prob.prof <- valid.counter.prob.prof + 
    length(prob.prof.paras.temp[,1])
}  

prob.prof.paras.check <- prob.prof.varying.paras(prob.prof.init.par = prob.prof.paras.table[,1],
                                         prob.prof.quad.temp = prob.prof.paras.table[,6],
                                         prob.prof.lin.temp = prob.prof.paras.table[,7])

prob.prof.paras <- prob.prof.paras.check[1:50000, c(1,6,7)]

# Parameters that are not time-dependent 
set.seed(874514)
RNGkind(sample.kind = "Rejection")
dm.paras <- runif(50000, 0.0001, 0.033)
init.cells.cols.paras <- runif(50000, 1, 5)
prob.death.paras <- runif(50000, 0.01, 0.1)

# Combine the parameters sampled into a single matrix. 
paras.table <- cbind(dn.paras[,1], gamma.paras[,1], rn.paras[,1], eta.paras[,1], 
                     dm.paras, alpha.paras[,1], init.cells.cols.paras, # 1 - 7
                     prob.death.paras, prob.prof.paras[,1], # 8 - 9
                     dn.paras[,2], dn.paras[,3], # 10 - 11 
                     gamma.paras[,2], gamma.paras[,3], # 12 - 13
                     rn.paras[,2], rn.paras[,3], # 14 - 15
                     eta.paras[,2], eta.paras[,3], # 16 - 17
                     alpha.paras[,2], alpha.paras[,3], # 18 - 19
                     prob.prof.paras[,2], prob.prof.paras[,3]) # 20 - 21))

colnames(paras.table) <- c("dn", "gamma", "rn", "eta", "dm", "alpha",
                        "init.cells.cols", "prob.death", "prob.prof",
                        "dn.quad", "dn.lin", "gamma.quad", "gamma.lin",
                        "rn.quad", "rn.lin", "eta.quad", "eta.lin",
                        "alpha.quad", "alpha.lin", "prob.prof.quad",
                        "prob.prof.lin")

# Store the parameters in a single .txt file. 
write.table(paras.table, "Round 1 initial time varying parameters.txt")



# Posterior results summary.R
# Author: Yunchen Xiao
# This .R file reads in the final simulation output from the 3 different runs 
# of the time-varying SCC simulation pattern and calculate their Bhattacharyya 
# distance (multivariate version) with the reference data. In addition, the 
# plots of actual ESS and bandwidth factors at the end of different rounds are 
# generated in this file. 

# Clear the workspace, load the necessary package and the file which contains
# the source functions. 
rm(list = ls())
library(readr)
library(fpc)
source("PDE 2D ABC functions adjusted SCC (time-dependent parameters).R")

# Read in the final simulation output from three different runs. 
paras.r8.run1.res <- read_rds("paras r8 average result run 1.rds")
paras.r10.run2.res <- read_rds("paras r10 average result run 2.rds")
paras.r9.run3.res <- read_rds("paras r9 average result run 3.rds")

# Read in the final parameter estimates from the three different runs,
# combine them into a matrix by row. 
paras.r8.run1 <- as.matrix(read.table("Round 8 parameters run 1.txt", 
                                      sep = "", header = TRUE))
paras.r10.run2 <- as.matrix(read.table("Round 10 parameters run 2.txt", 
                                      sep = "", header = TRUE))
paras.r9.run3 <- as.matrix(read.table("Round 9 parameters run 3.txt", 
                                      sep = "", header = TRUE))

paras.r8.run1.mean <- apply(paras.r8.run1, 2, mean)
paras.r10.run2.mean <- apply(paras.r10.run2, 2, mean)
paras.r9.run3.mean <- apply(paras.r9.run3, 2, mean)
paras.ests <- rbind(paras.r8.run1.mean, paras.r10.run2.mean, 
                    paras.r9.run3.mean)

# B-C distance calculations. 
run1.den.mat.d3 <- paras.r8.run1.res$den.mat.d3
run2.den.mat.d3 <- paras.r10.run2.res$den.mat.d3
run3.den.mat.d3 <- paras.r9.run3.res$den.mat.d3

run1.den.mat.d6 <- paras.r8.run1.res$den.mat.d6
run2.den.mat.d6 <- paras.r10.run2.res$den.mat.d6
run3.den.mat.d6 <- paras.r9.run3.res$den.mat.d6

run1.den.mat.d9 <- paras.r8.run1.res$den.mat.d9
run2.den.mat.d9 <- paras.r10.run2.res$den.mat.d9
run3.den.mat.d9 <- paras.r9.run3.res$den.mat.d9

run1.den.mat.d12 <- paras.r8.run1.res$den.mat.d12
run2.den.mat.d12 <- paras.r10.run2.res$den.mat.d12
run3.den.mat.d12 <- paras.r9.run3.res$den.mat.d12

run1.den.mat.d14 <- paras.r8.run1.res$den.mat.d14
run2.den.mat.d14 <- paras.r10.run2.res$den.mat.d14
run3.den.mat.d14 <- paras.r9.run3.res$den.mat.d14

# Vector of reference means
ref.mean <- c(mean(t3.ref.den), mean(t6.ref.den), mean(t9.ref.den),
              mean(t12.ref.den), mean(t14.ref.den))
# Covariance matrix
ref.data <- data.frame(t3.ref.den = as.numeric(t3.ref.den),
                       t6.ref.den = as.numeric(t6.ref.den),
                       t9.ref.den = as.numeric(t9.ref.den),
                       t12.ref.den = as.numeric(t12.ref.den),
                       t14.ref.den = as.numeric(t14.ref.den))
ref.cov <- cov(ref.data)

### Run 1 ###
run1.mean <- c(mean(run1.den.mat.d3), mean(run1.den.mat.d6),
               mean(run1.den.mat.d9), mean(run1.den.mat.d12),
               mean(run1.den.mat.d14))

run1.data <- data.frame(run1.den.mat.d3 = as.numeric(run1.den.mat.d3),
                        run1.den.mat.d6 = as.numeric(run1.den.mat.d6),
                        run1.den.mat.d9 = as.numeric(run1.den.mat.d9),
                        run1.den.mat.d12 = as.numeric(run1.den.mat.d12),
                        run1.den.mat.d14 = as.numeric(run1.den.mat.d14)
                        )
run1.cov <- cov(run1.data)

bcd.run1.ref <- bhattacharyya.dist(mu1 = ref.mean, mu2 = run1.mean, 
                                   Sigma1 = ref.cov, Sigma2 = run1.cov)
# 0.576

### Run 2 ###
run2.mean <- c(mean(run2.den.mat.d3), mean(run2.den.mat.d6),
               mean(run2.den.mat.d9), mean(run2.den.mat.d12),
               mean(run2.den.mat.d14))

run2.data <- data.frame(run2.den.mat.d3 = as.numeric(run2.den.mat.d3),
                        run2.den.mat.d6 = as.numeric(run2.den.mat.d6),
                        run2.den.mat.d9 = as.numeric(run2.den.mat.d9),
                        run2.den.mat.d12 = as.numeric(run2.den.mat.d12),
                        run2.den.mat.d14 = as.numeric(run2.den.mat.d14)
)
run2.cov <- cov(run2.data)

bcd.run2.ref <- bhattacharyya.dist(mu1 = ref.mean, mu2 = run2.mean, 
                                   Sigma1 = ref.cov, Sigma2 = run2.cov)
# 0.506

### Run 3 ###
run3.mean <- c(mean(run3.den.mat.d3), mean(run3.den.mat.d6),
               mean(run3.den.mat.d9), mean(run3.den.mat.d12),
               mean(run3.den.mat.d14))

run3.data <- data.frame(run3.den.mat.d3 = as.numeric(run3.den.mat.d3),
                        run3.den.mat.d6 = as.numeric(run3.den.mat.d6),
                        run3.den.mat.d9 = as.numeric(run3.den.mat.d9),
                        run3.den.mat.d12 = as.numeric(run3.den.mat.d12),
                        run3.den.mat.d14 = as.numeric(run3.den.mat.d14)
)
run3.cov <- cov(run3.data)

bcd.run3.ref <- bhattacharyya.dist(mu1 = ref.mean, mu2 = run3.mean, 
                                   Sigma1 = ref.cov, Sigma2 = run3.cov)
# 0.396

# Based on the bcd results, the final pattern produced at the end of run 3 should 
# be taken as the final result. 

write_rds(list(paras.ests.3runs = paras.ests, paras.final.est = paras.r9.run3.mean),
          "SCC time varying final parameter estimates non err calib.rds")

# Plots of actual ESS and corresponding bandwidth factors at the end of 
# each round. 
ess.bw.run1 <- as.matrix(read.table("Run 1 ESS BW.txt", sep = "", header = TRUE))
ess.bw.run2 <- as.matrix(read.table("Run 2 ESS BW.txt", sep = "", header = TRUE))
ess.bw.run3 <- as.matrix(read.table("Run 3 ESS BW.txt", sep = "", header = TRUE))


# Plot of actual ESS at the end of each round. (3 different runs)
plot(x = seq(1,7,by = 1), y = ess.bw.run1[,1], lwd = 2, xlab = "Rounds",
     ylab = "Effective Sample Size", xlim = c(1,9), ylim = c(4500,12600), 
     col = "purple", main = "Actual effective sample size at different rounds")
points(x = seq(1,9,by = 1), y = ess.bw.run2[,1], pch = 2, lwd = 2, col = "red")
points(x = seq(1,8,by = 1), y = ess.bw.run3[,1], pch = 10, lwd = 2, col = "blue")
abline(h = 7500, col = "black", lty = 2, lwd = 2)
legend(7.5, 2400, legend = c("Run 1", "Run 2", "Run 3"), 
       col = c("purple", "red", "blue", "black"), pch = c(1,2,10))

# Plot of corresponding bandwidth factors at the end of each round (3 
# different runs).
plot(x = seq(1,7,by = 1), y = ess.bw.run1[,2], lwd = 2, xlab = "Rounds",
     ylab = "Weight factor", xlim = c(1,9), ylim = c(0.15,12),
     col = "purple", main = "Weight factor at different rounds")
points(x = seq(1,9,by = 1), y = ess.bw.run2[,2], pch = 2, lwd = 2, col = "red")
points(x = seq(1,8,by = 1), y = ess.bw.run3[,2], pch = 10, lwd = 2, col = "blue")
legend(7.5, 10, legend = c("Run 1", "Run 2", "Run 3"), 
       col = c("purple", "red", "blue", "black"), pch = c(1,2,10))
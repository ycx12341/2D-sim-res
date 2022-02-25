# Posterior results summary.R
# Author: Yunchen Xiao
# This .R file reads in the final simulation output from the 3 different runs 
# of the T98G glioma simulation pattern and calculate their Bhattacharyya 
# distance (multivariate version) with the reference data. In addition, the 
# plots of actual ESS and bandwidth factors at the end of different rounds are 
# generated in this file. 

# Clear the workspace, load the necessary package and the file which contains
# the source functions. 
rm(list = ls())
library(readr)
library(fpc)
source("PDE 2D ABC functions adjusted.R")

# Read in the final simulation output from three different runs. 
paras.r7.run1.res <- read_rds("paras r7 average result run 1.rds")
paras.r7.run2.res <- read_rds("paras r7 average result run 2.rds")
paras.r7.run3.res <- read_rds("paras r7 average result run 3.rds")

# Read in the final parameter estimates from the three different runs,
# combine them into a matrix by row. 
paras.r7.run1 <- as.matrix(read.table("Round 7 parameters run 1.txt", 
                                      sep = "", header = TRUE))
paras.r7.run2 <- as.matrix(read.table("Round 7 parameters run 2.txt", 
                                      sep = "", header = TRUE))
paras.r7.run3 <- as.matrix(read.table("Round 7 parameters run 3.txt", 
                                      sep = "", header = TRUE))
paras.r7.run1.mean <- apply(paras.r7.run1, 2, mean)
paras.r7.run2.mean <- apply(paras.r7.run2, 2, mean)
paras.r7.run3.mean <- apply(paras.r7.run3, 2, mean)
paras.ests <- rbind(paras.r7.run1.mean, paras.r7.run2.mean, 
                    paras.r7.run3.mean)

# B-C distance calculations. 
run1.den.mat.d1 <- paras.r7.run1.res$den.mat.d1
run2.den.mat.d1 <- paras.r7.run2.res$den.mat.d1
run3.den.mat.d1 <- paras.r7.run3.res$den.mat.d1

run1.den.mat.d3 <- paras.r7.run1.res$den.mat.d3
run2.den.mat.d3 <- paras.r7.run2.res$den.mat.d3
run3.den.mat.d3 <- paras.r7.run3.res$den.mat.d3

# Vector of reference means
ref.mean <- c(mean(t1.ref.den), mean(t3.ref.den))
# Covariance matrix
ref.data <- data.frame(t1.ref.den = as.numeric(t1.ref.den),
                       t3.ref.den = as.numeric(t3.ref.den))
ref.cov <- cov(ref.data)

### Run 1 ###
run1.mean <- c(mean(run1.den.mat.d1), mean(run1.den.mat.d3))
run1.data <- data.frame(run1.den.mat.d1 = as.numeric(run1.den.mat.d1),
                        run1.den.mat.d3 = as.numeric(run1.den.mat.d3))
run1.cov <- cov(run1.data)

bcd.run1.ref <- bhattacharyya.dist(mu1 = ref.mean, mu2 = run1.mean, 
                                     Sigma1 = ref.cov, Sigma2 = run1.cov)
# 0.03984897

### Run 2 ###
run2.mean <- c(mean(run2.den.mat.d1), mean(run2.den.mat.d3))
run2.data <- data.frame(run2.den.mat.d1 = as.numeric(run2.den.mat.d1),
                        run2.den.mat.d3 = as.numeric(run2.den.mat.d3))
run2.cov <- cov(run2.data)
bcd.run2.ref <- bhattacharyya.dist(mu1 = ref.mean, mu2 = run2.mean,
                                     Sigma1 = ref.cov, Sigma2 = run2.cov)
# 0.0184112

### Run 3 ###
run3.mean <- c(mean(run3.den.mat.d1), mean(run3.den.mat.d3))
run3.data <- data.frame(run3.den.mat.d1 = as.numeric(run3.den.mat.d1),
                        run3.den.mat.d3 = as.numeric(run3.den.mat.d3))
run3.cov <- cov(run3.data)
bcd.run3.ref <- bhattacharyya.dist(mu1 = ref.mean, mu2 = run3.mean,
                                     Sigma1 = ref.cov, Sigma2 = run3.cov)
# 0.0371751

# Based on the bcd results, the final pattern produced at the end of run 2 should 
# be taken as the final result. 

write_rds(list(paras.ests.3runs = paras.ests, paras.final.est = paras.r7.run2.mean),
          "Glioma final parameter estimates non err calib.rds")

# Plots of actual ESS and corresponding bandwidth factors at the end of 
# each round. 
ess.bw.run1 <- as.matrix(read.table("Run 1 ESS BW.txt", sep = "", header = TRUE))
ess.bw.run2 <- as.matrix(read.table("Run 2 ESS BW.txt", sep = "", header = TRUE))
ess.bw.run3 <- as.matrix(read.table("Run 3 ESS BW.txt", sep = "", header = TRUE))


# Plot of actual ESS at the end of each round. (3 different runs)
plot(x = seq(1,6,by = 1), y = ess.bw.run1[,1], lwd = 2, xlab = "Rounds",
     ylab = "Effective Sample Size", xlim = c(1,9), ylim = c(1500,2600), 
     col = "purple", main = "Actual effective sample size at different rounds")
points(x = seq(1,6,by = 1), y = ess.bw.run2[,1], pch = 2, lwd = 2, col = "red")
points(x = seq(1,6,by = 1), y = ess.bw.run3[,1], pch = 10, lwd = 2, col = "blue")
abline(h = 1500, col = "black", lty = 2, lwd = 2)
legend(7.5, 2400, legend = c("Run 1", "Run 2", "Run 3"), 
       col = c("purple", "red", "blue", "black"), pch = c(1,2,10))

# Plot of corresponding bandwidth factors at the end of each round (3 
# different runs).
plot(x = seq(1,6,by = 1), y = ess.bw.run1[,2], lwd = 2, xlab = "Rounds",
     ylab = "Weight factor", xlim = c(1,9), ylim = c(1.5,10),
     col = "purple", main = "Weight factor at different rounds")
points(x = seq(1,6,by = 1), y = ess.bw.run2[,2], pch = 2, lwd = 2, col = "red")
points(x = seq(1,6,by = 1), y = ess.bw.run3[,2], pch = 10, lwd = 2, col = "blue")
legend(7.5, 10, legend = c("Run 1", "Run 2", "Run 3"), 
       col = c("purple", "red", "blue", "black"), pch = c(1,2,10))


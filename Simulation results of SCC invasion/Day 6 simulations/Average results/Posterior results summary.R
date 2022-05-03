# Posterior results summary.R
# Author: Yunchen Xiao
# This .R file reads in the final simulation output from the 3 different runs 
# of the post-day 6 SCC pattern and calculate their Bhattacharyya distance with 
# the reference data. The final output which yields the minimum B-C distance was
# chosen as the final result of the post-day 6 SCC pattern's simulation, and 
# the initial condition for the post-day 9 SCC pattern's simulation. In 
# addition, the plots of actual ESS and bandwidth factors at the end of 
# different rounds are generated in this file. 

# Clear the workspace, load the necessary package and the file which contains
# the source functions. 
rm(list = ls())
library(readr)
source("PDE 2D ABC functions adjusted SCC.R")

# Read in the final simulation output from three different runs. 
paras.r6.run1.res <- read_rds("paras r6 average result run 1.rds")
paras.r6.run2.res <- read_rds("paras r6 average result run 2.rds")
paras.r6.run3.res <- read_rds("paras r6 average result run 3.rds")

# Read in the final parameter estimates from the three different runs,
# combine them into a matrix by row. 
paras.r6.run1 <- as.matrix(read.table("Round 6 parameters run 1.txt", 
                                      sep = "", header = TRUE))
paras.r6.run2 <- as.matrix(read.table("Round 6 parameters run 2.txt", 
                                      sep = "", header = TRUE))
paras.r6.run3 <- as.matrix(read.table("Round 6 parameters run 3.txt", 
                                      sep = "", header = TRUE))
paras.r6.run1.mean <- apply(paras.r6.run1, 2, mean)
paras.r6.run2.mean <- apply(paras.r6.run2, 2, mean)
paras.r6.run3.mean <- apply(paras.r6.run3, 2, mean)
paras.day6 <- rbind(paras.r6.run1.mean, paras.r6.run2.mean, 
                    paras.r6.run3.mean)

# B-C distance calculations. 
run1.den.mat <- paras.r6.run1.res$den.mat.d6
run2.den.mat <- paras.r6.run2.res$den.mat.d6
run3.den.mat <- paras.r6.run3.res$den.mat.d6

bcd.run1.ref <- 0.25 * (log(0.25 * (var(c(run1.den.mat))/var(c(t6.ref.den)) + 
                                      var(c(t6.ref.den))/var(c(run1.den.mat)) + 2))) + 
  0.25 * (((mean(run1.den.mat) - mean(t6.ref.den))^2)/(var(c(t6.ref.den)) + var(c(run1.den.mat))))

# 2.242068e-05

bcd.run2.ref <- 0.25 * (log(0.25 * (var(c(run2.den.mat))/var(c(t6.ref.den)) + 
                                      var(c(t6.ref.den))/var(c(run2.den.mat)) + 2))) + 
  0.25 * (((mean(run2.den.mat) - mean(t6.ref.den))^2)/(var(c(t6.ref.den)) + var(c(run2.den.mat))))

# 1.332816e-05

bcd.run3.ref <- 0.25 * (log(0.25 * (var(c(run3.den.mat))/var(c(t6.ref.den)) + 
                                      var(c(t6.ref.den))/var(c(run3.den.mat)) + 2))) + 
  0.25 * (((mean(run3.den.mat) - mean(t6.ref.den))^2)/(var(c(t6.ref.den)) + var(c(run3.den.mat))))

# 0.0002641626

# Based on the bcd results, the final pattern produced at the end of run 2 should 
# be taken as the final result. 

write_rds(list(paras.day6.3runs = paras.day6, paras.day6.final.est = paras.r6.run2.mean),
          "Day 6 parameter estimates.rds")

# Plots of actual ESS and corresponding bandwidth factors at the end of 
# each round. 
ess.bw.run1 <- as.matrix(read.table("Run 1 ESS BW.txt", sep = "", header = TRUE))
ess.bw.run2 <- as.matrix(read.table("Run 2 ESS BW.txt", sep = "", header = TRUE))
ess.bw.run3 <- as.matrix(read.table("Run 3 ESS BW.txt", sep = "", header = TRUE))


# Plot of actual ESS at the end of each round. (3 different runs)
plot(x = seq(1,5,by = 1), y = ess.bw.run1[,1], lwd = 2, xlab = "Rounds",
     ylab = "Effective Sample Size", xlim = c(1,7), ylim = c(1400,2600), 
     col = "purple", 
     main = "Actual effective sample size at different rounds")
points(x = seq(1,5,by = 1), y = ess.bw.run2[,1], pch = 2, lwd = 2, 
       col = "red")
points(x = seq(1,5,by = 1), y = ess.bw.run3[,1], pch = 10, lwd = 2, 
       col = "blue")
abline(h = 1500, col = "black", lty = 2, lwd = 2)
legend(6, 2400, legend = c("Run 1", "Run 2", "Run 3"), 
       col = c("purple", "red", "blue", "black"), pch = c(1,2,10))

# Plot of corresponding bandwidth factors at the end of each round (3 
# different runs).
plot(x = seq(1,5,by = 1), y = ess.bw.run1[,2], lwd = 2, xlab = "Rounds",
     ylab = "Weight factor", xlim = c(1,7), ylim = c(0,6),
     col = "purple", main = "Weight factor at different rounds")
points(x = seq(1,5,by = 1), y = ess.bw.run2[,2], pch = 2, lwd = 2, col = "red")
points(x = seq(1,5,by = 1), y = ess.bw.run3[,2], pch = 10, lwd = 2, col = "blue")
legend(6, 6.0, legend = c("Run 1", "Run 2", "Run 3"), 
       col = c("purple", "red", "blue", "black"), pch = c(1,2,10))

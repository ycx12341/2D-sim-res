# Reference densities.R
# Author: Yunchen Xiao

# This .R file reads in the reference densities derived from the SCC patterns,
# reconstructs the reference density matrices and write them into a single .rds
# file. 

# Reset the current workspace and load the necessary packages. 
rm(list = ls())
library(readr)

# Post-day 3 reference density 
t3.dat.p1 <- read.csv("Summary pic day 3 p1.csv")
t3.dat.p2 <- read.csv("Summary pic day 3 p2.csv")
t3.dat.p3 <- read.csv("Summary pic day 3 p3.csv")
t3.dat.p4 <- read.csv("Summary pic day 3 p4.csv")
t3.dat.p5 <- read.csv("Summary pic day 3 p5.csv")
t3.dat.p6 <- read.csv("Summary pic day 3 p6.csv")
t3.dat.p7 <- read.csv("Summary pic day 3 p7.csv")

# Post-day 6 reference density
t6.dat.p1 <- read.csv("Summary pic day 6 p1.csv")
t6.dat.p2 <- read.csv("Summary pic day 6 p2.csv")
t6.dat.p3 <- read.csv("Summary pic day 6 p3.csv")
t6.dat.p4 <- read.csv("Summary pic day 6 p4.csv")
t6.dat.p5 <- read.csv("Summary pic day 6 p5.csv")
t6.dat.p6 <- read.csv("Summary pic day 6 p6.csv")
t6.dat.p7 <- read.csv("Summary pic day 6 p7.csv")

# Post-day 9 reference density
t9.dat.p1 <- read.csv("Summary pic day 9 p1.csv")
t9.dat.p2 <- read.csv("Summary pic day 9 p2.csv")
t9.dat.p3 <- read.csv("Summary pic day 9 p3.csv")
t9.dat.p4 <- read.csv("Summary pic day 9 p4.csv")
t9.dat.p5 <- read.csv("Summary pic day 9 p5.csv")
t9.dat.p6 <- read.csv("Summary pic day 9 p6.csv")
t9.dat.p7 <- read.csv("Summary pic day 9 p7.csv")

# Post-day 12 reference density
t12.dat.p1 <- read.csv("Summary pic day 12 p1.csv")
t12.dat.p2 <- read.csv("Summary pic day 12 p2.csv")
t12.dat.p3 <- read.csv("Summary pic day 12 p3.csv")
t12.dat.p4 <- read.csv("Summary pic day 12 p4.csv")
t12.dat.p5 <- read.csv("Summary pic day 12 p5.csv")
t12.dat.p6 <- read.csv("Summary pic day 12 p6.csv")
t12.dat.p7 <- read.csv("Summary pic day 12 p7.csv")

# Post-day 14 reference density
t14.dat.p1 <- read.csv("Summary pic day 14 p1.csv")
t14.dat.p2 <- read.csv("Summary pic day 14 p2.csv")
t14.dat.p3 <- read.csv("Summary pic day 14 p3.csv")
t14.dat.p4 <- read.csv("Summary pic day 14 p4.csv")
t14.dat.p5 <- read.csv("Summary pic day 14 p5.csv")
t14.dat.p6 <- read.csv("Summary pic day 14 p6.csv")
t14.dat.p7 <- read.csv("Summary pic day 14 p7.csv")

# Reconstruct the density matrices. 
t3.ref.den <- cbind(rev(t3.dat.p1[, 5]), rev(t3.dat.p2[, 5]), 
                    rev(t3.dat.p3[, 5]), rev(t3.dat.p4[, 5]), 
                    rev(t3.dat.p5[, 5]), rev(t3.dat.p6[, 5]),
                    rev(t3.dat.p7[, 5]))/100

t6.ref.den <- cbind(rev(t6.dat.p1[, 5]), rev(t6.dat.p2[, 5]), 
                    rev(t6.dat.p3[, 5]), rev(t6.dat.p4[, 5]), 
                    rev(t6.dat.p5[, 5]), rev(t6.dat.p6[, 5]),
                    rev(t6.dat.p7[, 5]))/100

t9.ref.den <- cbind(rev(t9.dat.p1[, 5]), rev(t9.dat.p2[, 5]), 
                    rev(t9.dat.p3[, 5]), rev(t9.dat.p4[, 5]), 
                    rev(t9.dat.p5[, 5]), rev(t9.dat.p6[, 5]),
                    rev(t9.dat.p7[, 5]))/100

t12.ref.den <- cbind(rev(t12.dat.p1[, 5]), rev(t12.dat.p2[, 5]), 
                     rev(t12.dat.p3[, 5]), rev(t12.dat.p4[, 5]), 
                     rev(t12.dat.p5[, 5]), rev(t12.dat.p6[, 5]),
                     rev(t12.dat.p7[, 5]))/100

t14.ref.den <- cbind(rev(t14.dat.p1[, 5]), rev(t14.dat.p2[, 5]), 
                     rev(t14.dat.p3[, 5]), rev(t14.dat.p4[, 5]), 
                     rev(t14.dat.p5[, 5]), rev(t14.dat.p6[, 5]),
                     rev(t14.dat.p7[, 5]))/100

# Write the recorded densities into a single .rds file. 
write_rds(list(t3.ref.den = t3.ref.den, t6.ref.den = t6.ref.den, 
               t9.ref.den = t9.ref.den, t12.ref.den = t12.ref.den, 
               t14.ref.den = t14.ref.den), "Reference densities SCC.rds")

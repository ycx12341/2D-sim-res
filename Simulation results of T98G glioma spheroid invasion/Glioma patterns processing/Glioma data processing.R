library(readr)

# Pattern d1 data 
glioma.d1.info.p1 <- read.csv("Summary pattern glioma t1 -40um pt1.csv")
glioma.d1.info.p2 <- read.csv("Summary pattern glioma t1 -40um pt2.csv")
glioma.d1.info.p3 <- read.csv("Summary pattern glioma t1 -40um pt3.csv")
glioma.d1.info.p4 <- read.csv("Summary pattern glioma t1 -40um pt4.csv")
glioma.d1.info.p5 <- read.csv("Summary pattern glioma t1 -40um pt5.csv")
glioma.d1.info.p6 <- read.csv("Summary pattern glioma t1 -40um pt6.csv")
glioma.d1.info.p7 <- read.csv("Summary pattern glioma t1 -40um pt7.csv")
glioma.d1.info.p8 <- read.csv("Summary pattern glioma t1 -40um pt8.csv")
glioma.d1.info.p9 <- read.csv("Summary pattern glioma t1 -40um pt9.csv")
glioma.d1.info.p10 <- read.csv("Summary pattern glioma t1 -40um pt10.csv")

glioma.d1.data <- rbind(glioma.d1.info.p1[,5]/100, glioma.d1.info.p2[,5]/100,
                        glioma.d1.info.p3[,5]/100, glioma.d1.info.p4[,5]/100,
                        glioma.d1.info.p5[,5]/100, glioma.d1.info.p6[,5]/100,
                        glioma.d1.info.p7[,5]/100, glioma.d1.info.p8[,5]/100,
                        glioma.d1.info.p9[,5]/100, glioma.d1.info.p10[,5]/100)

# Pattern d3 data
glioma.d3.info.p1 <- read.csv("Summary pattern glioma t3 -40um pt1.csv")
glioma.d3.info.p2 <- read.csv("Summary pattern glioma t3 -40um pt2.csv")
glioma.d3.info.p3 <- read.csv("Summary pattern glioma t3 -40um pt3.csv")
glioma.d3.info.p4 <- read.csv("Summary pattern glioma t3 -40um pt4.csv")
glioma.d3.info.p5 <- read.csv("Summary pattern glioma t3 -40um pt5.csv")
glioma.d3.info.p6 <- read.csv("Summary pattern glioma t3 -40um pt6.csv")
glioma.d3.info.p7 <- read.csv("Summary pattern glioma t3 -40um pt7.csv")
glioma.d3.info.p8 <- read.csv("Summary pattern glioma t3 -40um pt8.csv")
glioma.d3.info.p9 <- read.csv("Summary pattern glioma t3 -40um pt9.csv")
glioma.d3.info.p10 <- read.csv("Summary pattern glioma t3 -40um pt10.csv")

glioma.d3.data <- rbind(glioma.d3.info.p1[,5]/100, glioma.d3.info.p2[,5]/100,
                        glioma.d3.info.p3[,5]/100, glioma.d3.info.p4[,5]/100,
                        glioma.d3.info.p5[,5]/100, glioma.d3.info.p6[,5]/100,
                        glioma.d3.info.p7[,5]/100, glioma.d3.info.p8[,5]/100,
                        glioma.d3.info.p9[,5]/100, glioma.d3.info.p10[,5]/100)

glioma.ref.data <- list(t1.ref.den = glioma.d1.data, 
                        t3.ref.den = glioma.d3.data)

write_rds(glioma.ref.data, "Glioma reference data.rds")

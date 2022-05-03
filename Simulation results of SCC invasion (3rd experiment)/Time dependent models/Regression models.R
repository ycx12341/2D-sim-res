# Regression models.R 
# Author: Yunchen Xiao
# This .R file reads in the final parameter estimates for each of the 5 SCC
# patterns, constructs the regression models and decides which model gives the 
# best fit to the change of parameter estimates in different periods. 

# Reset the current workspace and load the necessary packages.
rm(list = ls())
library(readr)

# Read in the final parameter estimates for each pattern. 
paras.d3 <- read_rds("Day 3 parameter estimates.rds")
paras.d6 <- read_rds("Day 6 parameter estimates.rds")
paras.d9 <- read_rds("Day 9 parameter estimates.rds")
paras.d12 <- read_rds("Day 12 parameter estimates.rds")
paras.d14 <- read_rds("Day 14 parameter estimates.rds")
paras.d3.final <- paras.d3$paras.day3.final.est
paras.d3.final.no.init.cols <- paras.d3.final[-7]
paras.d6.final <- paras.d6$paras.day6.final.est
paras.d9.final <- paras.d9$paras.day9.final.est
paras.d12.final <- paras.d12$paras.day12.final.est
paras.d14.final <- paras.d14$paras.day14.final.est

# Sort these parameter estimates into a matrix.
post.paras.est.mat <- rbind(paras.d3.final.no.init.cols,
                            paras.d6.final, paras.d9.final,
                            paras.d12.final, paras.d14.final)
colnames(post.paras.est.mat) <- c("dn.ests", "gamma.ests", "rn.ests", 
                                  "eta.ests", "dm.ests", "alpha.ests", 
                                  "prop.death.ests", "prop.mitosis.ests")

# Set the explanatory variables to be the different periods, get ready for 
# constructing the regression models. 
periods <- c(1,2,3,4,5)
periods.2 <- periods ^ 2
post.paras.est.mat <- as.data.frame(cbind(periods, periods.2, 
                                          post.paras.est.mat))

# Regression models fitting. For each parameter's estimates, fit 3 different 
# kinds of regression models: constant, linear and quadratic. 

# dn
dn.smooth.linear <- lm(dn.ests ~ periods, data = post.paras.est.mat)
dn.smooth.inter.only <- lm(dn.ests ~ 1, data = post.paras.est.mat)
dn.smooth.quad <- lm(dn.ests ~ periods + periods.2, data = post.paras.est.mat)

summary(dn.smooth.linear) # Adj-R^2 = 0.6533, linear coefficient is weakly significant
summary(dn.smooth.inter.only) # Coefficient is not significant.
summary(dn.smooth.quad) # Adj-R^2 = 0.7407, none of the parameters' coefficient is significant. 

# Adjusted R-squared values and the p-values of coefficients suggest the linear fit without 
# intercept is better. 

AIC(dn.smooth.linear) # -47.95767
AIC(dn.smooth.inter.only) # -43.22211
AIC(dn.smooth.quad) # -49.43654

# AIC results suggest a quadratic fit is the best.

periods.cont <- seq(0,5, by = 0.01)
fitted.vals.predict.dn.linear <- rep(0, length = length(periods.cont))
fitted.vals.predict.dn.inter.only <- rep(0, length = length(periods.cont))
fitted.vals.predict.dn.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.dn.linear[i] <- predict(dn.smooth.linear,
                              newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.dn.inter.only[i] <- predict(dn.smooth.inter.only,
                              newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.dn.quad[i] <- predict(dn.smooth.quad,
                              newdata = data.frame(periods = periods.cont[i],
                              periods.2 = periods.cont[i]^2))
}

plot(periods, post.paras.est.mat$dn.ests, pch = 21, cex = 2, col = "purple",
     bg = "purple", xlab = "Periods", 
     ylab = "Estimates of dn", main = "Smoothing fit of dn estimates")
lines(periods.cont, fitted.vals.predict.dn.linear, col = "red", lwd = 2)
lines(periods.cont, fitted.vals.predict.dn.inter.only, col = "blue", lwd = 2)
lines(periods.cont, fitted.vals.predict.dn.quad, col = "black", lwd = 2)

# A quadratic fit should be the best option for dn! 

write_rds(dn.smooth.quad, "dn quadratic regression model.rds")

# gamma
gamma.smooth.linear <- lm(gamma.ests ~ periods, data = post.paras.est.mat)
gamma.smooth.inter.only <- lm(gamma.ests ~ 1, data = post.paras.est.mat)
gamma.smooth.quad <- lm(gamma.ests ~ periods + periods.2, 
                        data = post.paras.est.mat)

summary(gamma.smooth.linear) # Adj-R^2 = 0.4572, none of the coefficients is significant. 
summary(gamma.smooth.inter.only) # Moderate evidence that the coefficient is significant.
summary(gamma.smooth.quad) # Adj-R^2 = 0.8757, linear and quadratic coefficients are weakly significant. 

# Adjusted R-squared values and the p-values of the coefficients suggest we should either
# choose the linear fit or the quadratic fit. 

AIC(gamma.smooth.linear) # -14.15931
AIC(gamma.smooth.inter.only) # -11.66556
AIC(gamma.smooth.quad) # -21.55722

# AIC values suggest the quadratic fit is the best option. 

fitted.vals.predict.gamma.linear <- rep(0, length = length(periods.cont))
fitted.vals.predict.gamma.inter.only <- rep(0, length = length(periods.cont))
fitted.vals.predict.gamma.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.gamma.linear[i] <- predict(gamma.smooth.linear,
                               newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.gamma.inter.only[i] <- predict(gamma.smooth.inter.only,
                               newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.gamma.quad[i] <- predict(gamma.smooth.quad, 
                               newdata = data.frame(periods = periods.cont[i],
                               periods.2 = periods.cont[i]^2))
  
}


plot(periods, post.paras.est.mat$gamma.ests, xlab = "Periods", 
     ylab = "Estimates of gamma", pch = 21, cex = 2, bg = "purple", 
     main = "Smoothing fit of gamma estimates", ylim = c(0, 0.16))
lines(periods.cont, fitted.vals.predict.gamma.linear, col = "red", lwd = 2)
lines(periods.cont, fitted.vals.predict.gamma.inter.only, col = "blue", lwd = 2)
lines(periods.cont, fitted.vals.predict.gamma.quad, col = "black", lwd = 2)

# A quadratic fit should be the best option! 

write_rds(gamma.smooth.quad, "gamma quadratic regression model.rds")

# rn estimates
rn.smooth.linear <- lm(rn.ests ~ periods,
                       data = post.paras.est.mat)
rn.smooth.inter.only <- lm(rn.ests ~ 1, data = post.paras.est.mat)
rn.smooth.quad <- lm(rn.ests ~ periods + periods.2, data = post.paras.est.mat)

summary(rn.smooth.linear) # 0.3457, intercept is moderately significant.
summary(rn.smooth.inter.only) # Coefficient is strongly significant. 
summary(rn.smooth.quad) # 0.7227, none of the coefficient is significant.  

AIC(rn.smooth.linear) # -31.68584
AIC(rn.smooth.inter.only) # -30.12646
AIC(rn.smooth.quad) # -36.00461

# AIC results suggest quadratic fit is the best option. 

fitted.vals.predict.rn.linear <- rep(0, length = length(periods.cont))
fitted.vals.predict.rn.inter.only<- rep(0, length = length(periods.cont))
fitted.vals.predict.rn.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.rn.linear[i] <- predict(rn.smooth.linear,
                               newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.rn.inter.only[i] <- predict(rn.smooth.inter.only,
                               newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.rn.quad[i] <- predict(rn.smooth.quad, 
                               newdata = data.frame(periods = periods.cont[i],
                               periods.2 = periods.cont[i]^2))
}



plot(periods, post.paras.est.mat$rn.ests, xlab = "Periods", 
     ylab = "Estimates of rn", pch = 21, bg = "purple", cex = 2,
     main = "Smoothing fit of rn estimates")
lines(periods.cont, fitted.vals.predict.rn.linear, col = "red", lwd = 2)
lines(periods.cont, fitted.vals.predict.rn.inter.only, col = "blue", lwd = 2)
lines(periods.cont, fitted.vals.predict.rn.quad, col = "black", lwd = 2)

# A quadratic fit would be the best option. 

write_rds(rn.smooth.quad, "rn quadratic regression model.rds")

# eta 
eta.smooth.linear <- lm(eta.ests ~ periods,
                       data = post.paras.est.mat)
eta.smooth.inter.only <- lm(eta.ests ~ 1, data = post.paras.est.mat)
eta.smooth.quad <- lm(eta.ests ~ periods + periods.2, data = post.paras.est.mat)

summary(eta.smooth.linear) # 0.2614, intercept term is strongly significant.
summary(eta.smooth.inter.only) # Intercept is strongly significant. 
summary(eta.smooth.quad) # 0.5901, intercept term is weakly significant.

AIC(eta.smooth.linear) # 25.01037
AIC(eta.smooth.inter.only) # 25.96406
AIC(eta.smooth.quad) # 22.03963

# AIC results suggest the quadratic fit with an intercept term is the best fit. 
 
fitted.vals.predict.eta.linear <- rep(0, length = length(periods.cont))
fitted.vals.predict.eta.inter.only <- rep(0, length = length(periods.cont))
fitted.vals.predict.eta.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.eta.linear[i] <- predict(eta.smooth.linear,
                              newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.eta.inter.only[i] <- predict(eta.smooth.inter.only,
                              newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.eta.quad[i] <- predict(eta.smooth.quad, 
                                    data.frame(periods = periods.cont[i],
                                    periods.2 = periods.cont[i]^2))
}

plot(periods, post.paras.est.mat$eta.ests, xlab = "Periods", 
     ylim = c(11, 18), ylab = "Estimates of eta",
     pch = 21, cex = 2, bg = "purple",
     main = "Smoothing fit of eta estimates")
lines(periods.cont, fitted.vals.predict.eta.linear, col = "red", lwd = 2)
lines(periods.cont, fitted.vals.predict.eta.inter.only, col = "blue", 
      lwd = 2)
lines(periods.cont, fitted.vals.predict.eta.quad, col = "black", lwd = 2)

# Quadratic fit should be considered as the best option. 

write_rds(eta.smooth.quad, "eta quadratic regression model.rds")

# dm
dm.smooth.linear <- lm(dm.ests ~ periods,
                        data = post.paras.est.mat)
dm.smooth.inter.only <- lm(dm.ests ~ 1, data = post.paras.est.mat)
dm.smooth.quad <- lm(dm.ests ~ periods + periods.2, data = post.paras.est.mat)

summary(dm.smooth.linear) # -0.2934, intercept term is strongly significant. 
summary(dm.smooth.inter.only) # Strong evidence that the coefficient is significant.
summary(dm.smooth.quad) # -0.826, intercept term is weakly significant.  

AIC(dm.smooth.linear) # -48.00336
AIC(dm.smooth.inter.only) # -49.85147
AIC(dm.smooth.quad) # -46.30652

# AIC results suggest a constant fit is the best.  

fitted.vals.predict.dm.linear <- rep(0, length = length(periods.cont))
fitted.vals.predict.dm.inter.only <- rep(0, length = length(periods.cont))
fitted.vals.predict.dm.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.dm.linear[i] <- predict(dm.smooth.linear,
                              newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.dm.inter.only[i] <- predict(dm.smooth.inter.only,
                              newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.dm.quad[i] <- predict(dm.smooth.quad, 
                                newdata = data.frame(periods = periods.cont[i],
                                periods.2 = periods.cont[i]^2))
}

plot(periods, post.paras.est.mat$dm.ests, xlab = "Periods", 
     ylab = "Estimates of dm", cex = 2, bg = "purple", 
     pch = 21, 
     main = "Smoothing fit of dm estimates")
lines(periods.cont, fitted.vals.predict.dm.linear, col = "red", lwd = 2)
lines(periods.cont, fitted.vals.predict.dm.inter.only, col = "blue", 
      lwd = 2)
lines(periods.cont, fitted.vals.predict.dm.quad, col = "black", lwd = 2)

# Constant fit should be the best option. 

write_rds(dm.smooth.inter.only, "dm constant regression model.rds")

# alpha
alpha.smooth.linear <- lm(alpha.ests ~ periods,
                       data = post.paras.est.mat)
alpha.smooth.inter.only <- lm(alpha.ests ~ 1, data = post.paras.est.mat)
alpha.smooth.quad <- lm(alpha.ests ~ periods + periods.2, 
                        data = post.paras.est.mat)


summary(alpha.smooth.linear) # 0.3332, strong evidence that the intercept term is significant.
summary(alpha.smooth.inter.only) # Strong evidence that the intercept term is significant.
summary(alpha.smooth.quad) # 0.7791, the intercept term is strongly significant, the linea term is weakly significant. 

# Adjusted R-squared values and p-values suggest we should choose either the 
# quadratic fit or the constant fit. 

AIC(alpha.smooth.linear) # -23.12725
AIC(alpha.smooth.inter.only) # -21.66243
AIC(alpha.smooth.quad) # -28.67872

# AIC results suggest the quadratic fit is the best among all 3.

fitted.vals.predict.alpha.linear <- rep(0, length = length(periods.cont))
fitted.vals.predict.alpha.inter.only <- rep(0, length = length(periods.cont))
fitted.vals.predict.alpha.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.alpha.linear[i] <- predict(alpha.smooth.linear,
                               newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.alpha.inter.only[i] <- predict(alpha.smooth.inter.only,
                               newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.alpha.quad[i] <- predict(alpha.smooth.quad, 
                                newdata = data.frame(periods = periods.cont[i],
                                                periods.2 = periods.cont[i]^2))
}

plot(periods, post.paras.est.mat$alpha.ests, xlab = "Periods", 
     ylab = "Estimates of alpha", ylim = c(0.09, 0.18), cex = 2, bg = "purple",
     pch = 21, main = "Smoothing fit of alpha estimates")
lines(periods.cont, fitted.vals.predict.alpha.linear, col = "red", lwd = 2)
lines(periods.cont, fitted.vals.predict.alpha.inter.only, col = "blue", 
      lwd = 2)
lines(periods.cont, fitted.vals.predict.alpha.quad, col = "black", lwd = 2)

# Seems like quadratic is the best.

write_rds(alpha.smooth.quad, "alpha quadratic regression model.rds")

# Prop.extinct
prop.ex.smooth.linear <- lm(prop.death.ests ~ periods,
                          data = post.paras.est.mat)
prop.ex.smooth.inter.only <- lm(prop.death.ests ~ 1,
                                   data = post.paras.est.mat)
prop.ex.smooth.quad <- lm(prop.death.ests ~ periods + periods.2, 
                          data = post.paras.est.mat)

summary(prop.ex.smooth.linear) # -0.2809, intercept term is moderately significant.
summary(prop.ex.smooth.inter.only) # Strong evidence that the intercept term is significant.
summary(prop.ex.smooth.quad) # -0.6349, none of the coefficient is significant.

AIC(prop.ex.smooth.linear) # -28.65545
AIC(prop.ex.smooth.inter.only) # -30.45477
AIC(prop.ex.smooth.quad) # -27.46268

# AIC results suggest we should choose the constant fit 

fitted.vals.predict.prop.ex.linear <- rep(0, length = length(periods.cont))
fitted.vals.predict.prop.ex.inter.only <- rep(0, 
                                              length = length(periods.cont))
fitted.vals.predict.prop.ex.quad <- rep(0, length = length(periods.cont))

for (i in 1:length(periods.cont)) {
  fitted.vals.predict.prop.ex.linear[i] <- predict(prop.ex.smooth.linear,
                              newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.prop.ex.inter.only[i] <- predict(prop.ex.smooth.inter.only,
                              newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.prop.ex.quad[i] <- predict(prop.ex.smooth.quad, 
                                newdata = data.frame(periods = periods.cont[i],
                                periods.2 = periods.cont[i]^2))
  
}

plot(periods, post.paras.est.mat$prop.death.ests, xlab = "Periods", 
     ylab = "Estimates of prop.extinct", pch = 21, cex = 2, bg = "purple",
     main = "Smoothing fit of prop.extinct estimates")
lines(periods.cont, fitted.vals.predict.prop.ex.linear, col = "red", lwd = 2)
lines(periods.cont, fitted.vals.predict.prop.ex.inter.only, col = "blue", 
      lwd = 2)
lines(periods.cont, fitted.vals.predict.prop.ex.quad, col = "black", lwd = 2)

# A constant fit would be the best! 

write_rds(prop.ex.smooth.inter.only, "prob.death constant regression model.rds")

# Prop.prof.

prop.prof.smooth.linear <- lm(prop.mitosis.ests ~ periods,
                            data = post.paras.est.mat)
prop.prof.smooth.inter.only <- lm(prop.mitosis.ests ~ 1,
                                     data = post.paras.est.mat)
prop.prof.smooth.quad <- lm(prop.mitosis.ests ~ periods + periods.2, 
                            data = post.paras.est.mat)


summary(prop.prof.smooth.linear) # -0.2773, intercept term is weakly significant.
summary(prop.prof.smooth.inter.only) # Coefficient is strongly significant. 
summary(prop.prof.smooth.quad) # 0.9274, moderate evidence that the linear and quadratic terms are significant. 


AIC(prop.prof.smooth.linear) # 4.709535
AIC(prop.prof.smooth.inter.only) # 2.924205
AIC(prop.prof.smooth.quad) # -9.6583

# AIC results suggest the quadratic fit is the most suitable one among all 3. 

fitted.vals.predict.prop.prof.linear <- rep(0, length = length(periods.cont))
fitted.vals.predict.prop.prof.inter.only <- rep(0, 
                                                length = length(periods.cont))
fitted.vals.predict.prop.prof.quad <- rep(0, length = length(periods.cont))


for (i in 1:length(periods.cont)) {
  fitted.vals.predict.prop.prof.linear[i] <- predict(prop.prof.smooth.linear,
                                                   newdata = 
                                         data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.prop.prof.inter.only[i] <- predict(prop.prof.smooth.inter.only,
                               newdata = data.frame(periods = periods.cont[i]))
  
  fitted.vals.predict.prop.prof.quad[i] <- predict(prop.prof.smooth.quad, 
                                newdata = data.frame(periods = periods.cont[i],
                                periods.2 = periods.cont[i]^2))
}

plot(periods, post.paras.est.mat$prop.mitosis.ests, xlab = "Periods", 
     ylab = "Estimates of prop.mitosis", pch = 21, cex = 2, bg = "purple",
     main = "Smoothing fit of prop.mitosis estimates")
lines(periods.cont, fitted.vals.predict.prop.prof.linear, col = "red", lwd = 2)
lines(periods.cont, fitted.vals.predict.prop.prof.inter.only, col = "blue", 
      lwd = 2)
lines(periods.cont, fitted.vals.predict.prop.prof.quad, col = "black", 
      lwd = 2)

write_rds(prop.prof.smooth.quad, "prob.prof quadratic regression model.rds")

# Determining sGARCH for normal distribution using AIC
aic_garch_norm = matrix(0,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      garch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                    model = "sGARCH"), 
                              mean.model = list(armaOrder=c(0,0),
                                                include.mean = TRUE),
                              distribution.model = "norm")
      garch_fit = ugarchfit(spec=garch_spec, data=returns, 
                            solver = "gosolnp",
                            solver.control=list(trace = 1))
      aic_garch_norm[i,j] = infocriteria(garch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

aic_garch_norm
which(min(aic_garch_norm) == aic_garch_norm)

par(mar=c(5, 5, 4, 5))
plot(aic_garch_norm, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, GARCH, N(0,1)",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off


# Model specification
garchspec_norm = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                include.mean = TRUE), 
                              variance.model = list(garchOrder = c(1,1), 
                                                    model = "sGARCH"),
                              distribution.model = "nig")
garchfit_norm = ugarchfit(data = returns, spec = garchspec_norm, solver = "gosolnp",
                            solver.control=list(trace = 1))
garchfit_norm

plot(garchfit_norm)

# Model validation
autoplot(acf(garchfit_norm@fit$z), main = "ACF of the standardized residuals") + # acf of the standardized residuals
  autoplot(acf(garchfit_norm@fit$z^2), main = "ACF of the standardized squared residuals") # acf of (standardized residuals)^2 

par(mfrow = c(1,2))
plot(garchfit_norm, which = 9) # QQ plot of the standardized residuals
plot(garchfit_norm, which = 8)
dev.off


#-------------------------------------------------------------------------------

# Determining sGARCH for sstd distribution using AIC
aic_garch_sstd = matrix(0,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      garch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                    model = "sGARCH"), 
                              mean.model = list(armaOrder=c(0,0),
                                                include.mean = FALSE),
                              distribution.model = "sstd")
      garch_fit = ugarchfit(spec=garch_spec, data=returns, 
                            solver = "gosolnp",
                            solver.control=list(trace = 1))
      aic_garch_sstd[i,j] = infocriteria(garch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

aic_garch_sstd
which(min(aic_garch_sstd) == aic_garch_sstd)

par(mar=c(5, 5, 4, 5))
plot(aic_garch_sstd, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, GARCH, SSTD(0,1)",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off


# Model specification
garchspec_sstd = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                include.mean = TRUE), 
                              variance.model = list(garchOrder = c(1,1), 
                                                    model = "sGARCH"),
                              distribution.model = "norm")
garchfit_sstd = ugarchfit(data = returns, spec = garchspec_sstd, solver = "gosolnp",
                            solver.control=list(trace = 1))
garchfit_sstd


garchfit_sstd

Weighted.Box.test((garchfit_sstd@fit$z)^2, lag = 1, type = "Ljung-Box", fitdf = 0)


# Model validation
autoplot(acf(garchfit_sstd@fit$z), main = "ACF of the standardized residuals") + # acf of the standardized residuals
  autoplot(acf(garchfit_sstd@fit$z^2), main = "ACF of the standardized squared residuals") # acf of (standardized residuals)^2 

par(mfrow = c(1,2))
plot(garch22fit_sstd, which = 9) # QQ plot of the standardized residuals
plot(garch22fit_sstd, which = 8)
dev.off
#-------------------------------------------------------------------------------

# Determining sGARCH for ged distribution using AIC
aic_garch_ged = matrix(0,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      garch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                    model = "sGARCH"), 
                              mean.model = list(armaOrder=c(0,0),
                                                include.mean = FALSE),
                              distribution.model = "ged")
      garch_fit = ugarchfit(spec=garch_spec, data=returns, 
                            solver = "gosolnp",
                            solver.control=list(trace = 1))
      aic_garch_ged[i,j] = infocriteria(garch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

aic_garch_ged
which(min(aic_garch_ged) == aic_garch_ged)

par(mar=c(5, 5, 4, 5))
plot(aic_garch_ged, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, GARCH models with generalized error distribution",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off


# Model specification
garchspec_ged = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                               include.mean = TRUE), 
                             variance.model = list(garchOrder = c(1,1), 
                                                   model = "sGARCH"),
                             distribution.model = "ged")
garchfit_ged = ugarchfit(data = returns, spec = garchspec_ged, solver = "gosolnp",
                           solver.control=list(trace = 1))
garchfit_ged

# Model validation
par(mfrow = c(1,2))
acf(garchfit_ged@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(garchfit_ged@fit$z^2,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
dev.off
par(mfrow = c(1,2))
plot(garchfit_ged, which = 9) # QQ plot of the standardized residuals
plot(garchfit_ged, which = 8)
dev.off
# --------------------------------------------------------------------------
# Leverage effect
ccf_norm = ccf(garch22fit_norm@fit$sigma, as.ts(returns))
ccf_sstd = ccf(garchfit_sstd@fit$sigma, as.ts(returns))
ccf_ged = ccf(garch22fit_ged@fit$sigma, as.ts(returns))

layout(matrix(c(0,1,1,0,2,2,3,3), 2, 4, byrow = TRUE))
plot(ccf_norm[ccf_norm$lag[36:71],], main = expression(paste(hat(rho)(sigma[t], epsilon[t-h]), " where " , eta[t], " is standard normal")), las =1)
autoplot(plot(ccf_sstd[ccf_sstd$lag[36:71],]), main = expression(paste(hat(rho)(sigma[t], epsilon[t-h]), " where " , eta[t], " is SSTD")), las =1)
plot(ccf_ged[ccf_ged$lag[36:71],], main = expression(paste(hat(rho)(sigma[t], epsilon[t-h]), " where " , eta[t], " is GED")), las =1)

autoplot(ccf_sstd, main = "CCF between the volatility and the returns") + xlim(c(0,36))

# Simulations
# SSTD
sim = ugarchsim(fit = garch22fit_sstd, n.sim = length(returns), n.start = 0, m.sim = 1000, rseed = NA)
kurt = kurtosis(sim@simulation$seriesSim)
mean(t(kurt))

# Normal
sim = ugarchsim(fit = garch22fit_norm, n.sim = length(returns), n.start = 0, m.sim = 1000, rseed = NA)
kurt = kurtosis(sim@simulation$seriesSim)
mean(t(kurt))

# GED
sim = ugarchsim(fit = garch22fit_ged, n.sim = length(returns), n.start = 0, m.sim = 1000, rseed = NA)
kurt = kurtosis(sim@simulation$seriesSim)
mean(t(kurt))


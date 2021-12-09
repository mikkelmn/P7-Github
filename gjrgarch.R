# Determining eGARCH for normal distribution using AIC
bic_gjrgarch_norm = matrix(0,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      gjrgarch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                     model = "gjrGARCH"), 
                               mean.model = list(armaOrder=c(0,0),
                                                 include.mean = TRUE),
                               distribution.model = "norm")
      gjrgarch_fit = ugarchfit(spec=gjrgarch_spec, data=returns, 
                             solver = "gosolnp",
                             solver.control=list(trace = 1))
      bic_gjrgarch_norm[i,j] = infocriteria(gjrgarch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

bic_gjrgarch_norm
which(min(bic_gjrgarch_norm) == bic_gjrgarch_norm)

par(mar=c(5, 5, 4, 5))
plot(bic_gjrgarch_norm, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, GJR-GARCH, N(0,1)",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off

# Model specification
gjrgarchspec_norm = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                 include.mean = TRUE), 
                               variance.model = list(garchOrder = c(1,1), 
                                                     model = "gjrGARCH"),
                               distribution.model = "norm")
gjrgarchfit_norm = ugarchfit(data = returns, spec = gjrgarchspec_norm, solver = "gosolnp",
                             solver.control=list(trace = 1))
gjrgarchfit_norm

# Model validation
par(mfrow = c(1,2))
acf(gjrgarchfit_norm@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(gjrgarchfit_norm@fit$z,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
dev.off
par(mfrow = c(1,2))
plot(gjrgarchfit_norm, which = 9) # QQ plot of the standardized residuals
plot(gjrgarchfit_norm, which = 8)
dev.off

#-------------------------------------------------------------------------------

# Determining GJR-GARCH for sstd distribution using AIC
bic_gjrgarch_sstd = matrix(0,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      gjrgarch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                       model = "gjrGARCH"), 
                                 mean.model = list(armaOrder=c(0,0),
                                                   include.mean = TRUE),
                                 distribution.model = "sstd")
      gjrgarch_fit = ugarchfit(spec=gjrgarch_spec, data=returns, 
                               solver = "gosolnp",
                               solver.control=list(trace = 1))
      bic_gjrgarch_sstd[i,j] = infocriteria(gjrgarch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

bic_gjrgarch_sstd
which(min(bic_gjrgarch_sstd) == bic_gjrgarch_sstd)

par(mar=c(5, 5, 4, 5))
plot(bic_gjrgarch_sstd, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, GJR-GARCH, SSTD(0,1)",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off



# Model specification
gjrgarchspec_sstd = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                 include.mean = TRUE), 
                               variance.model = list(garchOrder = c(1,1), 
                                                     model = "gjrGARCH"),
                               distribution.model = "sstd")
gjrgarchfit_sstd = ugarchfit(data = returns, spec = gjrgarchspec_sstd, solver = "gosolnp",
                             solver.control=list(trace = 1))
gjrgarchfit_sstd

# Model validation
par(mfrow = c(1,2))
acf(gjrgarchfit_sstd@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(gjrgarchfit_sstd@fit$z,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
dev.off
par(mfrow = c(1,2))
plot(gjrgarchfit_sstd, which = 9) # QQ plot of the standardized residuals
plot(gjrgarchfit_sstd, which = 8)
dev.off

#-------------------------------------------------------------------------------

# Determining GJR-GARCH for ged distribution using BIC-----
bic_gjrgarch_ged = matrix(0,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      gjrgarch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                     model = "gjrGARCH"), 
                               mean.model = list(armaOrder=c(0,0),
                                                 include.mean = TRUE),
                               distribution.model = "ged")
      gjrgarch_fit = ugarchfit(spec=gjrgarch_spec, data=returns, 
                             solver = "gosolnp",
                             solver.control=list(trace = 1))
      bic_gjrgarch_ged[i,j] = infocriteria(gjrgarch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

bic_gjrgarch_ged
which(min(bic_gjrgarch_ged) == bic_gjrgarch_ged)

par(mar=c(5, 5, 4, 5))
plot(bic_gjrgarch_ged, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, GJR-GARCH models with GED",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off


# Model specification-------
gjrgarchspec_ged = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                include.mean = TRUE), 
                              variance.model = list(garchOrder = c(1,1), 
                                                    model = "gjrGARCH"),
                              distribution.model = "ged")
gjrgarchfit_ged = ugarchfit(data = returns, spec = gjrgarchspec_ged, solver = "gosolnp",
                            solver.control=list(trace = 1))
gjrgarchfit_ged
plot.ts(gjrgarchfit_ged@fit$z)
# Model validation
par(mfrow = c(1,2))
acf(gjrgarchfit_ged@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(gjrgarchfit_ged@fit$z,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
dev.off
par(mfrow = c(1,2))
plot(gjrgarchfit_ged, which = 9) # QQ plot of the standardized residuals
plot(gjrgarchfit_ged, which = 8)
dev.off

# Simulations ----
# SSTD
sim = ugarchsim(fit = gjrgarch31fit_sstd, n.sim = length(returns), n.start = 0, m.sim = 1000, rseed = NA)
kurt = kurtosis(sim@simulation$seriesSim)
mean(t(kurt))

# Normal
sim = ugarchsim(fit = gjrgarch31fit_norm, n.sim = length(returns), n.start = 0, m.sim = 1000, rseed = NA)
kurt = kurtosis(sim@simulation$seriesSim)
mean(t(kurt))

# GED
sim = ugarchsim(fit = gjrgarch31fit_ged, n.sim = length(returns), n.start = 0, m.sim = 1000, rseed = NA)
kurt = kurtosis(sim@simulation$seriesSim)
mean(t(kurt))




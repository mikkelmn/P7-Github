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
     main = "BIC, GJR-GARCH models with Gaussian Distribution",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off

# Model specification
gjrgarch31spec_norm = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                 include.mean = TRUE), 
                               variance.model = list(garchOrder = c(3,1), 
                                                     model = "gjrGARCH"),
                               distribution.model = "norm")
gjrgarch31fit_norm = ugarchfit(data = returns, spec = gjrgarch31spec_norm, solver = "gosolnp",
                             solver.control=list(trace = 1))
gjrgarch31fit_norm

# Model validation
par(mfrow = c(1,2))
acf(gjrgarch31fit_norm@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(gjrgarch31fit_norm@fit$z,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
par(mfrow = c(1,2))
plot(gjrgarch31fit_norm, which = 9) # QQ plot of the standardized residuals
plot(gjrgarch31fit_norm, which = 8)

#-------------------------------------------------------------------------------

# Determining GJR-GARCH for sstd distribution using AIC
bic_gjrgarch_norm = matrix(0,5,5)
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
      bic_gjrgarch_norm[i,j] = infocriteria(gjrgarch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

bic_gjrgarch_norm
which(min(bic_gjrgarch_norm) == bic_gjrgarch_norm)

par(mar=c(5, 5, 4, 5))
plot(bic_gjrgarch_norm, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, GJR-GARCH models with SSTD",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off



# Model specification
gjrgarch31spec_sstd = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                 include.mean = TRUE), 
                               variance.model = list(garchOrder = c(3,1), 
                                                     model = "gjrGARCH"),
                               distribution.model = "sstd")
gjrgarch31fit_sstd = ugarchfit(data = returns, spec = gjrgarch31spec_sstd, solver = "gosolnp",
                             solver.control=list(trace = 1))
gjrgarch31fit_sstd

# Model validation
par(mfrow = c(1,2))
acf(gjrgarch31fit_sstd@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(gjrgarch31fit_sstd@fit$z,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
par(mfrow = c(1,2))
plot(gjrgarch31fit_sstd, which = 9) # QQ plot of the standardized residuals
plot(gjrgarch31fit_sstd, which = 8)

#-------------------------------------------------------------------------------

# Determining GJR-GARCH for ged distribution using AIC
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


# Model specification
gjrgarch31spec_ged = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                include.mean = TRUE), 
                              variance.model = list(garchOrder = c(3,1), 
                                                    model = "gjrGARCH"),
                              distribution.model = "ged")
gjrgarch31fit_ged = ugarchfit(data = returns, spec = gjrgarch31spec_ged, solver = "gosolnp",
                            solver.control=list(trace = 1))
gjrgarch31fit_ged

# Model validation
par(mfrow = c(1,2))
acf(gjrgarch31fit_ged@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(gjrgarch31fit_ged@fit$z,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
par(mfrow = c(1,2))
plot(gjrgarch31fit_ged, which = 9) # QQ plot of the standardized residuals
plot(gjrgarch31fit_ged, which = 8)


# Simulations
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




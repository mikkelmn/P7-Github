# Determining tGARCH for normal distribution using AIC
bic_tgarch_norm = matrix(0,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      tgarch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                     model = "fGARCH", submodel = "TGARCH"), 
                               mean.model = list(armaOrder=c(0,0),
                                                 include.mean = TRUE),
                               distribution.model = "norm")
      tgarch_fit = ugarchfit(spec=tgarch_spec, data=returns, 
                             solver = "gosolnp",
                             solver.control=list(trace = 1))
      bic_tgarch_norm[i,j] = infocriteria(tgarch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

bic_tgarch_norm
which(min(bic_tgarch_norm) == bic_tgarch_norm)

par(mar=c(5, 5, 4, 5))
plot(bic_tgarch_norm, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, TGARCH, N(0,1)",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off

# Model specification
tgarchspec_norm = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                 include.mean = TRUE), 
                               variance.model = list(garchOrder = c(1,1), 
                                                     model = "fGARCH", submodel = "TGARCH"),
                               distribution.model = "norm")
tgarchfit_norm = ugarchfit(data = returns, spec = tgarchspec_norm, solver = "gosolnp",
                             solver.control=list(trace = 1))
tgarchfit_norm
plot(tgarchfit_norm)

# Model validation
par(mfrow = c(1,2))
acf(tgarchfit_norm@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(tgarchfit_norm@fit$z,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
dev.off
par(mfrow = c(1,2))
plot(tgarchfit_norm, which = 9) # QQ plot of the standardized residuals
plot(tgarchfit_norm, which = 8)
dev.off
#-------------------------------------------------------------------------------

# Determining tGARCH for sstd distribution using BIC
bic_tgarch_sstd = matrix(0,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      tgarch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                     model = "fGARCH", submodel = "TGARCH"), 
                               mean.model = list(armaOrder=c(0,0),
                                                 include.mean = TRUE),
                               distribution.model = "sstd")
      tgarch_fit = ugarchfit(spec=tgarch_spec, data=returns, 
                             solver = "gosolnp",
                             solver.control=list(trace = 1))
      bic_tgarch_sstd[i,j] = infocriteria(tgarch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

bic_tgarch_sstd
which(min(bic_tgarch_sstd) == bic_tgarch_sstd)

par(mar=c(5, 5, 4, 5))
plot(bic_tgarch_sstd, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, TGARCH, SSTD(0,1)",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off


# Model specification
tgarchspec_sstd = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                 include.mean = TRUE), 
                               variance.model = list(garchOrder = c(1,1), 
                                                     model = "fGARCH", submodel = "TGARCH"),
                               distribution.model = "sstd")
tgarchfit_sstd = ugarchfit(data = returns, spec = tgarchspec_sstd, solver = "gosolnp",
                             solver.control=list(trace = 1))
tgarchfit_sstd

plot(tgarchfit_sstd)
# Model validation
par(mfrow = c(1,2))
acf(tgarchfit_sstd@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(tgarchfit_sstd@fit$z,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
dev.off
par(mfrow = c(1,2))
plot(tgarchfit_sstd, which = 9) # QQ plot of the standardized residuals
plot(tgarchfit_sstd, which = 8)
dev.off
#-------------------------------------------------------------------------------

# Determining eGARCH for ged distribution using BIC
bic_tgarch_ged = matrix(0,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      tgarch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                     model = "fGARCH", submodel = "TGARCH"), 
                               mean.model = list(armaOrder=c(0,0),
                                                 include.mean = TRUE),
                               distribution.model = "ged")
      tgarch_fit = ugarchfit(spec=tgarch_spec, data=returns, 
                             solver = "gosolnp",
                             solver.control=list(trace = 1))
      bic_tgarch_ged[i,j] = infocriteria(tgarch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

bic_tgarch_ged
which(min(bic_tgarch_ged) == bic_tgarch_ged)

par(mar=c(5, 5, 4, 5))
plot(bic_tgarch_ged, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, TGARCH models with generalized error distribution",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off


# Model specification
tgarchspec_ged = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                include.mean = TRUE), 
                              variance.model = list(garchOrder = c(1,1), 
                                                    model = "fGARCH", submodel = "TGARCH"),
                              distribution.model = "ged")
tgarchfit_ged = ugarchfit(data = returns, spec = tgarchspec_ged, solver = "gosolnp",
                            solver.control=list(trace = 1))
tgarchfit_ged

# Model validation
par(mfrow = c(1,2))
acf(tgarchfit_ged@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(tgarchfit_ged@fit$z,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
dev.off
par(mfrow = c(1,2))
plot(tgarchfit_ged, which = 9) # QQ plot of the standardized residuals
plot(tgarchfit_ged, which = 8)
dev.off

# Simulations
# SSTD
sim = ugarchsim(fit = tgarch11fit_sstd, n.sim = length(returns), n.start = 0, m.sim = 1000, rseed = NA)
kurt = kurtosis(sim@simulation$seriesSim)
mean(t(kurt))

# Normal
sim = ugarchsim(fit = tgarch41fit_norm, n.sim = length(returns), n.start = 0, m.sim = 1000, rseed = NA)
kurt = kurtosis(sim@simulation$seriesSim)
mean(t(kurt))

# GED
sim = ugarchsim(fit = tgarch21fit_ged, n.sim = length(returns), n.start = 0, m.sim = 1000, rseed = NA)
kurt = kurtosis(sim@simulation$seriesSim)
mean(t(kurt))

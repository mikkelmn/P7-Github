# Determining eGARCH for normal distribution using AIC
bic_egarch_norm = matrix(0,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      egarch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                    model = "eGARCH"), 
                              mean.model = list(armaOrder=c(0,0),
                                                include.mean = TRUE),
                              distribution.model = "norm")
      egarch_fit = ugarchfit(spec=egarch_spec, data=returns, 
                            solver = "gosolnp",
                            solver.control=list(trace = 1))
      bic_egarch_norm[i,j] = infocriteria(egarch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

bic_egarch_norm
which(min(bic_egarch_norm) == bic_egarch_norm)

par(mar=c(5, 5, 4, 5))
plot(bic_egarch_norm, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, EGARCH, N(0,1)",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off

# Model specification
egarch33spec_norm = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                include.mean = TRUE), 
                              variance.model = list(garchOrder = c(1,1), 
                                                    model = "eGARCH"),
                              distribution.model = "norm")
egarch33fit_norm = ugarchfit(data = returns, spec = egarch33spec_norm, solver = "gosolnp",
                            solver.control=list(trace = 1))
plot(as.xts(egarch33fit_norm@fit$sigma))
egarch33fit_norm@fit$sigma

vix1 = data.frame(VIX)
vix = vix1 %>% transmute(gjrgarchfit_sstd@fit$sigma)
sigma1 = as.xts(vix)
plot(sigma1, main = "GJR-GARCH, SSTD(0,1)", xlab = "", yaxis.right = FALSE, ylim = c(0.00,0.085), format.labels="%b\n%Y")
dev.off
par(mfrow = c(4,2))

mean(garchfit_norm@fit$sigma-garchfit_sstd@fit$sigma)
mean(egarch33fit_norm@fit$sigma-egarchfit_sstd@fit$sigma)
mean(tgarchfit_norm@fit$sigma-tgarchfit_sstd@fit$sigma)
mean(gjrgarchfit_norm@fit$sigma-gjrgarchfit_sstd@fit$sigma)

# Model validation
par(mfrow = c(1,2))
acf(egarch33fit_norm@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(egarch33fit_norm@fit$z,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
dev.off
par(mfrow = c(1,2))
plot(egarch33fit_norm, which = 9) # QQ plot of the standardized residuals
plot(egarch33fit_norm, which = 8)
dev.off

#-------------------------------------------------------------------------------

# Determining eGARCH for sstd distribution using AIC
bic_egarch_sstd = matrix(0,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      egarch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                    model = "eGARCH"), 
                              mean.model = list(armaOrder=c(0,0),
                                                include.mean = TRUE),
                              distribution.model = "sstd")
      egarch_fit = ugarchfit(spec=egarch_spec, data=returns, 
                            solver = "gosolnp",
                            solver.control=list(trace = 1))
      bic_egarch_sstd[i,j] = infocriteria(egarch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

bic_egarch_sstd
which(min(bic_egarch_sstd) == bic_egarch_sstd)

par(mar=c(5, 5, 4, 5))
plot(bic_egarch_sstd, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, EGARCH, SSTD(0,1)",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off


# Model specification
egarchspec_sstd = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                                include.mean = TRUE), 
                              variance.model = list(garchOrder = c(1,1), 
                                                    model = "eGARCH"),
                              distribution.model = "sstd")
egarchfit_sstd = ugarchfit(data = returns, spec = egarchspec_sstd, solver = "gosolnp",
                            solver.control=list(trace = 1))
egarchfit_sstd
plot.ts(egarchfit_sstd@fit$z)
# Model validation
par(mfrow = c(1,2))
acf(egarchfit_sstd@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(egarchfit_sstd@fit$z,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
dev.off
par(mfrow = c(1,2))
plot(egarchfit_sstd, which = 9) # QQ plot of the standardized residuals
plot(egarchfit_sstd, which = 8)
dev.off

#-------------------------------------------------------------------------------

# Determining eGARCH for ged distribution using BIC
bic_egarch_ged = matrix(0,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      egarch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                    model = "eGARCH"), 
                              mean.model = list(armaOrder=c(0,0),
                                                include.mean = TRUE),
                              distribution.model = "ged")
      egarch_fit = ugarchfit(spec=egarch_spec, data=returns, 
                            solver = "gosolnp",
                            solver.control=list(trace = 1))
      bic_egarch_ged[i,j] = infocriteria(egarch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

bic_egarch_ged
which(min(bic_egarch_ged) == bic_egarch_ged)

par(mar=c(5, 5, 4, 5))
plot(bic_egarch_ged, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "BIC, GARCH models with skewed student's t-distribution",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
dev.off


# Model specification
egarchspec_ged = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                               include.mean = TRUE), 
                             variance.model = list(garchOrder = c(1,1), 
                                                   model = "eGARCH"),
                             distribution.model = "ged")
egarchfit_ged = ugarchfit(data = returns, spec = egarchspec_ged, solver = "gosolnp",
                           solver.control=list(trace = 1))
egarchfit_ged

# Model validation
par(mfrow = c(1,2))
acf(egarchfit_ged@fit$z, main = "ACF of the standardized residuals") # acf of the standardized residuals
acf(egarchfit_ged@fit$z,  main = "ACF of the squared standardized residuals") # acf of (standardized residuals)^2 
dev.off
par(mfrow = c(1,2))
plot(egarchfit_ged, which = 9) # QQ plot of the standardized residuals
plot(egarchfit_ged, which = 8)
dev.off

#------------------------------------------------------------------------------
# Simulations
# SSTD
sim = ugarchsim(fit = egarch43fit_sstd, n.sim = length(returns), n.start = 0, m.sim = 1000, rseed = NA)
kurt = kurtosis(sim@simulation$seriesSim)
mean(t(kurt))

# Normal
sim = ugarchsim(fit = egarch33fit_norm, n.sim = length(returns), n.start = 0, m.sim = 1000, rseed = NA)
kurt = kurtosis(sim@simulation$seriesSim)
mean(t(kurt))

# GED
sim = ugarchsim(fit = egarch33fit_ged, n.sim = length(returns), n.start = 0, m.sim = 1000, rseed = NA)
kurt = kurtosis(sim@simulation$seriesSim)
mean(t(kurt))


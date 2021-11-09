library(tseries)
library(forecast)
library(rugarch)
library(MLmetrics)
library(plot.matrix)

VIX=get.hist.quote(instrument = "^VIX",provider = "yahoo",
                   quote = "Close", end = "2021-10-15")
plot(VIX)

SP=get.hist.quote(instrument = "^GSPC",provider = "yahoo",
                  quote = "Close", end = "2021-10-15")
returns=diff(log(SP))
plot(returns)

hist(returns, breaks = 100, xlim = c(-0.1,0.1), col = "lightblue", 
     border = "blue3", freq = F, xlab = "S&P 500 returns", 
     main = "Histogram of S&P 500 returns")
curve(dnorm(x, mean = mean(returns), sd = sd(returns)), add = TRUE, col = "red")

# --------------------------------------------------------------

# Creating loop for fitting the best GARCH model according to AIC and BIC

# distribution assumed to be Gaussian --------------------------
aic_garch = matrix(0,5,5)
bic_garch = matrix(0,5,5)

for (i in 1:5) {
  for (j in 1:5) {
    garch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                  model = "sGARCH"), 
                            mean.model = list(armaOrder=c(0,0),
                                              include.mean = FALSE),
                            distribution.model = "norm")
    garch_fit = ugarchfit(spec=garch_spec, data=returns,
                          solver.control=list(trace = 1))
    aic_garch[i,j] = infocriteria(garch_fit)[1]
    bic_garch[i,j] = infocriteria(garch_fit)[2]
  }
}
aic_garch
bic_garch
which(min(aic_garch) == aic_garch)
which(min(bic_garch) == bic_garch)

par(mar=c(5, 6, 4, 8))
plot(aic_garch, col=heat.colors(n = 10, alpha = 0.9))
plot(bic_garch, col=heat.colors(n = 10, alpha = 0.9))
dev.off

# distribution assumed to be student T -------------------------
aic_garch_t = matrix(0,5,5)
bic_garch_t = matrix(0,5,5)

for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
    garch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                  model = "sGARCH"), 
                            mean.model = list(armaOrder=c(0,0),
                                              include.mean = FALSE),
                            distribution.model = "std")
    garch_fit = ugarchfit(spec=garch_spec, data=returns,
                          solver.control=list(trace = 1))
    aic_garch_t[i,j] = infocriteria(garch_fit)[1]
    bic_garch_t[i,j] = infocriteria(garch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}
aic_garch_t
bic_garch_t
which(min(aic_garch_t) == aic_garch_t)
which(min(bic_garch_t) == bic_garch_t)

aic_garch_t[aic_garch_t == 0] = NA
bic_garch_t[bic_garch_t == 0] = NA

par(mar=c(5, 6, 4, 8))
plot(aic_garch_t, col=heat.colors(n = 10, alpha = 0.9))
plot(bic_garch_t, col=heat.colors(n = 10, alpha = 0.9))
dev.off

# assuming Generalized Error Distribution ----------------------
aic_garch_ged = matrix(0,5,5)
bic_garch_ged = matrix(0,5,5)

for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      garch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                    model = "sGARCH"), 
                              mean.model = list(armaOrder=c(0,0),
                                                include.mean = FALSE),
                              distribution.model = "ged")
      garch_fit = ugarchfit(spec=garch_spec, data=returns,
                            solver.control=list(trace = 1))
      aic_garch_ged[i,j] = infocriteria(garch_fit)[1]
      bic_garch_ged[i,j] = infocriteria(garch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}
aic_garch_ged
bic_garch_ged
which(min(aic_garch_ged) == aic_garch_ged)
which(min(bic_garch_ged) == bic_garch_ged)

par(mar=c(5, 6, 4, 8))
plot(aic_garch_ged, col=heat.colors(n = 10, alpha = 0.9))
plot(bic_garch_ged, col=heat.colors(n = 10, alpha = 0.9))
dev.off


# --------------------------------------------------------------

# fitting GARCH(1,1) model to returns
garch11spec = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                           include.mean = FALSE), 
                       variance.model = list(garchOrder = c(1,1), 
                                             model = "sGARCH"),
                       distribution.model = "norm")
garch11fit = ugarchfit(data = returns, spec = garch11spec)
vol_ret_garch11 = sigma(garch11fit)
plot(vol_ret_garch11)

# fitting GARCH(2,2) model to returns
garch22spec = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                           include.mean = FALSE), 
                         variance.model = list(garchOrder = c(2,2), 
                                               model = "sGARCH"),
                         distribution.model = "norm")
garch22fit = ugarchfit(data = returns, spec = garch22spec)
vol_ret_garch22 = sigma(garch22fit)
plot(vol_ret_garch22)

# model validation
garch11fit
res = residuals(garch11fit)
res_standardized = residuals(garch11fit, standardize = T)
plot(res)
plot(res_standardized) # should approach normal distribution
acf(x = res_standardized^2) # should not be autocorrelated
hist(res, breaks = 100, xlim = c(-0.1,0.1), col = "lightblue", 
     border = "blue3", freq = T)
qqnorm(res_standardized)
qqline(res_standardized)
shapiro.test(res[1:500])

plot(garch11fit, which = 3) # Conditional SD
plot(garch11fit, which = 8) # Density of Stand. Res.
plot(garch11fit, which = 9) # QQ-PLot of Stand. Res.
plot(garch11fit, which = 10) # ACF of Stand. Res.
plot(garch11fit, which = 11) # ACF of Squared Stand. Res.

plot(garch22fit, which = 3) # Conditional SD
plot(garch22fit, which = 8) # Density of Stand. Res.
plot(garch22fit, which = 9) # QQ-PLot of Stand. Res.
plot(garch22fit, which = 10) # ACF of Stand. Res.
plot(garch22fit, which = 11) # ACF of Squared Stand. Res.


# fitting GARCH(1,1) model to returns using EGARCH
egarchspec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                        variance.model = list(model = "eGARCH"),
                        distribution.model = "norm")
egarchfit = ugarchfit(data = returns, spec = egarchspec)
e.vol.ret = sigma(egarchfit)
plot(e.vol.ret)
lines(vol.ret, col = "red")

#modelling CBOE VIX using returns
mod1=lm((VIX[-1]) ~ vol.ret)
plot(VIX[-1])
lines(fitted(mod1), col = "red")
MSE.g.ret = MSE(y_pred = fitted(mod1), y_true = VIX)
MSE.g.ret # 14.03579

mod2=lm((VIX[-1]) ~ e.vol.ret)
plot(VIX[-1])
lines(fitted(mod2), col = "red")
MSE.eg.ret = MSE(y_pred = fitted(mod2), y_true = VIX)
MSE.eg.ret # 13.3777



#modelling VIX using GARCH
garchspec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                       variance.model = list(model = "sGARCH"),
                       distribution.model = "norm")
garchfit = ugarchfit(data = VIX, spec = garchspec)
vol.vix = sigma(garchfit)
plot(vol.vix)

#modelling VIX using EGARCH (does not work)
egarchspec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                        variance.model = list(model = "eGARCH"),
                        distribution.model = "norm")
egarchfit = ugarchfit(data = VIX, spec = egarchspec, 
                      solver = "gosolnp", 
                      solver.control = list(tol = 10))
e.vol.vix = sigma(egarchfit)
plot(e.vol.vix)


#modelling CBOE VIX using VIX
mod3=lm((VIX) ~ vol.vix)
plot(VIX)
lines(fitted(mod3), col = "red")
MSE.g.vix = MSE(y_pred = fitted(mod3), y_true = VIX)
MSE.g.vix # 9.448859


#modelling CBOE VIX using returns and VIX
mod4=lm((VIX[-1]) ~ vol.ret + vol.vix[-1])
plot(VIX[-1])
lines(fitted(mod4), col = "red")
MSE.g.retvix = MSE(y_pred = fitted(mod4), y_true = VIX)
MSE.g.retvix # 8.125895



# --------------------------------------------------------------
# rolling forecast

train = returns[1:7656]
garchspec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                       variance.model = list(model = "sGARCH"),
                       distribution.model = "norm")
garchroll=ugarchroll(spec = garchspec, data = train, n.ahead = 1, 
                     n.start = 3000, refit.every = 300, 
                     refit.window = "recursive")
pred = as.data.frame(garchroll)
plot.ts(pred$Sigma)
plot(garchroll)

# forecast
garchfitfore = ugarchfit(data = returns, spec = garchspec, out.sample = 1000)
foremod = ugarchforecast(fitORspec = garchfitfore, data = returns, 
                         n.ahead = 100, out.sample = 1000, n.roll = 10)
plot(foremod)

library(tseries)
library(quantmod)
library(forecast)
library(rugarch)
library(MLmetrics)
library(plot.matrix)
library(tictoc)
library(furrr)
library(tidyverse)

#VIX = get.hist.quote(instrument = "^VIX",provider = "yahoo",
#                   quote = "Close", end = "2021-10-15")
getSymbols("^VIX", src = "yahoo", from = "1990-01-01", to = "2020-12-31")
VIX = VIX$VIX.Close[-1]
plot(VIX)

#SP = get.hist.quote(instrument = "^GSPC",provider = "yahoo",
#                  quote = "Close", end = "2021-10-15")
#returns=diff(log(SP))
getSymbols("^GSPC", src = "yahoo", from = "1990-01-01", to = "2020-12-31")
par(mfrow=c(2,1))
plot(GSPC$GSPC.Close, main = "S&P 500 index")
returns = diff(log(GSPC$GSPC.Close))[-1]
plot(returns, main = "Returns of the S&P 500")

hist(returns, breaks = 100, xlim = c(-0.1,0.1), col = "lightblue", 
     border = "blue3", freq = F, xlab = "S&P 500 returns", 
     main = "Histogram of S&P 500 returns")
curve(dnorm(x, mean = mean(returns), sd = sd(returns)), 
      add = TRUE, col = "red")

# --------------------------------------------------------------

# Creating loop in loop for finding GARCH and ARMA orders 
# according to AIC and BIC, where distribution is assumed to be 
# Gaussian -----------------------------------------------------

aic_arma_order = matrix(0, nrow = 12, ncol = 12)
bic_arma_order = matrix(0, nrow = 12, ncol = 12)
plan(multisession, workers = 8)
set.seed(91)
tic()
for (k in 0:2) {
  for (l in 0:2) {
    for (i in 1:4) {
      for (j in 1:4) {
        tryCatch({
          garch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                        model = "sGARCH"), 
                                  mean.model = list(armaOrder=c(k,l),
                                                    include.mean = FALSE),
                                  distribution.model = "norm")
          garch_fit = ugarchfit(spec=garch_spec, data=returns,
                                solver.control = list(trace = 1, n.restarts = 1),
                                solver = "gosolnp",
                                parallel = T, cores = 8)
          aic_arma_order[(k*4+i),(l*4+j)] = infocriteria(garch_fit)[1]
          bic_arma_order[(k*4+i),(l*4+j)] = infocriteria(garch_fit)[2]
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
  }
}
toc() # 700 sec

aic_arma_order
bic_arma_order

aic_arma_order[aic_arma_order == 0] = NA
bic_arma_order[bic_arma_order == 0] = NA

par(mar=c(5, 4, 4, 5))
plot(aic_arma_order, col=heat.colors(n = 10, alpha = 0.9), digits = 4, 
     cex = 0.8)
plot(bic_arma_order, col=heat.colors(n = 10, alpha = 0.9), digits = 4, 
     cex = 0.8)
which.min(aic_arma_order)
which.min(bic_arma_order)


# --------------------------------------------------------------

# Creating loop for finding GARCH orders according to AIC and BIC

# distribution assumed to be Gaussian --------------------------
aic_garch_norm = matrix(0,5,5)
bic_garch_norm = matrix(0,5,5)

tic()
for (i in 1:5) {
  for (j in 1:5) {
    tryCatch({
      garch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                    model = "sGARCH"), 
                              mean.model = list(armaOrder=c(0,0),
                                                include.mean = FALSE),
                              distribution.model = "norm")
      garch_fit = ugarchfit(spec=garch_spec, data=returns, 
                            solver = "gosolnp",
                            solver.control=list(trace = 1))
      aic_garch_norm[i,j] = infocriteria(garch_fit)[1]
      bic_garch_norm[i,j] = infocriteria(garch_fit)[2]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}
toc()
aic_garch_norm
bic_garch_norm
which(min(aic_garch_norm) == aic_garch_norm)
which(min(bic_garch_norm) == bic_garch_norm)

par(mar=c(5, 5, 4, 5))
plot(aic_garch_norm, col=heat.colors(n = 10, alpha = 0.9), digits = 4,
     main = "AIC, GARCH models with Gaussian Distribution",
     xlab = expression(paste("Column, ", italic("q"), " order")), 
     ylab = expression(paste("Row, ", italic("p"), " order")))
plot(bic_garch_norm, col=heat.colors(n = 10, alpha = 0.9), digits = 4, 
     main = "BIC, GARCH models with Gaussian Distribution")
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
                          solver = "gosolnp",
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

par(mar=c(5, 6, 4, 8))
plot(aic_garch_t, col=heat.colors(n = 10, alpha = 0.9), 
     main = "AIC, GARCH models with Student T Distribution")
plot(bic_garch_t, col=heat.colors(n = 10, alpha = 0.9), 
     main = "BIC, GARCH models with Student T Distribution")
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
plot(aic_garch_ged, col=heat.colors(n = 10, alpha = 0.9), 
     main = "AIC, GARCH models with GGD")
plot(bic_garch_ged, col=heat.colors(n = 10, alpha = 0.9), 
     main = "BIC, GARCH models with GGD")
dev.off


# --------------------------------------------------------------

# checking which distribution seems to be most accurate --------
min(aic_garch_norm)
min(aic_garch_t)
min(aic_garch_ged)
min(bic_garch_norm)
min(bic_garch_t)
min(bic_garch_ged)

# fitting GARCH models -----------------------------------------

# fitting GARCH(1,1) model to returns
garch11spec_norm = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                           include.mean = FALSE), 
                       variance.model = list(garchOrder = c(1,1), 
                                             model = "sGARCH"),
                       distribution.model = "norm")
garch11fit_norm = ugarchfit(data = returns, spec = garch11spec_norm)
plot(sigma(garch11fit_norm))
plot(garch11fit_norm)
garch11fit_norm
garch11fit_norm@fit$z*garch11fit_norm@fit$sigma - returns
plot(garch11fit_norm@fit$z)


# fitting GARCH(2,2) model to returns
garch22spec_norm = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                           include.mean = FALSE), 
                         variance.model = list(garchOrder = c(2,2), 
                                               model = "sGARCH"),
                         distribution.model = "norm")
garch22fit_norm = ugarchfit(data = returns, spec = garch22spec_norm)
plot(sigma(garch22fit_norm))

# fitting GARCH(2,2) model to returns assuming Student T Distribution
garch22spec_std = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                           include.mean = FALSE), 
                         variance.model = list(garchOrder = c(2,2), 
                                               model = "sGARCH"),
                         distribution.model = "std")
garch22fit_std = ugarchfit(data = returns, spec = garch22spec_std, 
                           solver="gosolnp")
plot(sigma(garch22fit_std))
plot.ts(garch22fit_std@fit$z)

# fitting GARCH(2,2) model to returns using Generalized Gaussian Distribution
garch22spec_ged = ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                               include.mean = FALSE), 
                             variance.model = list(garchOrder = c(2,2), 
                                                   model = "sGARCH"),
                             distribution.model = "ged")
garch22fit_ged = ugarchfit(data = returns, spec = garch22spec_ged)
plot(sigma(garch22fit_ged))
r

# model validation ---------------------------------------------

garch11fit_norm
res = residuals(garch11fit_norm)
res_standardized = residuals(garch11fit_norm, standardize = T)
plot(res)
plot(res_standardized) # should approach normal distribution
acf(x = res_standardized^2) # should not be autocorrelated
hist(res, breaks = 100, xlim = c(-0.1,0.1), col = "lightblue", 
     border = "blue3", freq = T)
qqnorm(res_standardized)
qqline(res_standardized)
shapiro.test(res[1:500])

plot(garch11fit_norm, which = 3) # Conditional SD
plot(garch11fit_norm, which = 8) # Density of Stand. Res.
plot(garch11fit_norm, which = 9) # QQ-PLot of Stand. Res.
plot(garch11fit_norm, which = 10) # ACF of Stand. Res.
plot(garch11fit_norm, which = 11) # ACF of Squared Stand. Res.

plot(garch22fit_norm, which = 3) # Conditional SD
plot(garch22fit_norm, which = 8) # Density of Stand. Res.
plot(garch22fit_norm, which = 9) # QQ-PLot of Stand. Res.
plot(garch22fit_norm, which = 10) # ACF of Stand. Res.
plot(garch22fit_norm, which = 11) # ACF of Squared Stand. Res.

plot(garch22fit_std, which = 3) # Conditional SD
plot(garch22fit_std, which = 8) # Density of Stand. Res.
plot(garch22fit_std, which = 9) # QQ-PLot of Stand. Res.
plot(garch22fit_std, which = 10) # ACF of Stand. Res.
plot(garch22fit_std, which = 11) # ACF of Squared Stand. Res.

plot(garch22fit_ged, which = 3) # Conditional SD
plot(garch22fit_ged, which = 8) # Density of Stand. Res.
plot(garch22fit_ged, which = 9) # QQ-PLot of Stand. Res.
plot(garch22fit_ged, which = 10) # ACF of Stand. Res.
plot(garch22fit_ged, which = 11) # ACF of Squared Stand. Res.


# --------------------------------------------------------------

# LM model -----------------------------------------------------
# fitting EGARCH(1,1) model to returns using EGARCH
egarchspec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                        variance.model = list(model = "eGARCH"),
                        distribution.model = "norm")
egarchfit = ugarchfit(data = returns, spec = egarchspec)
e.vol.ret = sigma(egarchfit)
plot(e.vol.ret)
lines(vol.ret, col = "red")

# modelling CBOE VIX using returns
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



# fitting GARCH(1,1) model to VIX
garchspec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                       variance.model = list(model = "sGARCH"),
                       distribution.model = "norm")
garchfit = ugarchfit(data = VIX, spec = garchspec)
vol.vix = sigma(garchfit)
plot(vol.vix)

# fitting EGARCH(1,1) model to VIX (does not work)
egarchspec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                        variance.model = list(model = "eGARCH"),
                        distribution.model = "norm")
egarchfit = ugarchfit(data = VIX, spec = egarchspec, 
                      solver = "gosolnp", 
                      solver.control = list(tol = 10))
e.vol.vix = sigma(egarchfit)
plot(e.vol.vix)


# modelling CBOE VIX using VIX
mod3=lm((VIX) ~ vol.vix)
plot(VIX)
lines(fitted(mod3), col = "red")
MSE.g.vix = MSE(y_pred = fitted(mod3), y_true = VIX)
MSE.g.vix # 9.448859


# modelling CBOE VIX using returns and VIX
mod4=lm((VIX[-1]) ~ vol.ret + vol.vix[-1])
plot(VIX[-1])
lines(fitted(mod4), col = "red")
MSE.g.retvix = MSE(y_pred = fitted(mod4), y_true = VIX)
MSE.g.retvix # 8.125895



# --------------------------------------------------------------
# rolling forecast

train = returns["/2020-12-31"]
test = returns["2021-01-01/2021-10-31"]
garchspec = ugarchspec(variance.model = list(garchOrder=c(1,1), 
                                             model = "eGARCH"), 
                       mean.model = list(armaOrder=c(0,0),
                                         include.mean = T),
                       distribution.model = "ged")
garchroll=ugarchroll(spec = garchspec, data = returns, n.ahead = 1, 
                     n.start = length(train), forecast.length = length(test), 
                     refit.every = 30, refit.window = "recursive",
                     solver = "solnp")
pred = as.data.frame(garchroll)
my_dates <- strptime(rownames(pred), format="%Y-%m-%d")
ggplot(data = pred) + labs(x = "Date") +
  geom_line(aes(x = as.Date(rownames(pred)), y = Sigma), col = "steelblue")

plot(garchroll, which = 2)

# forecast roll and static
garchfitfore = ugarchfit(data = returns, spec = garchspec, 
                         out.sample = length(test))
foremod = ugarchforecast(fitORspec = garchfitfore, data = returns, 
                         n.ahead = 10, n.roll = length(test), 
                         external.forecasts = list())
plot(foremod, which = 4)

plot.ts(foremod@model$modeldata$sigma)
plot.ts(foremod@forecast$sigmaFor[1,])

lm(VIX["/2020-12-31"] ~ sigma(garchfitfore))
plot.ts(VIX["2021-01-01/2021-10-31"], ylim = c(10,40))
lines(5.86+1379.83*foremod@forecast$sigmaFor[1,], col = "red")



#-------------------------------------------------------------------------------
#----------------- The actual procedure (sGARCH) ----------------------------------------
#-------------------------------------------------------------------------------

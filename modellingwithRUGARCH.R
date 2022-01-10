library(tseries)
library(quantmod)
library(forecast)
library(rugarch)
library(MLmetrics)
library(plot.matrix)
library(tictoc)
library(furrr)
library(tidyverse)
library(patchwork)

#VIX = get.hist.quote(instrument = "^VIX",provider = "yahoo",
#                   quote = "Close", end = "2021-10-15")
getSymbols("^VIX", src = "yahoo", from = "1990-01-01", to = "2020-12-31")
VIX = VIX$VIX.Close[-1]
plot(VIX, main = "CBOE VIX", format.labels="%b\n%Y")

#SP = get.hist.quote(instrument = "^GSPC",provider = "yahoo",
#                  quote = "Close", end = "2021-10-15")
#returns=diff(log(SP))
getSymbols("^GSPC", src = "yahoo", from = "1990-01-01", to = "2021-11-30")
par(mfrow=c(2,1))
plot(GSPC$GSPC.Close, main = "S&P 500 index", format.labels="%b\n%Y")
returns = diff(log(GSPC$GSPC.Close))[-1]
(autoplot(GSPC$GSPC.Close, main = "S&P 500 index") + xlab("") + ylab("Price")) /
  (autoplot(returns, main = "Returns of the S&P 500 index") +  xlab("") + ylab("Return"))




hist(returns, breaks = 100, xlim = c(-0.1,0.1), col = "lightblue", 
     border = "blue3", freq = F, xlab = "S&P 500 returns", 
     main = "Histogram of S&P 500 returns")
curve(dnorm(x, mean = mean(returns), sd = sd(returns)), 
      add = TRUE, col = "red")

# --------------------------------------------------------------

# Creating loop in loop for finding GARCH and ARMA orders 
# according to AIC and BIC, where distribution is assumed to be 
# Gaussian -----------------------------------------------------

bic_arma_order = matrix(0, nrow = 16, ncol = 16)
plan(multisession, workers = 8)
set.seed(91)
tic()
for (k in 0:3) {
  for (l in 0:3) {
    for (i in 1:4) {
      for (j in 1:4) {
        tryCatch({
          garch_spec = ugarchspec(variance.model = list(garchOrder=c(i,j), 
                                                        model = "sGARCH"), 
                                  mean.model = list(armaOrder=c(k,l),
                                                    include.mean = FALSE),
                                  distribution.model = "sstd")
          garch_fit = ugarchfit(spec=garch_spec, data=returns,
                                solver.control = list(trace = 1, n.restarts = 1),
                                solver = "gosolnp",
                                parallel = T, cores = 8)
          # aic_arma_order[(k*4+i),(l*4+j)] = infocriteria(garch_fit)[1]
          bic_arma_order[(k*4+i),(l*4+j)] = infocriteria(garch_fit)[2]
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
  }
}
toc() # 700 sec

bic_arma_order

bic_arma_order[bic_arma_order == 0] = NA

colnames(bic_arma_order) = rep(" ",16)
rownames(bic_arma_order) = rep(" ",16)

xtick = c(expression(paste(q[A],"=0, ",q[G],"=1")),
          expression(paste(q[A],"=0, ",q[G],"=2")),
          expression(paste(q[A],"=0, ",q[G],"=3")),
          expression(paste(q[A],"=0, ",q[G],"=4")),
          expression(paste(q[A],"=1, ",q[G],"=1")),
          expression(paste(q[A],"=1, ",q[G],"=2")),
          expression(paste(q[A],"=1, ",q[G],"=3")),
          expression(paste(q[A],"=1, ",q[G],"=4")),
          expression(paste(q[A],"=2, ",q[G],"=1")),
          expression(paste(q[A],"=2, ",q[G],"=2")),
          expression(paste(q[A],"=2, ",q[G],"=3")),
          expression(paste(q[A],"=2, ",q[G],"=4")),
          expression(paste(q[A],"=3, ",q[G],"=1")),
          expression(paste(q[A],"=3, ",q[G],"=2")),
          expression(paste(q[A],"=3, ",q[G],"=3")),
          expression(paste(q[A],"=3, ",q[G],"=4"))
          )
ytick = c(expression(paste(p[A],"=0, ",p[G],"=1")),
          expression(paste(p[A],"=0, ",p[G],"=2")),
          expression(paste(p[A],"=0, ",p[G],"=3")),
          expression(paste(p[A],"=0, ",p[G],"=4")),
          expression(paste(p[A],"=1, ",p[G],"=1")),
          expression(paste(p[A],"=1, ",p[G],"=2")),
          expression(paste(p[A],"=1, ",p[G],"=3")),
          expression(paste(p[A],"=1, ",p[G],"=4")),
          expression(paste(p[A],"=2, ",p[G],"=1")),
          expression(paste(p[A],"=2, ",p[G],"=2")),
          expression(paste(p[A],"=2, ",p[G],"=3")),
          expression(paste(p[A],"=2, ",p[G],"=4")),
          expression(paste(p[A],"=3, ",p[G],"=1")),
          expression(paste(p[A],"=3, ",p[G],"=2")),
          expression(paste(p[A],"=3, ",p[G],"=3")),
          expression(paste(p[A],"=3, ",p[G],"=4"))
          )
par(mar=c(2, 4, 4, 4.5))
plot(bic_arma_order, col=heat.colors(n = 10, alpha = 0.9), digits = 4, 
     cex = 0.8, main = "", xlab='', ylab='', axis.col=3)
axis(side = 3, at = seq(1,16,1), labels = xtick, las = 3, cex.axis = 0.9)
axis(side = 2, at = seq(1,16,1), labels = ytick, las = 2, cex.axis = 0.9)

which.min(bic_arma_order)
which(is.na(bic_arma_order))

ytick = c("P=3,p=4","P=3,p=3","P=3,p=2","P=3,p=1",
          "P=2,p=4","P=2,p=3","P=2,p=2","P=2,p=1",
          "P=1,p=4","P=1,p=3","P=1,p=2","P=1,p=1",
          "P=0,p=4","P=0,p=3","P=0,p=2","P=0,p=1")
xtick = c("Q=0,q=1","Q=0,q=2","Q=0,q=3","Q=0,q=4",
          "Q=1,q=1","Q=1,q=2","Q=1,q=3","Q=1,q=4",
          "Q=2,q=1","Q=2,q=2","Q=2,q=3","Q=2,q=4",
          "Q=3,q=1","Q=3,q=2","Q=3,q=3","Q=3,q=4")

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

getSymbols("DTB3", src = "FRED", from = "1990-01-01", to = "2021-11-30")

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
#----------------- Misc plots---------- ----------------------------------------
#-------------------------------------------------------------------------------
vix1 = data.frame(VIX)
vix = vix1 %>% transmute(garchfit_sstd@fit$sigma)
sigma1 = as.xts(vix)

vix2 = data.frame(VIX)
vix = vix2 %>% transmute(egarchfit_sstd@fit$sigma)
sigma2 = as.xts(vix)

vix3 = data.frame(VIX)
vix = vix3 %>% transmute(tgarchfit_sstd@fit$sigma)
sigma3 = as.xts(vix)

vix4 = data.frame(VIX)
vix = vix4 %>% transmute(gjrgarchfit_sstd@fit$sigma)
sigma4 = as.xts(vix)

# 2x2 plot
((autoplot(sigma1, main = "GARCH, SSTD(0,1)")+ ylim(0.003,0.075) + ylab(expression(sigma[t])) + xlab("")) + 
  (autoplot(sigma2, main = "EGARCH, SSTD(0,1)")+ ylim(0.003,0.075) + ylab(expression(sigma[t])) + xlab(""))) /
  ((autoplot(sigma3, main = "TGARCH, SSTD(0,1)")+ ylim(0.003,0.075) + ylab(expression(sigma[t])) + xlab("")) + 
   (autoplot(sigma4, main = "GJR-GARCH, SSTD(0,1)")+ ylim(0.003,0.075) + ylab(expression(sigma[t])) + xlab("")))


((autoplot(sigma1, main = "GARCH, SSTD(0,1)")+ ylim(0.003,0.075) + ylab(expression(sigma[t])) + xlab(""))) /
((autoplot(sigma4, main = "GJR-GARCH, SSTD(0,1)")+ ylim(0.003,0.075) + ylab(expression(sigma[t])) + xlab("")))


# 4x1 plot
(autoplot(sigma1, main = "GARCH, SSTD(0,1)")+ ylim(0.003,0.075) + ylab(expression(sigma[t])) + xlab("")) / 
    (autoplot(sigma2, main = "EGARCH, SSTD(0,1)")+ ylim(0.003,0.075) + ylab(expression(sigma[t])) + xlab("")) /
  (autoplot(sigma3, main = "TGARCH, SSTD(0,1)")+ ylim(0.003,0.075) + ylab(expression(sigma[t])) + xlab("")) /
     (autoplot(sigma4, main = "GJR-GARCH, SSTD(0,1)")+ ylim(0.003,0.075) + ylab(expression(sigma[t])) + xlab(""))

autoplot(VIX, main = "CBOE VIX") + ylab("") + xlab("")

autoplot(garchfit_sstd@fit$sigma)

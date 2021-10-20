library(tseries)
library(forecast)

VIX=get.hist.quote(instrument = "^VIX",provider = "yahoo",quote = "Close", 
                   end = "2021-10-15")
VIXret=diff(log(VIX))
plot(VIX)

SP=get.hist.quote(instrument = "^GSPC",provider = "yahoo",quote = "Close", 
                  end = "2021-10-15")
returns=diff(log(SP))
plot(returns)

library(rugarch)
library(MLmetrics)

#modelling returns using GARCH
garchspec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                       variance.model = list(model = "sGARCH"),
                       distribution.model = "norm")
garchfit = ugarchfit(data = returns, spec = garchspec)
vol.ret = sigma(garchfit)
plot(vol.ret)

#modelling returns using EGARCH
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
MSE.g.ret = MSE(y_pred = fitted(mod1), y_true = VIX) ; MSE.g.ret # 14.03579

mod2=lm((VIX[-1]) ~ e.vol.ret)
plot(VIX[-1])
lines(fitted(mod2), col = "red")
MSE.eg.ret = MSE(y_pred = fitted(mod2), y_true = VIX) ; MSE.eg.ret # 13.3777



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
                      solver = "gosolnp", solver.control = list(tol = 100))
e.vol.vix = sigma(egarchfit)
plot(e.vol.vix)
lines(vol.vix, col = "red")


#modelling CBOE VIX using VIX
mod3=lm((VIX) ~ vol.vix)
plot(VIX)
lines(fitted(mod3), col = "red")
MSE.g.vix = MSE(y_pred = fitted(mod3), y_true = VIX) ; MSE.g.vix # 9.448859


#modelling CBOE VIX using returns and VIX
mod4=lm((VIX[-1]) ~ vol.ret + vol.vix[-1])
plot(VIX[-1])
lines(fitted(mod4), col = "red")
MSE.g.retvix = MSE(y_pred = fitted(mod4), y_true = VIX) ; MSE.g.retvix # 8.125895


#__________________________________________________________
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
foremod = ugarchforecast(fitORspec = garchfitfore, data = returns, n.ahead = 100, 
                         out.sample = 1000, n.roll = 10)
plot(foremod)

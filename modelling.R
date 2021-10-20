library(tseries)

VIX=get.hist.quote(instrument = "^VIX",provider = "yahoo",quote = "Close")
VIX=diff(VIX)
plot(diff(VIX))

SP=get.hist.quote(instrument = "^GSPC",provider = "yahoo",quote = "Close")
returns=diff(log(SP))
plot(returns)

library(fGarch)

mod=garchFit(x = returns, order = c(1,1))
summary(mod)
plot(mod)

mod1=garchFit(x = VIX, order = c(1,1))
summary(mod1)
plot(mod1)

plot(mod1@h.t)

plot.ts(mod@sigma.t)
plot(VIX)

library(forecast)

pred=forecast(mod1@h.t, h = 10)
autoplot(pred, xlim = c(1800,2000))


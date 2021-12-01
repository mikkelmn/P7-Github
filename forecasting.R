# packages ----------------
library(tidyverse)
library(quantmod)
library(rugarch)
library(xts)
library(timetk)
library(forecast)

# data -------------------
getSymbols("^VIX", src = "yahoo", from = "1990-01-01", to = "2021-10-31")
VIX = VIX$VIX.Close[-1]
getSymbols("^GSPC", src = "yahoo", from = "1990-01-01", to = "2021-10-31")
returns = diff(log(GSPC$GSPC.Close))[-1]
SP500_Volume = GSPC$GSPC.Volume[-1]
getSymbols("DGS3MO", src = "FRED")
TBill3m = na.approx(DGS3MO["1990-01-01/2021-10-31"]) #no [-1], it is already na

# rolling forecast ---------------

train = returns["/2020-12-31"]
test = returns["2021-01-01/2021-10-31"]
garchspec = ugarchspec(variance.model = list(garchOrder=c(1,1), 
                                             model = "sGARCH"), 
                       mean.model = list(armaOrder=c(0,0),
                                         include.mean = T),
                       distribution.model = "norm")
garchroll=ugarchroll(spec = garchspec, data = returns, n.ahead = 1, 
                     n.start = length(train), forecast.length = length(test), 
                     refit.every = 30, refit.window = "recursive",
                     solver = "solnp")
pred = as.data.frame(garchroll)
my_dates <- strptime(rownames(pred), format="%Y-%m-%d")
ggplot(data = pred) + labs(x = "Date") +
  geom_line(aes(x = as.Date(rownames(pred)), y = Sigma), col = "steelblue")

plot(garchroll, which = 2)

# rolling forecast of the S&P500 returns ----------------

test10y = returns["2011-11-01/2021-10-31"]
test21 = returns["2021-01-01/2021-10-31"]
specs = ugarchspec(variance.model = list(garchOrder=c(1,1), model = "eGARCH"),
                       mean.model = list(armaOrder=c(0,0), include.mean = T),
                       distribution.model = "sstd")
garchfitfore = ugarchfit(data = returns, spec = specs, 
                         out.sample = length(test21))
foremod = ugarchforecast(fitORspec = garchfitfore, n.ahead = length(test21), 
                         n.roll = length(test21))
predval = foremod@forecast$sigmaFor[1,]

plot(foremod, which = "all")
plot.ts(foremod@model$modeldata$sigma)
plot.ts(predval)

# rolling forecast of VIX returns ----------------

VIXreturns = diff(log(VIX))
plot(VIXreturns)
specs = ugarchspec(variance.model = list(garchOrder=c(1,1), model = "sGARCH"), 
                   mean.model = list(armaOrder=c(0,0), include.mean = T),
                   distribution.model = "norm")
fit_VIXreturns = ugarchfit(data = VIXreturns[-1], spec = specs, 
                         out.sample = length(test21))
fc_VIXreturns = ugarchforecast(fitORspec = fit_VIXreturns, data = VIXreturns[-1], 
                         n.ahead = 10, n.roll = length(test21))
predval_VIXreturns = fc_VIXreturns@forecast$sigmaFor[1,]
plot(fc_VIXreturns, which = 4)
plot.ts(sigma(fit_VIXreturns))

# forecast VIX using volatility of forecasted returns ---------------
lm(VIX["/2011-10-31"] ~ sigma(garchfitfore))
plot.ts(VIX["2011-11-01/2021-10-31"], ylim = c(5,90))
lines(6.212 + 1365.612 * foremod@forecast$sigmaFor[1,], col = "steelblue")

tib = tibble(fortify.zoo(VIX)) %>% 
  rename(VIX = VIX.Close)
df = tib %>% filter(Index >= "2011-11-01") %>% 
    mutate(RollForecast = 6.212 + 1365.612 * predval[-1])

ggplot(data = df) + geom_line(aes(x = Index, y = VIX, col = "Actual VIX")) +
  geom_line(aes(x = Index, y = RollForecast, col = "Rolling Forecast")) +
  labs(x = "Date", y = "Value")
(MSE1 = mean((df$VIX - df$RollForecast)^2))

# determining AR(p) of VIX -----------------

auto.arima(VIX["/2021-01-01"], max.q = 0)
acf(VIX["/2011-10-31"])
pacf(VIX["/2011-10-31"])
acf(diff(VIX["/2011-10-31"]), na.action = na.pass)
pacf(diff(VIX["/2011-10-31"]), na.action = na.pass)

# forecast VIX using volatility of returns and lags of VIX -----------

tib = tibble(fortify.zoo(VIX)) %>% 
  left_join(x = ., y = fortify.zoo(sigma(garchfitfore)), by = "Index") %>% 
  left_join(x = ., y = fortify.zoo(SP500_Volume), by = "Index") %>% 
  left_join(x = ., y = fortify.zoo(TBill3m), by = "Index") %>%
  dplyr::rename(
    VIX = VIX.Close,
    SigmaReturns = "sigma(garchfitfore)",
    TradingVolume = GSPC.Volume,
    RiskFreeRate = DGS3MO
    ) %>% 
  timetk::tk_augment_lags(.value = c(TradingVolume, RiskFreeRate), .lags = 1) %>% 
#  timetk::tk_augment_lags(.value = VIX, .lags = c(1:5), .names = "auto") %>% 
  dplyr::filter(Index > "1990-01-10")

mod = lm(VIX ~ SigmaReturns + TradingVolume_lag1 + RiskFreeRate_lag1, 
#         + VIX_lag1 + VIX_lag2 + VIX_lag3 + VIX_lag4 + VIX_lag5, 
         data = tib %>% filter(Index <= "2021-01-01"))
summary(mod)

df = tib %>% filter(Index >= "2021-01-01") %>% 
  mutate(RollForecast = coef(mod)[1] + coef(mod)[2] * predval[-1] 
         + coef(mod)[3] * TradingVolume_lag1 + coef(mod)[4] * RiskFreeRate_lag1)
#         + coef(mod)[5] * VIX_lag1 + coef(mod)[6] * VIX_lag2 
#         + coef(mod)[7] * VIX_lag3 + coef(mod)[8] * VIX_lag4
#         + coef(mod)[9] * VIX_lag5)

ggplot(data = df) + geom_line(aes(x = Index, y = VIX, col = "Actual VIX")) +
  geom_line(aes(x = Index, y = RollForecast, col = "Rolling Forecast")) +
  labs(x = "Date", y = "Value")
MSE_OutOfSample[4] = mean((df$VIX - df$RollForecast)^2)

# in sample
tib = tib %>% filter(Index <= "2021-01-01")
ggplot(data = tib) + geom_line(aes(x = Index, y = VIX, col = "Actual VIX")) +
  geom_line(aes(x = Index, y = fitted(mod), col = "Fitted Values")) +
  labs(x = "Date", y = "Value")
MSE_InSample[4] = mean((tib$VIX - fitted(mod))^2)

ModelName = c("GARCH normal", "GARCH skewed T", 
              "E-GARCH normal", "E-GARCH skewed T", 
              "T-GARCH normal", "T-GARCH skewed T", 
              "GJR-GARCH normal", "GJR-GARCH skewed T", 
              "Without Sigma", "Only Lags")

options(pillar.sigfig = 5) %>% 
  tibble(ModelName,MSE_InSample,MSE_OutOfSample,
         MSE_InSample_withlags,MSE_OutOfSample_withlags)




mean((df$VIX - df$VIX_lag1)^2)



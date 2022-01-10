# data -------------------
getSymbols("^VIX", src = "yahoo", from = "1990-01-01", to = "2021-11-30")
VIX = VIX$VIX.Close[-1]
getSymbols("^GSPC", src = "yahoo", from = "1990-01-01", to = "2021-11-30")
returns = diff(log(GSPC$GSPC.Close))[-1]
SP500_Volume = GSPC$GSPC.Volume[-1]
getSymbols("DTB3", src = "FRED")
TBill3m = na.approx(DTB3["1990-01-01/2021-11-30"])

# rolling forecast of the S&P500 returns ----------------

test10y = returns["2011-11-01/2021-11-30"]
test21 = returns["2021-01-01/2021-11-30"]
specs = ugarchspec(variance.model = list(garchOrder=c(1,1), model = "sGARCH"),
                   mean.model = list(armaOrder=c(0,0), include.mean = T),
                   distribution.model = "sstd")
fit = ugarchfit(data = returns, spec = specs, out.sample = length(test21))
foremod = ugarchforecast(fitORspec = fit, n.ahead = length(test21), 
                         n.roll = length(test21))
predval = foremod@forecast$sigmaFor[1,]

# determining AR(p) of VIX -----------------

plot(VIX)
auto.arima(VIX["/2020-12-31"], max.q = 0)
acf(VIX["/2020-12-31"])
pacf(VIX["/2020-12-31"])
acf(diff(VIX["/2020-12-31"]), na.action = na.pass)
pacf(diff(VIX["/2020-12-31"]), na.action = na.pass)

autoplot(VIX, main = "VIX") /
(autoplot(acf(VIX["/2020-12-31"], lag.max = 25, plot = FALSE), main = "ACF of VIX")+
autoplot(pacf(VIX["/2020-12-31"], lag.max = 25, plot = FALSE), main = "PACF of VIX"))
autoplot(diff(VIX), main = "1st difference of VIX") /
(autoplot(acf(diff(VIX["/2020-12-31"]), lag.max = 25, na.action = na.pass, plot = FALSE), main = "ACF of differenced VIX")+
autoplot(pacf(diff(VIX["/2020-12-31"]), lag.max = 25, na.action = na.pass, plot = FALSE), main = "PACF of differenced VIX"))

conf.level <- 0.95
ciline <- qnorm((1 - conf.level)/2)/sqrt(length(VIX["/2020-12-31"]))
bacf <- acf(VIX["/2020-12-31"], plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))
ggplot(data=bacfdf, mapping=aes(x=lag, y=acf)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_hline(aes(yintercept = ciline), linetype = 2, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = 2, color = 'darkblue')
bpacf <- pacf(VIX["/2020-12-31"], plot = FALSE)
bpacfdf <- with(bpacf, data.frame(lag, acf))
ggplot(data=bpacfdf, mapping=aes(x=lag, y=acf)) +
  geom_bar(stat = "identity", position = "identity") + labs(y = "pacf") +
  geom_hline(aes(yintercept = ciline), linetype = 3, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'darkblue')

# forecast VIX using volatility of returns and lags of VIX -----------

tib = tibble(fortify.zoo(VIX)) %>% 
  left_join(x = ., y = fortify.zoo(sigma(fit)), by = "Index") %>% 
  left_join(x = ., y = fortify.zoo(SP500_Volume), by = "Index") %>% 
  left_join(x = ., y = fortify.zoo(TBill3m), by = "Index") %>%
  dplyr::rename(
    VIX = VIX.Close,
    SigmaReturns = "sigma(fit)",
    TradingVolume = GSPC.Volume,
    RiskFreeRate = DTB3
  ) %>% 
  timetk::tk_augment_lags(.value = c(TradingVolume, RiskFreeRate), .lags = 1) %>% 
  timetk::tk_augment_lags(.value = SigmaReturns, .lags = c(1:5), .names = "auto") %>% 
  timetk::tk_augment_lags(.value = VIX, .lags = c(1:5), .names = "auto") %>% 
  dplyr::filter(Index > "1990-01-10")

# mod with VIXlags ------------------

mod = lm(VIX ~ SigmaReturns + RiskFreeRate_lag1
         + VIX_lag1 + VIX_lag2 + VIX_lag3 + VIX_lag4 + VIX_lag5,
         data = tib %>% filter(Index <= "2021-01-01"))
summary(mod)

# sigma = c(tib$SigmaReturns[1:(nrow(tib)-length(predval[-1]))],predval[-1])

df = tib %>% 
  mutate(SigmaForecast = sigma) %>% 
  select(Index, VIX, RiskFreeRate_lag1, VIX_lag1:SigmaForecast) %>% 
  # timetk::tk_augment_lags(.value = SigmaForecast, .lags = c(1:5), .names = "auto") %>% 
  filter(Index >= "2021-01-01") %>%
  mutate(RollForecast = coef(mod)[1] + coef(mod)[2] * SigmaForecast
         + coef(mod)[3] * RiskFreeRate_lag1
# + coef(mod)[4] * SigmaForecast_lag1 + coef(mod)[5] * SigmaForecast_lag2
# + coef(mod)[6] * SigmaForecast_lag3 + coef(mod)[7] * SigmaForecast_lag4
# + coef(mod)[8] * SigmaForecast_lag5) 
         + coef(mod)[4] * VIX_lag1 + coef(mod)[5] * VIX_lag2
         + coef(mod)[6] * VIX_lag3 + coef(mod)[7] * VIX_lag4
         + coef(mod)[8] * VIX_lag5)

gg_S_for = ggplot(data = df) + geom_line(aes(x = Index, y = VIX, col = "Actual VIX")) +
  geom_line(aes(x = Index, y = RollForecast, col = "Rolling Forecast")) +
  labs(x = "Date", y = "Value", color = "", title = "Forecast only with VIX lags") + 
  theme(legend.position="top")
MSE_OutOfSample_withlags[1] = mean((df$VIX - df$RollForecast)^2)

# in sample
tib = tib %>% filter(Index <= "2021-01-01")
gg_S_fit = ggplot(data = tib) + geom_line(aes(x = Index, y = VIX, col = "Actual VIX")) +
  geom_line(aes(x = Index, y = fitted(mod), col = "Fitted Values")) +
  labs(x = "Date", y = "Value", color = "", title = "Fit only with VIX lags") + 
  theme(legend.position="top")
MSE_InSample_withlags[6] = mean((tib$VIX - fitted(mod))^2)

options(pillar.sigfig = 5) %>% 
  tibble(ModelName,MSE_InSample_withlags,MSE_OutOfSample_withlags) %>% View(.)






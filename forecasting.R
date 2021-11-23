library(tidyverse)
library(rugarch)
library(xts)

# rolling forecast

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

# forecast roll and static
test10y = returns["2011-11-01/2021-10-31"]
specs = ugarchspec(variance.model = list(garchOrder=c(1,1), model = "sGARCH"), 
                       mean.model = list(armaOrder=c(0,0), include.mean = T),
                       distribution.model = "norm")
garchfitfore = ugarchfit(data = returns, spec = specs, 
                         out.sample = length(test10y))
foremod = ugarchforecast(fitORspec = garchfitfore, data = returns, 
                         n.ahead = 10, n.roll = length(test10y), 
                         external.forecasts = list())
plot(foremod, which = 4)

plot.ts(foremod@model$modeldata$sigma)
plot.ts(foremod@forecast$sigmaFor[1,])

lm(VIX["/2011-10-31"] ~ sigma(garchfitfore))
plot.ts(VIX["2011-11-01/2021-10-31"], ylim = c(5,90))
lines(6.212 + 1365.612 * foremod@forecast$sigmaFor[1,], col = "steelblue")

predval = foremod@forecast$sigmaFor[1,]
tib = tibble(fortify.zoo(VIX))
tib = tib %>% filter(Index >= "2011-11-01") %>% 
    mutate(RollForecast = 6.212 + 1365.612 * predval[-1])

ggplot(data = tib) + geom_line(aes(x = Index, y = VIX.Close, col = "Actual VIX")) +
  geom_line(aes(x = Index, y = RollForecast, col = "Rolling Forecast")) +
  labs(x = "Date", y = "Value")



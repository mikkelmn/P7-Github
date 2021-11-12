library(tseries)

VIX = get.hist.quote('^GSPC', provider = 'yahoo', end = "2021-10-06")
plot(VIX$Close, main="S&P 500 ", xlab = "Date", ylab = "Price")

VIX.return = diff(log(VIX$Close))

plot(VIX.return, main = "S&P 500 returns",xlab = "Date", ylab = "Return")

acf(VIX.return, na.action = na.exclude, main = "ACF of S&P 500 returns")
pacf(VIX.return, na.action = na.exclude, main = "ACF of S&P 500 returns")

layout(matrix(c(1,2), nr=2, byrow=T))
acf(abs(VIX.return), na.action = na.exclude, main = "ACF of absolute S&P 500 returns")
acf(VIX.return^2, na.action = na.exclude, main = "ACF of squared S&P 500 returns")

m <- mean(VIX.return)
s <- sd(VIX.return)
hist(VIX.return, prob=TRUE, breaks=200, col="light blue", xlab = "S&P 500 returns", main = "Histogram of S&P 500 returns",
     xlim = c(-0.1, 0.1))
curve(dnorm(x, mean = m, sd = s), add = TRUE, lwd = 2, col = "red")


hist(VIX.return, prob=TRUE, breaks=200, col="light blue", xlab = "S&P 500 returns", main = "Histogram of S&P 500 returns", 
     xlim = c(-0.10, 0), ylim = c(0, 10))
curve(dnorm(x, mean = m, sd = s), add = TRUE, lwd = 2, col = "red")

e_plus = pmax(VIX.return, numeric(length(VIX.return)))
plot(e_plus)

e_min = pmin(VIX.return, numeric(length(VIX.return)))
plot(e_min)

ccf(abs(VIX.return), e_plus, na.action = na.remove)
ccf(abs(VIX.return), -e_min, na.action = na.remove)
lags = ccf(abs(VIX.return), VIX.return, na.action = na.remove)
str(lags)
lags$lag
plot(lags[lags$lag[36:71],])



# Kurtosis
library(DistributionUtils)
library(e1071)
kurtosis(VIX.return)




garch_spec = ugarchspec(variance.model = list(garchOrder=c(1,1), 
                                              model = "fGARCH", submodel = "TGARCH"), 
                        mean.model = list(armaOrder=c(0,0),
                                          include.mean = FALSE),
                        distribution.model = "std")
garch_fit = ugarchfit(spec=garch_spec, data=VIX.return,
                      solver.control=list(trace = 1))
plot(garch_fit)
garch_fit@fit$residuals

k = c()
for (i in 1:1000) {
  x <- rnorm(n = length(garch_fit@fit$sigma), sd = 1, mean = 0)
  eps = garch_fit@fit$sigma*x
  k[i] = kurtosis(eps)
}
mean(k)
x <- rnorm(n = length(garch_fit@fit$sigma), sd = 1, mean = 0)



adf.test(VIX.return, k=0)


rm(list = ls())
num <- 20
u <- runif(num)
l <- 2
y <- sqrt(-log(1-u) / l)
sumy2 <- sum(y^2)

# This sets the plot margins
par(mai=c(1.4, 1.4, 1, 0.5))

pdfy <- function(x,l) {
  2 * x * l * exp(-x*x*l)
}

mu    = 1/2 * sqrt(pi/l)    ## mean of pdfy
sigma = (4-pi)/(4*l)        ## Standard deviation

lambda = seq(0.0, 10.0, sigma/20)

plot(lambda, pdfy(y[1],lambda), type='l', lwd=2,
     main="Individual likelihoods per point",
     ylab=expression(paste("L(",lambda,"; y"[i],")")),
     xlab=expression(lambda))

for(i in 2:num)
  lines(lambda, pdfy(y[i],lambda), lwd=2)

## Likelihood for all the data
lik = 1
for(i in 1:num)
  lik = lik * pdfy(y[i],lambda)

plot(lambda, lik, type='l', lwd=2,
     main="Likelihood function",
     ylab=expression(paste("L(",lambda,"; y",")")),
     xlab=expression(lambda))
abline(v=num/(sumy2),col='blue',lwd=3,lty=2)

# Log likelihood for all the data
log.lik = 0
for(i in 1:num)
  log.lik = log.lik + log(pdfy(y[i],lambda))

plot(lambda, log.lik, type='l', lwd = 2,
     main = "Log-Likelihood Function",
     ylab = expression(paste("l(", lambda, " ; y)")),
     xlab = expression(lambda))
abline(v=num/(sumy2),col='blue',lwd=3,lty=2)
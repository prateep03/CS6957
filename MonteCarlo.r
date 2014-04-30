## Exponential distribution sampling, using Inverse Transform Method
my.rexp = function(n, lambda)
{
  u = runif(n)
  return(-log(u) / lambda)
}

## Test against pdf
lambda = 3
x = my.rexp(10000, lambda)
hist(x, freq=FALSE, main = "Exponential Distribution Sampling")
t = seq(0, max(x), 0.01)
lines(t, dexp(t, lambda), lwd=3, col='red')

## Gamma distribution sampling, using transform method
## k must be an integer
## We could generalize to non-integer k using acceptance-rejection, see
## http://en.wikipedia.org/wiki/Gamma_distribution
my.rgamma = function(n, k, theta)
{
  x = matrix(my.rexp(n * k, 1/theta), k, n)
  return(apply(x, MARGIN = 2, FUN = sum))
}

## Test against pdf
k = 5
theta = 0.5
x = my.rgamma(10000, k, theta)
hist(x, freq=FALSE, main = "Gamma Distribution Sampling")
t = seq(0, max(x), 0.01)
lines(t, dgamma(t, shape = k, scale = theta), lwd=3, col='red')

## Beta distribution sampling, using transform method
## alpha and beta must be integers (we are using the Gamma sampler above)
my.rbeta = function(n, alpha, beta)
{
  u = my.rgamma(n, alpha, 1)
  v = my.rgamma(n, beta, 1)

  return(u / (u + v))
}

## Test against pdf
a = 2
b = 5
x = my.rbeta(10000, a, b)
hist(x, freq=FALSE, main = "Beta Distribution Sampling")
t = seq(0, 1, 0.01)
lines(t, dbeta(t, a, b), lwd=3, col='red')

## Beta distribution with Acceptance-Rejection Sampling
## Works for a >= 1, b >= 1
## Note this is less efficient than transform method above
my.rbeta.rej = function(n, a, b)
{
  y = numeric(n)
  num.accepted = 0
  num.trials = 0
  while(num.accepted < n)
  {
    num.trials = num.trials + 1
    u = runif(1)
    x = runif(1)
    if(x^(a-1) * (1 - x)^(b-1) > u)
    {
      num.accepted = num.accepted + 1
      y[num.accepted] = x
    }
  }

  cat(paste("Acceptance rate =", n / num.trials, "\n"))
  return(y)
}

x = my.rbeta.rej(10000, a, b)
hist(x, freq=FALSE, main = "Beta Distribution Sampling (Acceptance-Rejection)")
t = seq(0, 1, 0.01)
lines(t, dbeta(t, a, b), lwd=3, col='red')

## Buffon's Needle Experiment to compute pi
## n = number of drops
## d = distance between lines
## L = length of needle
## Here we assume L <= d
## See: http://en.wikipedia.org/wiki/Buffon's_needle
buffon.needle = function(n, d, L)
{
  x = runif(n, min=0, max=(d/2))
  theta = runif(n, min=0, max=pi/2)

  k = sum( as.integer(x <= (L/2)*sin(theta)) )

  2 * L * n / (d * k)
}

## Area ratio experiment to compute pi
## n = number of points uniformly in square
area.pi = function(n)
{
  x = runif(n, min=-1, max=1)
  y = runif(n, min=-1, max=1)

  k = sum( as.integer( x^2 + y^2 <= 1 ) )

  4 * k / n
}

n = 10000

## Simple Monte Carlo Integration of x^2, from 0..1
## Correct answer is 1/3
u = runif(n)
mean(u^2)

## Compare to deterministic integration
u = seq(0, 1, 1 / n)
sum(u^2) / n

## Simple Monte Carlo Integration of x^2, from 2..4
## Correct answer is 18 + 2/3
u = runif(n, min=2, max=4)
2 * mean(u^2)

## Compare to deterministic integration
u = seq(2, 4, 2 / n)
2 * sum(u^2) / n

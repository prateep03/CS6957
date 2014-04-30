options(digits=6)

## Compute the power of a matrix
mpow = function(A, n)
{
  result = diag(1, nrow(A), ncol(A))
  while(n > 0)
  {
    if(n %% 2 != 0)
    {
      result = result %*% A
      n = n - 1
    }
    else
    {
      A = A %*% A
      n = n / 2
    }
  }

  result
}

## Prints the sequence T^k * p, for k = 1..n
markov.seq = function(T, n = 10, p = runif(3))
{
  p = p / sum(p)

  cat("p = "); cat(p); cat("\n")

  for(k in 1:n)
  {
    cat(paste("T^", k, " * p = ", sep=''))
    cat(mpow(T, k) %*% p)
    cat("\n")
  }
}

## Prints eigenanalysis of transition kernel T
markov.eig = function(T)
{
  e = eigen(T)

  ## 1st eigenvector is stationary distribution
  v = e$vectors[,1]
  cat("Stationary dist = ")
  cat(abs(v) / sum(abs(v)))
  cat("\n")

  ## 1st eigenvalue is always 1
  cat("Eigenvalues = ")
  cat(e$values)
  cat("\n")

  ## 2nd eigenvalue (absolute value) is rate of convergence
  cat("Convergence rate = ")
  cat(abs(e$values[2]))
  cat("\n")
}

## Transition kernel
T1 = t(matrix(c(
  0.5,  0.25, 0.25,
  0.25, 0.5,  0.25,
  0.25, 0.25, 0.5
  ), 3, 3))

## Different initial probabilities converge to same stationary distribution
markov.seq(T1)
markov.seq(T1)
markov.seq(T1)

## Eigenanalysis
markov.eig(T1)



## Example in Andrieu, et al.
T2 = t(matrix(c(
  0, 0,   0.6,
  1, 0.1, 0.4,
  0, 0.9, 0
  ), 3, 3))

## Different initial probabilities converge to same stationary distribution
## Notice the convergence is slower than T1
markov.seq(T2, 40)
markov.seq(T2, 40)
markov.seq(T2, 40)

## Eigenanalysis
markov.eig(T2)



## Example of periodic transition
T3 = t(matrix(c(
  0, 0, 1,
  1, 0, 0,
  0, 1, 0
  ), 3, 3))

## Notice matrix powers cycle
mpow(T3, 10)
mpow(T3, 11)
mpow(T3, 12)
mpow(T3, 13)

## Eigenanalysis
markov.eig(T3)


## Continuous example: AR(1)
## Generate a random sequence of n AR(1)
rar = function(n, theta = 0, sigma = 1, x1 = 0)
{
  x = numeric(n)
  x[1] = x1
  for(i in 2:n)
    x[i] = x[i - 1] * theta + rnorm(1, 0, sigma)

  return(x)
}

theta = 0.5
sigma = 1

x = rar(200, theta = theta, sigma = sigma)
plot(x, type='l', main="AR(1) Sequence")

## AR(1) sequence starting from end of last one
x = rar(10000, theta = theta, sigma = sigma, x1 = x[200])

## Print stats
cat(paste("Mean =", mean(x), "\n"))
cat(paste("Variance =", var(x), "\n"))
stat.var = sigma^2 / (1 - theta^2)
cat(paste("Stationary variance =", stat.var, "\n"))

## Compare histogram to stationary distribution
hist(x, freq=FALSE)
s = sqrt(stat.var)
t = seq(-3*s,3*s,0.01*s)
lines(t, dnorm(t, 0, s), col='red')

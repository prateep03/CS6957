require(animation)

## Code from Radford Neal's paper (slightly modified, adds animation)
HMC = function (U, grad_U, epsilon, L, current_q, anim = FALSE)
{
  q = current_q
  p = rnorm(length(q),0,1) # independent standard normal variates
  current_p = p

  qList = q
  pList = p

  ## Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2

  ## Alternate full steps for position and momentum
  for(i in 1:(L-1))
  {
    ## Make a full step for the position
    q = q + epsilon * p

    ## Make a full step for the momentum, except at end of trajectory
    p = p - epsilon * grad_U(q)

    if(anim)
    {
      H = function(x, y) { U(x) + 0.5 * y^2 }

      t = seq(-2, 2, 0.01)
      Hvals = outer(t, t, H)
      contour(t, t, Hvals, levels=seq(0, 3, 0.1), asp=1,
              drawlabel=FALSE, xlab="q", ylab="p")
      lines(qList, pList, lwd=2, col='red')
      ani.record()
      pList = c(pList, p)
      qList = c(qList, q)
    }
  }

  q = q + epsilon * p

  ## Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2

  ## Negate momentum at end of trajectory to make the proposal symmetric
  p = -p

  ## Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2

  ## Accept or reject the state at end of trajectory, returning either
  ## the position at the end of the trajectory or the initial position
  if(current_U + current_K > proposed_U + proposed_K)
    return(q) # reject
  else if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
    return (q) # accept
  else
    return (current_q) # reject
}

## This is the same unormalized pdf we sampled from in MetropolisHastings.r
f = function(x)
{
  exp(-x^8 + 4 * x^4 - 3 * x^2)
}

clamp = function(x, a, b) { (x<a)*a + (x>b)*b + (x>=a)*(x<=b)*x }

## U(q) and its derivative
u = function(q) { q^8 - 4 * q^4 + 3 * q^2 }
du = function(q) { clamp(8 * q^7 - 16 * q^3 + 6 * q, -1e8, 1e8) }


trace.plot = function(q)
{
  plot(1:length(q), q, type='l')
}

hist.plot = function(q)
{
  t = seq(-2, 2, 0.01)

  ## Approximate the normalizing constant (for purposes of plotting the pdf)
  C = sum(f(t)) * 0.01

  hist(q, freq=FALSE, breaks=80)
  t = seq(-2, 2, 0.01)
  lines(t, f(t) / C, lwd = 3, col = 'red')
}

crazy.example = function(numBurn = 10, numSamples = 1000, epsilon = 0.05, L = 40)
{
  q = 0

  ## Burn-in
  for(i in 1:numBurn)
    q = HMC(u, du, epsilon, L, q)

  rejectCount = 0
  ## Sampling
  for(i in 1:(numSamples-1))
  {
    q = c(q, HMC(u, du, epsilon, L, q[i]))
    if(q[i] == q[i+1])
      rejectCount = rejectCount + 1
  }

  cat(paste("Acceptance rate =", 1 - rejectCount / numSamples, "\n"))

  q
}

crazy.anim = function(epsilon = 0.05, L = 100, q = 0.5)
{
  ## Now we'll animate the Hamiltonian integration
  ani.options(interval = 0.1)
  ani.start()
  HMC(u, du, epsilon, L, q, anim=TRUE)
  ani.stop()
}

# U(\beta)
U = function(X, Y, Be, sigma)
{
  U = sum((1-Y) * (X %*% Be) + log(1 + exp(-X %*% Be))) + sum(Be * Be) / (2*sigma*sigma);
  return(U);
}

# \nabla U(\beta)
nablaU = function(X, Y, Be, sigma)
{
  res = colSums((1-Y)*X + (c(exp(-X %*% Be)) * -X) / c(1+exp(-X %*% Be)) ) + Be / (sigma^2);
  return(res);
}
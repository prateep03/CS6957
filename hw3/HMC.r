require(animation)
require(package=numDeriv)

## Code from Radford Neal's paper (slightly modified, adds animation)
HMC = function (X, Y, epsilon, L, current_q, sigma, anim = FALSE)
{
  q = current_q
  p = rnorm(length(q),0,1) # independent standard normal variates
  current_p = p

  qList = q
  pList = p

  epsilon1 = hess.inverse(X,Y,q,sigma);
  #cat(paste("Ep = ", epsilon1, "\n"));
  
  ## Make a half step for momentum at the beginning
  p = p - epsilon1 %*% grad_U(X, Y, q, sigma) / 2

  ## Alternate full steps for position and momentum
  for(i in 1:(L-1))
  {
    ## Make a full step for the position
    q = q + epsilon * p

    ## Make a full step for the momentum, except at end of trajectory
    p = p - epsilon1 %*% grad_U(X,Y,q,sigma)

    if(anim)
    {
      H = function(x, y) { U(X, Y, x, sigma) + 0.5 * y^2 }

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
  p = p - epsilon1 %*% grad_U(X, Y, q, sigma) / 2

  ## Negate momentum at end of trajectory to make the proposal symmetric
  p = -p

  ## Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(X, Y, current_q, sigma)
  current_K = sum(current_p^2) / 2
  proposed_U = U(X, Y, q, sigma)
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

crazy.anim = function(epsilon = 0.05, L = 100, q = 0.5)
{
  ## Now we'll animate the Hamiltonian integration
  ani.options(interval = 0.1)
  ani.start()
  HMC(u, du, epsilon, L, q, anim=TRUE)
  ani.stop()
}

clamp = function(x, a, b) { (x<a)*a + (x>b)*b + (x>=a)*(x<=b)*x }

# \hess U(\beta)
hess.inverse = function(X,Y,Be,sigma)
{
  #    res = sum( (c(exp(-X%*%Be)) * X*X) / c((1+exp(-X %*% Be))^2)) + 1/(sigma^2);
  #   res = (c(exp(-X%*%Be)) * X*X) / c((1+exp(-X %*% Be))^2);
  res = hessian(func = U, x = Be, method.args=as.list(X = X,Y = Y,Be = Be,sigma = sigma));
  #   res = res + 1/(sigma^2);
  return(chol2inv(res));
}

# U(\beta)
U = function(X, Y, Be, sigma)
{
  res = sum((1-Y) * (X %*% Be) + log(1 + exp(-X %*% Be))) + sum(Be * Be) / (2*sigma*sigma);
  return(res);
}

# \grad U(\beta)
grad_U = function(X, Y, Be, sigma)
{
  res = colSums((1-Y)*X + (c(exp(-X %*% Be)) * -X) / c(1+exp(-X %*% Be)) ) + Be / (sigma^2);
  return(clamp(res,-1e8,1e8));
}

## Utility functions
trace.plot = function(q,col_name)
{
  plot(1:length(q), q, type='l',xlab ="Number of Samples",ylab=col_name);
}
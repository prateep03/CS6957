# Function to compute covariance
covariance = function(x1, x2, lamda, eps)
{
  l1 = length(x1);
  l2 = length(x2);
  res = matrix(0.0, l1, l2);
  
  for(i in 1:l1)
  {
    tmp = x2 - x1[i];
    res[i,] = exp(-0.5*tmp*tmp/lamda) + eps*1.0*(abs(tmp)==0);
  }
  return(res)
}

# Draw sample from data
drawsample = function(x,lamda, eps)
{
  num = length(x);
  Sigma = covariance(x, x, lamda, eps);
  lt = t(chol(Sigma));
  u = rnorm(1:num, 0, 1);
  res = lt %*% u;
  return(res)
}
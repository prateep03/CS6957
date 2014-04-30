## Function to create a conditional probability table
## varnames: vector of variable names (strings)
## probs: vector of probabilities for the flattened probability table
## levelsList: a list containing a vector of levels (outcomes) for each variable
## See the BayesNetExamples.r file for examples of how this function works
createCPT = function(varnames, probs, levelsList)
{
  ## Check dimensions agree
  if(length(probs) != prod(sapply(levelsList, FUN=length)))
    return(NULL)

  ## Set up table with appropriate dimensions
  m = length(probs)
  n = length(varnames)
  g = matrix(0, m, n)

  ## Convert table to data frame (with column labels)
  g = as.data.frame(g)
  names(g) = varnames

  ## This for loop fills in the entries of the variable values
  k = 1
  for(i in n:1)
  {
    levs = levelsList[[i]]
    g[,i] = rep(levs, each = k, times = m / (k * length(levs)))
    k = k * length(levs)
  }

  return(data.frame(probs = probs, g))
}

## Build a CPT from a data frame
## Constructs a conditional probability table as above, but uses frequencies
## from a data frame of data to generate the probabilities.
createCPT.fromData = function(x, varnames)
{
  levelsList = list()

  for(i in 1:length(varnames))
  {
    name = varnames[i]
    levelsList[[i]] = sort(unique(x[,name]))
  }

  m = prod(sapply(levelsList, FUN=length))
  n = length(varnames)
  g = matrix(0, m, n)

  ## Convert table to data frame (with column labels)
  g = as.data.frame(g)
  names(g) = varnames

  ## This for loop fills in the entries of the variable values
  k = 1
  for(i in n:1)
  {
    levs = levelsList[[i]]
    g[,i] = rep(levs, each = k, times = m / (k * length(levs)))
    k = k * length(levs)
  }

  ## This is the conditional probability column
  probs = numeric(m)
  numLevels = length(levelsList[[1]])
  skip = m / numLevels

  ## This chunk of code creates the vector "fact" to index into probs using
  ## matrix multiplication with the data frame x
  fact = numeric(ncol(x))
  lastfact = 1
  for(i in length(varnames):1)
  {
    j = which(names(x) == varnames[i])
    fact[j] = lastfact
    lastfact = lastfact * length(levelsList[[i]])
  }
  ## Compute unnormalized counts of subjects that satisfy all conditions
  a = as.matrix(x - 1) %*% fact + 1
  for(i in 1:m)
    probs[i] = sum(a == i)

  ## Now normalize the conditional probabilities
  for(i in 1:skip)
  {
    denom = 0 ## This is the normalization
    for(j in seq(i, m, skip))
      denom = denom + probs[j]
    for(j in seq(i, m, skip))
    {
      if(denom != 0)
        probs[j] = probs[j] / denom
    }
  }

  return(data.frame(probs = probs, g))
}

## Product of two factors
## A, B: two factor tables
##
## Should return a factor table that is the product of A and B.
## You can assume that the product of A and B is a valid operation.
productFactor = function(A, B)
{
  v <- intersect(names(A), names(B));
  u <- union(names(A), names(B));
  vvars <- v[v != "probs"];
  f <- merge(A,B,by=vvars);
  fvars <- names(f);
  uvars <- u[u != "probs"];
  vars <- subset(f, select=c(uvars));
  probvars <- setdiff(fvars, uvars);
  pvars <- f[,c(probvars)];
  pvars <- transform(pvars, n = pvars[,1] * pvars[,2]);
  prodProb <- cbind(pvars$n, vars);
  names(prodProb) <- c("probs", uvars);
  return(prodProb);
}

marginalizeFactor = function(X, margVar)
{
  varNames = setdiff(names(X), margVar)
  # generate model <- probs + X\margVar
  f = as.formula(paste("probs ~ ", paste(varNames, collapse= "+")))
  # sum over the patterns in X according to model f
  marginal.X = aggregate(f, data=X, FUN=sum)
  # bring probs at the beginning
  marginal.X = marginal.X[, length(marginal.X):1]
  return(marginal.X)
}

# marginalize : I discussed this with Sayan Dey
marginalize = function(net, margVars)
{
  if(length(margVars) != 0)
  {
	  for(i in 1:length(margVars))
	  {
		  colref = lapply(net, function(x) names(x))
		  cf = lapply(colref, function(x) grepl(margVars[i],x))
		  cfl=lapply(cf, function(x) sum(x)>0)
      
		  sub_net=net[which(cfl==TRUE)]
		  net[which(cfl==TRUE)]<-NULL	
		  prod_factor=data.frame(sub_net[[1]])
		  if(length(sub_net)>1)
		  {
			  for(j in 2:length(sub_net))
			  {
				  prod_factor=productFactor(prod_factor,data.frame(sub_net[[j]]))	
			  }
		  }
		  final_factor=marginalizeFactor(prod_factor,margVars[i])
		  net[[length(net)+1]]=data.frame(final_factor)				
	  }
	  val=lapply(net,function(y) sum(as.numeric(apply(y,1,function(x) (x[1])))))
   
    val_1=lapply(net,function(x) nrow(x))
	  check=list()
    for(i in 1:length(val))
    {
      if(val[[i]]==val_1[[i]])
      {
        check[i]=TRUE
      }
      else
        check[i]=FALSE
    }
   
	  net=net[which(check==FALSE)]
    
  }
  
  return(net)  
}

observe = function(bayesnet, obsVars, obsVals)
{
  len_bayes = length(bayesnet)
  for (i in 1:len_bayes) { # number of factors
    commonVar = intersect(obsVars, names(as.data.frame(bayesnet[i])))
    len_Var = length(commonVar)
    if (len_Var != 0 ) {
      for (j in 1:len_Var) {
        obserV = c(commonVar[j], obsVals[which(obsVars == commonVar[j])])
        obserV = paste(obserV, collapse="==")
        bayesnet[i] = list(subset(as.data.frame(bayesnet[i]),eval(parse(text = obserV))))
      }      
    } 
  }
  return(bayesnet)
}

infer = function(bayesnet, margVars, obsVars, obsVals)
{
  bayesnet = marginalize(bayesnet, margVars)
  bayesnet = observe(bayesnet, obsVars, obsVals)
  # renormalizing and return the joint and probability
  for (i in 1:length(bayesnet)) {
    if (i == 1) {ret.joint = as.data.frame(bayesnet[i])}
    else {ret.joint = productFactor(ret.joint, as.data.frame(bayesnet[i]))}
  }
  ret.joint[, 1] = ret.joint[, 1]/sum(ret.joint[, 1])
  return(ret.joint)
}

##
## Q3. Function to generate health outcome given income, and marginalize
## rest of the RVs.
income.healthOutcome = function(bayesnet, margVars, obsVars)
{
  res = 1:8
  for(i in 1:8) {
    r = infer(bayesnet, margVars, obsVars, c(i));
    res[i] = r[1,1];
  }
  return(res);   
}
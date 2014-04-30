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
	header_A<-names(A);	
	header_B<-names(B);
	value<-intersect(header_A,header_B);
	value<-value[-c(1)];
	final_header<-union(header_A,header_B);
	vars<-final_header[-c(1)];
	merge_tables<-merge(A,B,by=value);

	var_table<-subset(merge_tables,select=c(vars));
	names_merge<-names(merge_tables);
	names_need<-setdiff(names_merge,vars);	
	merge_vals<-merge_tables[,c(names_need)];	
	merge_vals<-transform(merge_vals,new=merge_vals[,1]*merge_vals[,2]);
	final_table<-cbind(merge_vals$new,var_table);
	names(final_table)<-c("probs",vars);
	
	return(final_table);
}

## Marginalize a variable from a factor
## A: a factor table
## margVar: a string of the variable name to marginalize
##
## Should return a factor table that marginalizes margVar out of A.
## You can assume that margVar is on the left side of the conditional.
marginalizeFactor = function(X, margVar)
{
	
	B<-margVar
	names_var<-names(X)
	names_need<-setdiff(names_var,B)
	form<-paste(names_need, collapse= "+")
	form<-paste("probs ~ ",form)
	fmla = as.formula(form)
  marg = aggregate(fmla,FUN=sum,data=X)
  marg = marg[,length(marg):1]
	return(marg);

}

## Marginalize a list of variables
## bayesnet: a list of factor tables
## margVars: a vector of variable names (as strings) to be marginalized
##
## Should return a Bayesian network (list of factor tables) that results
## when the list of variables in margVars is marginalized out of bayesnet.
marginalize = function(bayesnet, margVars)
{
  if(!is.null(margVars))
  {
	  for(i in 1:length(margVars))
	  {
		  colref <- lapply(bayesnet, function(x) names(x))
		  cf <- lapply(colref, function(x) grepl(margVars[i],x))
		  cfl <- lapply(cf, function(x) sum(x)>0)
      
		  sub_bayesnet<-bayesnet[which(cfl==TRUE)]
		  bayesnet[which(cfl==TRUE)]<-NULL	
		  prod_factor<-data.frame(sub_bayesnet[[1]])
		  if(length(sub_bayesnet)>1)
		  {
			  for(j in 2:length(sub_bayesnet))
			  {
				  prod_factor<-productFactor(prod_factor,data.frame(sub_bayesnet[[j]]))	
			  }
		  }
		  final_factor<-marginalizeFactor(prod_factor,margVars[i])
		  bayesnet[[length(bayesnet)+1]]<-data.frame(final_factor)				
	  }
	  val<-lapply(bayesnet,function(y) sum(as.numeric(apply(y,1,function(x) (x[1])))))
   
    val_1<-lapply(bayesnet,function(x) nrow(x))
	  check<-list()
    for(i in 1:length(val))
    {
      if(val[[i]]==val_1[[i]])
      {
        check[i]<-TRUE
      }
      else
        check[i]<-FALSE
    }
   
	  bayesnet<-bayesnet[which(check==FALSE)]
    
  }
  
  return(bayesnet)  
}

## Observe values for a set of variables
## bayesnet: a list of factor tables
## obsVars: a vector of variable names (as strings) to be observed
## obsVals: a vector of values for corresponding variables (in the same order)
##
## Set the values of the observed variables. Other values for the variables
## should be removed from the tables. You do not need to normalize the factors
## to be probability mass functions.
observe = function(bayesnet, obsVars, obsVals)
{
  
  if(!is.null(obsVars))
  {
    for(i in 1:length(obsVars))
    {
      colref <- lapply(bayesnet, function(x) names(x))
      cf <- lapply(colref, function(x) grepl(obsVars[i],x))
      cfl <- lapply(cf, function(x) sum(x)>0)
      sub_bayesnet<-bayesnet[which(cfl==TRUE)]
      bayesnet[which(cfl==TRUE)]<-NULL	
      for(j in 1:length(sub_bayesnet))
      {
        sub_bayes<-sub_bayesnet[[j]];
        cf <- lapply(sub_bayes[obsVars[i]],function(x) grepl(obsVals[i],x))
        sub<-lapply(cf,function(x) sub_bayes[which(x==TRUE),])
        bayesnet[length(bayesnet)+1]<-sub
      }
    }
  }
  return(bayesnet)
}

## Run inference on a Bayesian network
## bayesnet: a list of factor tables
## margVars: a vector of variable names to marginalize
## obsVars: a vector of variable names to observe
## obsVals: a vector of values for corresponding variables (in the same order)
##
## This function should run marginalization and observation of the sets of
## variables. In the end, it should return a single joint probability table. The
## variables that are marginalized should not appear in the table. The variables
## that are observed should appear in the table, but only with the single
## observed value. The variables that are not marginalized or observed should
## appear in the table with all of their possible values. The probabilities
## should be normalized to sum to one.
infer = function(bayesnet, margVars, obsVars, obsVals)
{
  bayesnet<-marginalize(bayesnet,margVars)
  bayesnet<-observe(bayesnet,obsVars, obsVals)
  #print(bayesnet)
  prod_factor<-data.frame(bayesnet[[1]])
  if(length(bayesnet)>1)
  {
    for(j in 2:length(bayesnet))
    {
      prod_factor<-productFactor(prod_factor,data.frame(bayesnet[[j]]))	
    }
  }  
  sum_val<-sum(prod_factor$probs)
  new_factor<-lapply(prod_factor$probs,function(x) x/sum_val)
  prod_factor$probs<- new_factor
  return(prod_factor)
}

##
## Q3. Function to generate health outcome given income, and marginalize
## rest of the RVs.
income.healthOutcome = function(bayesnet, margVars, obsVars)
{
  res = array(1:8);
  for(i in 1:8) {
    r = infer(bayesnet, margVars, obsVars, c(i));
    res[i] = r[1,1];
  }
  return(res);   
}
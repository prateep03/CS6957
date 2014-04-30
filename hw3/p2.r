rm(list = ls())
plot.new()
source('C:/Users/prateep/Documents/CS6957/hw3/HMC.r');

## X
data.versi = as.data.frame.matrix(subset(iris, iris$Species==c("versicolor")));
data.virginica = as.data.frame.matrix(subset(iris, iris$Species==c("virginica")));
data.train <- rbind(data.versi[1:30,], data.virginica[1:30,]);
data.test <- rbind(data.versi[31:50,], data.virginica[31:50,]);

numtrains <- nrow(data.train);
numtests <- nrow(data.test);

## Y
label.train = 1:numtrains;
label.train[which(data.train$Species == "versicolor")] = 0;
label.train[which(data.train$Species == "virginica")] = 1;
label.test = 1:numtests;
label.test[which(data.test$Species == "versicolor")] = 0;
label.test[which(data.test$Species == "virginica")] = 1;

# prepare X
data.train = data.train[-grep('Species',colnames(data.train))];
data.train = cbind(data.train, matrix(1,numtrains,1));
data.train = as.matrix(data.train);

data.test = data.test[-grep('Species', colnames(data.test))];
data.test = cbind(data.test, matrix(1,numtests,1));
data.test = as.matrix(data.test);

## HMC
L = 10;
numBurn = 100;
dim = 5;
current_q = matrix(0, dim, 1);
epsilon = 0.05;
sigma = 1.0;
numSamples = 1000;

q = c(-1,0,-1,0,-1)
my=matrix(rep(0,5),1,5)

## Burn-in
alp=matrix(rep(0,5),numSamples,5)
for(i in 1:numBurn)
  q = HMC(data.train, label.train, epsilon, L, q, sigma)

my=q
my1=t(my)
rejectCount = 0 
## Sampling
alp[1,]=my1
for(i in 1:(numSamples-1))
{ 
  q=  HMC(data.train, label.train, epsilon, L, alp[i,], sigma)
  my=q
  my1=t(my)
  alp[i+1,] =my1
  
  if(all(alp[i,] == alp[i+1,]))
    rejectCount = rejectCount + 1
}

#print(alp)
cat(paste("Acceptance rate =", 1 - rejectCount / numSamples, "\n"))

## Computing predictive posterior probability
y.predict = matrix(-1, numtests, 1);
for(i in 1:numtests)
{
  sum.predict = 0;
  for(j in 1:numSamples)
  {
    q = as.matrix(alp[j,]);
    x = data.test[i,];
    de = x%*%q;
    pr = 1/(1+exp(-de));
    sum.predict = sum.predict + pr;
  }
  prob.predict = sum.predict / numSamples;
  #cat(paste("Posterior predictive probability for test sample", i, " is ", prob.predict, "\n" ));
  if(prob.predict >= 0.5) 
  {
    y.predict[i]=1;
  } 
  else
  {
    y.predict[i] = 0;
  }
}

colname = c("Sepal Length","Sepal Width","Petal Length","Petal Width","Constant");
for(d in 1:dim)
{
  q <- alp[,d];
  c <- colname[d];
  trace.plot(q,c);
}

## average error rate
error = 0;
for(i in 1:numtests)
{
  if(y.predict[i] != label.test[i])
  {
    print(i)
    error = error + 1;
  }
}

cat(paste("Average error rate =", error/numtests, "\n"))
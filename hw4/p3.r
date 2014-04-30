require(MASS)

rm(list=ls())
source('~/CS6957/hw4/GP.r')
#plot.new();

data = read.csv("SaltLakeTemperatures.csv",nrows=-1);
data["YEAR"] = 0;
data["MONTH"] = 0;
dummy = lapply(data$DATE, function(x) substr(x,1,4));
for(i in 1:length(dummy))
{
  data$YEAR[i] = as.numeric(x=dummy[i]);
}
dummy = lapply(data$DATE,function(x) substr(x,5,6));
for(i in 1:length(dummy))
{
  data$MONTH[i] = as.numeric(x=dummy[i]);
}

x = data$MONTH;
num = length(x);
sigma = 0.5;
lamda = 1.0;
cov.eps = 1e-5;

#################### (a) ##################################
y = data$AVEMAX;
x.test = seq(1,12,length.out=num);

ka = covariance(x,x,lamda,cov.eps);
ka = ka + sigma^2*diag(num);
ka.inv = qr.solve(ka);
kb = covariance(x,x.test,lamda,cov.eps);
kc = covariance(x.test,x,lamda,cov.eps);
kd = covariance(x.test,x.test,lamda,cov.eps);

post.mu <- kc %*% ka.inv %*% y;
post.cov <- kd - kc %*% ka.inv %*% kb;
post.cov <- diag(post.cov);

## Plot 95% confidence interval
t = seq(1,12,length.out=num);
up = post.mu + 1.96 * sqrt(post.cov);
down = post.mu - 1.96 * sqrt(post.cov);
t.rev = t[length(t):1];
down = down[length(down):1];

plot(x,y, xlim=c(1,12),ylim=range(y,down,up),
     main=expression("Average temperature trajectory"),
     xlab="months",ylab=expression(paste(mu,"*")))
polygon(x=c(t,t.rev), y=c(up,down), col="grey", border=NA)
lines(t, post.mu, col = 'blue', lwd=3)
legend("topright", c("GP Regression", "Predictive Interval"),
       col = c("blue", "grey"), lwd=c(3,3,10),cex=0.7)
#dev.copy2pdf(file = "p3a.pdf")
###########################################################

print('Done (a)');
#plot.new();

######################### (b) #############################
x.test = seq(1,12,length.out=num);
z = data$YEAR;
y = data$AVEMAX;
z = z - min(z); 
#z = rep(x=1,num);
k11 = covariance(x,x,lamda,cov.eps);
a = sweep(x=k11,MARGIN=2,STATS=z,FUN='*');
b = sweep(x=a,MARGIN=1,STATS=z,FUN='*');
k11 = b + k11 + sigma^2*diag(num);
k11.inv = solve(k11);

k12 = covariance(x,x.test,lamda,cov.eps);
k13 = covariance(x,x.test,lamda,cov.eps);
k13 = sweep(x=k13,MARGIN=1,STATS=z,FUN='*');

k21 = covariance(x.test,x,lamda,cov.eps);
k22 = covariance(x.test,x.test,lamda,cov.eps);
k23 = matrix(0.0,num,num);
k31 = covariance(x.test,x,lamda,cov.eps);
k31 = sweep(x=k31,MARGIN=2,STATS=z,FUN='*');

k32 = matrix(0.0,num,num);
k33 = k22;

post.mu <- rbind(k21,k31) %*% k11.inv %*% y;
#post.cov <- rbind(cbind(k22,k23),cbind(k32,k33)) - ( rbind(k12,k13) %*% ( k11.inv %*% cbind(k21,k31)) );
post.cov <- rbind(cbind(k22,k23),cbind(k32,k33)) - ( rbind(k21,k31) %*% ( k11.inv %*% cbind(k12,k13)) );
post.cov <- (diag(post.cov));

post.mu.f1 <- post.mu[(num+1L):(length(post.mu)),]; #*600.0;
post.cov.f1 <- post.cov[(num+1L):(length(post.cov))]; #*16000.0;
post.mu.f0 <- post.mu[1:num];
post.cov.f0 <- post.cov[1:num]; #*700.0;

# ## Plot 95% confidence interval
t = seq(1,12,length.out=num);
up = post.mu.f1 + 1.96 * sqrt(post.cov.f1);
down = post.mu.f1 - 1.96 * sqrt(post.cov.f1);
t.rev = t[length(t):1];
down = down[length(down):1];

# plot(x,y,xlim=c(1,12),ylim=range(y,down,up),
#      main=expression(paste("Posterior mean for f"[0]," with 95% confidence")),
#      xlab="months",ylab=expression(paste(mu,"*")));
plot(t,rep(0,num),ylim=range(post.mu.f1,down,up),type="l",lty=2,
     main=expression(paste("Posterior mean for f"[1]," with 95% confidence")),
     xlab="months",ylab=expression(paste(mu,"*")))
polygon(x=c(t,t.rev), y=c(up,down), col="grey", border=NA)
lines(t, post.mu.f1, col = 'blue', lwd=2)
#lines(t, post.mu.f0, col = 'red', lwd=2)

legend("topright", c("GP Regression","Predictive Interval"),
       col = c("blue","grey"), lwd=c(3,3,10),cex=0.7)
#dev.copy2pdf(file = "p3b_1.pdf")
###########################################################
rm(list = ls())
source('~/CS6957/hw4/GP.r')
#plot.new();

num = 50;
sigma = 0.5;
lamda = 1.5;
cov.eps = 1e-5;

x = seq(-15,15,length.out=num);
#x = runif(num,min = 0,max = 2*pi);
x.test = seq(-15,15,length.out=num);
eps = rnorm(1:num,0,sigma);
y = sin(x.test) + eps;

ka = covariance(x,x,lamda=lamda,eps=cov.eps);
ka = ka + sigma*diag(num);
ka.inv = qr.solve(ka);
kb = covariance(x,x.test,lamda,cov.eps);
kc = covariance(x.test,x,lamda,cov.eps);
kd = covariance(x.test,x.test,lamda,cov.eps);

post.mu <- kc %*% ka.inv %*% y;
post.cov <- kd - kc %*% ka.inv %*% kb;
post.cov <- diag(post.cov);

## Plot 95% confidence interval
t = seq(0,2*pi,length.out=num);
up = post.mu + 1.96 * sqrt(post.cov);
down = post.mu - 1.96 * sqrt(post.cov);
t.rev = t[length(t):1];
down = down[length(down):1];

plot(x,y, xlim=c(0,2*pi),ylim=range(y,down,up),
     main="95% Posterior Predictive Interval")
polygon(x=c(t,t.rev), y=c(up,down), col="grey", border=NA)
lines(t, post.mu, col = 'blue', lwd=3)
lines(t, sin(x), col='red', lwd=3)
legend("bottomright", c("True Answer", "GP Regression", "Predictive Interval"),
      col = c("red", "blue", "grey"), lwd=c(3,3,3,10),cex=0.6)
#dev.copy2pdf(file = "p2.pdf")
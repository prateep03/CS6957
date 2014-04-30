rm(list = ls())
source('~/CS6957/hw4/GP.r')
#plot.new();

#lamda = 0.1;
#plot(x,y)
eps = 1e-5;
x = seq(-5.0, 5.0, by=0.01);
y = drawsample(x, 0.01, eps);
plot(x,y,type='l', lwd=3,col="green",
     ylim =c(-3,3),
     main=c("Gaussian process samples"))
y = drawsample(x, 0.01, eps);
lines(x,y,col="red",lwd=3)
y = drawsample(x, 1, eps);
lines(x,y,col="blue",lwd=3)
y = drawsample(x, 10, eps);
lines(x,y,col="black",lwd=3)
#y = drawsample(x, 100, eps);
#lines(x,y,col="magenta",lwd=3)

legend('bottomright',
       c(expression(paste(lambda,"= 0.01")),
         expression(paste(lambda,"= 0.1")),
         expression(paste(lambda,"= 1")),
         expression(paste(lambda,"= 10"))),
       col=c("green","red","blue","black"),lwd=3,cex=0.6)
#dev.copy2pdf(file = "p1.pdf")
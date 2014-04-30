# program 4
rm(list = ls()) # clear workspace
num <- 10000
u <- runif(num)
l <- 2
y <- sqrt(-log(1-u) / l)

yHist <- hist(y, freq = F) 
hist(y, freq = F, border = "dark blue", col = "light blue",
              main = "Histogram of random variable y", xlab = "y")
x <- seq(0.0, 3.0, length.out = num)
pdfy <- 2 * exp(-x*x*l) * x * l
#pdfy <- pdfy / max(pdfy)
lines(x, pdfy, col = "red", lwd = 2)
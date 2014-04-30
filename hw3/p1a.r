rm(list = ls())
plot.new()
source('C:/Users/prateep/Documents/CS6957/hw3/MRF.r');

# Original image
filename = 'C:./noisy-message.png';
oim = read.image(filename);
D = dim(oim);

# Gibbs sampling with Ising model
alpha = 0.01
beta = 0.85
iter = 200
nbrhood = nbrMask(4);

M = D[1];
N = D[2];
# crazy.anim(oim, M, N, nbrhood, iter, alpha, beta, 1.00);

##random image
# M = 30;
# N = 40;
# im = (matrix(runif(M*N),M,N) > 0.5)*2-1;
# display.image(im);
rim = ising.prior(oim, M, N, nbrhood, iter, alpha, beta);
display.image(rim);

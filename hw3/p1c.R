rm(list = ls())
#plot.new()
source('C:/Users/prateep/Documents/CS6957/hw3/MRF.r');

# Original image
filename = 'C:./noisy-message.png';
oim = read.image(filename);
D = dim(oim);

# Gibbs sampling with Ising model
alpha = 0.01
beta = 0.85
iter = 100
sigma = 1.0
burnIn = 20
nbrhood = nbrMask(4);

# M = D[1];
# N = D[2];
# crazy.anim(oim, M, N, nbrhood, iter, alpha, beta, sigma, burnIn);

## random image
M = D[1];
N = D[2];
dim = ising.posterior(oim, M, N, nbrhood, iter, alpha, beta, sigma, burnIn);
display.image(dim);
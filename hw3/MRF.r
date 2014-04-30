require(png)
require(R.basic)
require(package=animation)

## Reads a PNG image (as a simple matrix), which for some reason needs to be
## transposed and flipped to be displayed correctly
read.image = function(filename)
{
  y = readPNG(filename)
  
  ## This is just how I encoded pixel intensities
  ## After this, the labels will be black = -1 and white = +1
  y = y * 20 - 10
  
  n = nrow(y)
  
  ## transpose and flip image
  t(y)[,n:1]
}

## Displays an image
display.image = function(x, col=gray(seq(0,1,1/256)))
{
  w = dim(x)[1]
  h = dim(x)[2]
  par(mai = c(0,0,0,0))
  image(x, asp=h/w, col=col)
}

## Generate neighbor mask
nbrMask = function(n)
{      
  x = matrix(0,3,3);
  x[2,c(1,3)] = 1;
  x[c(1,3),2] = 1; 
  return(x)
}

## Get neighbors sum
getNbrValues = function(im, M, N, nbrMask)
{
  imPad = matrix(0,M+2,N+2);
  imPad[2:(M+1),2:(N+1)] = im;
  
  imPad[1,] = imPad[M+1,]; 
  imPad[M+2,] = imPad[2,];
  
  imPad[,1] = imPad[,N+1];
  imPad[,N+2] = imPad[,2];
  
  nbrs = apply2d(imPad, nbrMask, FUN=sum);
  #nbrs = matrix(nbrs, M+2, N+2);
  nbrs = nbrs[2:(M+1), 2:(N+1)];
  return(nbrs);
}

## Ising model(with no posterior)
ising.prior = function(im, M, N, nbrMask, iterations, al, be)
{
  flip = matrix(c(1,0),M,N);
  for(i in 1:iterations)
  {
    nbrhood = getNbrValues(im, M, N, nbrMask);
    U = -al*im - be*(im*nbrhood);
    p = exp(-U); 
    transp = (matrix(runif(M*N),M,N) > p)*flip*(-2)+1;
    im = im*transp;
    flip = 1-flip;;
    print(paste("Iteration ", i));
    #display.image(im);
  }
  return(im);
}

ising.posterior = function(trueIm, M, N, nbrMask, iterations, al, be, sigma, burnIn)
{
  res = matrix(0.0,M,N);
  im = (matrix(runif(M*N),M,N) > 0.5)*2 -1;
  flip = matrix(c(1,0), M,N,byrow=T);
  vec_of_sigma <- numeric();
  for(i in 1:iterations)
  {
    nbrhood = getNbrValues(im,M,N,nbrMask);
    U = -al*im-be*(im*nbrhood) + 0.5/(sigma^2)*(im-trueIm)^2;
    p = exp(-U);
    transp = (matrix(runif(M*N),M,N) > p)*flip*(-2)+1;
    im = im*transp;
    flip = 1 - flip;
    print(paste("Iteration ", i));
    if(i > burnIn) { res = res + im; }
    sigma=sqrt(sum((im-trueIm)^2)/(M*N));
    print(sigma);
    vec_of_sigma = c(vec_of_sigma, sigma);
  #  display.image(res);
  }
  plot(1:iterations,vec_of_sigma,type='o',lwd=1,
       ylim=c(1.2,1.71),
       main=expression(bold(paste(sigma, " for each iteration"))),
       ylab=expression(paste(sigma)),
       xlab=paste("iterations"));
  abline(h=1.217,col='black',lwd=3,lty=3);
  return(res/(iterations-burnIn));
}

crazy.anim = function(oim, M, N, nbrhood, iter, al,be,sig, burnIn = 0 )
{
  ## Now we'll animate the Hamiltonian integration
  ani.options(interval = 0.1)
  ani.start()
  #random_binary_image(al,be,sig);
  #ising.prior(oim,M,N, nbrhood, iter, al, be );
  ising.posterior(oim,M,N,nbrhood, iter,al,be,sig,burnIn);
  ani.stop()
}
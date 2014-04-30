source("/home/raj/Documents/Probabilistic_Modelling/Assignment4/code/MRF.r")
require(animation)

energy_func=function(M,sum_product,alpha,beta,y,nc,nr,sigma)
{
  #E=-(alpha*sum(M))-(beta*sum_product);
  sum=0;
  for(i in 1:nr)
  {
    for(j in 1:nc)
    {
      sum=sum+(M[i,j]-y[i,j])^2;
    }
  }
  E=-(alpha*sum(M))-(beta*sum_product)+((1/(2*sigma^2))*sum);
  return(E);
}
neighbor_addition=function(M,x,y,xplus,yplus,xmin,ymin,nr,nc)
{
  product=0;
  #1
  if(xmin==0 && ymin==0)
  {
    #print("IN CASE 1");
    product=M[y,x]*M[yplus,x]+M[y,x]*M[y,xplus];
  }
  #2
  else if(xplus>nc && ymin<1)
  {
    #print("IN CASE 2");
    product=M[y,x]*M[yplus,x]+M[y,x]*M[y,xmin];
  }
  #3
  else if(xplus>nc && yplus>nr)
  {
    #print("IN CASE 3");
    product=M[y,x]*M[ymin,x]+M[y,x]*M[y,xmin];
    #print(product);
  }
  #4
  else if(xmin==0 && yplus>nr)
  {
    #print("IN CASE 4");
    product=M[y,x]*M[ymin,x]+M[y,x]*M[y,xplus];
  }
  #5
  else if(xmin==0)
  {
    #print("IN CASE 5");
    product=M[y,x]*M[ymin,x]+M[y,x]*M[y,xplus]+M[y,x]*M[yplus,x];
  }
  #6
  else if(xplus>nc)
  {
    #print("IN CASE 6");
    product=M[y,x]*M[ymin,x]+M[y,x]*M[y,xmin]+M[y,x]*M[yplus,x];
  }
  #8
  else if(ymin==0)
  {
    #print("IN CASE 8");
    #print(M);
    #print(M[x,y]*M[xplus,y]);
    product=M[y,x]*M[y,xmin]+M[y,x]*M[y,xplus]+M[y,x]*M[yplus,x];
  }
  #9
  else if(yplus>nr)
  {
    #print("IN CASE 9");
    product=M[y,x]*M[y,xmin]+M[y,x]*M[y,xplus]+M[y,x]*M[ymin,x];
  }
  #7
  else
  {
    #print("IN CASE 7");
    product=M[y,x]*M[ymin,x]+M[y,x]*M[yplus,x]+M[y,x]*M[y,xmin]+M[y,x]*M[y,xplus];
  }
  return(product);
}
random_binary_image = function(al,be,sig)
{
  no_of_iter=50;
  orig=read.image('/home/raj/Documents/Probabilistic_Modelling/Assignment4/code/noisy-message.png');
  nc=dim(orig)[2];
  nr=dim(orig)[1];
  #nc=5;
  #nr=3;
  M = matrix(sample(c(-1, 1), nr * nc, replace = TRUE), nrow = nr)
  #print(M[1,1])
  #display.image(M)
  num_rows=nr;
  num_cols=nc;
  alpha=al;
  beta=be;
  sigma=sig;
  U=0;
  prev_U=0;
  rel_prob=0;
  delta_U=0;
  y1=M;
  M=orig;
 
  for(n in 1:no_of_iter)
  {
    print("ITERARTION");
    print(n);
    prev_U=U;
    sum_product=0;
    #print(M);
    if(n>1)
    {
      a<-sample(1:nr, 1, replace=T);
      #print(a);
      b<-sample(1:nc, 1, replace=T);
      #print(b);
      #print(M);
     
      M[a,b]<--M[a,b];
     
    }
    for(i in 1:nr)
    {
      for(j in 1:nc)
      {
        #### Calculate the neighbors of a pixel
        x=j;
        #print("X:");
        #print(x);
        y=i;
        #print("Y:");
        #print(y);
        yplus=(y+1);
        #print("Y_PLUS:");
        #print(yplus);
        xplus=(x+1);
        #print("X_PLUS:");
        #print(xplus);
        ymin=(y-1);
        #print("Y_MIN:");
        #print(ymin);
        xmin=(x-1);
        #print("X_MIN:");
        #print(xmin);
        product=neighbor_addition(M,x,y,xplus,yplus,xmin,ymin,nr,nc);
        #print("PRODUCT");
        #print(product);
        sum_product=sum_product+product;
        #print("SUM-PRODUCT");
        #print(sum_product);
        #print("-------------------------------------------")
      }
    }
    U=energy_func(y1,sum_product,alpha,beta,M,nc,nr,sigma);
    #print("U:");
    #print(U);
    if(n>1)
    {
      if(U<=prev_U)
      {
        delta_U=prev_U-U;
        
        M[a,b]<-M[a,b];
       
      }
      else
      {
        delta_U=U-prev_U;
        rel_prob=exp(-delta_U*beta); 
        num<-runif(1);
        if(rel_prob>=num)
        {
          M[a,b]<-M[a,b];
        }
        else
        {
          M[a,b]<--M[a,b];
        }
      }
    }
    display.image(M);
    
  }
}

crazy.anim = function(al,be,sig)
{
  ## Now we'll animate the Hamiltonian integration
  ani.options(interval = 0.1)
  ani.start()
  random_binary_image(al,be,sig);
  ani.stop()
}
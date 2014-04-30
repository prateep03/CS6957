rm(list = ls())

#options(digits=22)
library(package='LaplacesDemon')
## Read in the cross-sectional OASIS data, including hippocampal volumes
x = read.csv("oasis_cross-sectional.csv")
y = read.csv("oasis_hippocampus.csv")
hippo = merge(x, y, by = "ID")

## Let's look at only elderly subjects
hippo = hippo[hippo$Age >= 60,]
rownames(hippo) = NULL

## Create three categories of subjects according to dementia level (from CDR)
## cntrls: CDR = 0, Mild Dementia: CDR = 0.5, Dementia: CDR >= 1.0
hippo$Group[hippo$CDR == 0.0] = "cntrl"
hippo$Group[hippo$CDR == 0.5] = "Mild"
hippo$Group[hippo$CDR >= 1.0] = "Dementia"
hippo$Group = factor(hippo$Group, levels=c("cntrl", "Mild", "Dementia"))

group1<-subset(hippo,hippo$Group == "cntrl",select=c(RightHippoVol));
group2<-subset(hippo,hippo$Group == "Mild",select=c(RightHippoVol));
group3<-subset(hippo,hippo$Group == "Dementia",select=c(RightHippoVol));

mean_group1<-mean(group1$RightHippoVol);
sd_group1<-var(group1$RightHippoVol);

mean_group2<-mean(group2$RightHippoVol);
sd_group2<-var(group2$RightHippoVol);

mean_group3<-mean(group3$RightHippoVol);
sd_group3<-var(group3$RightHippoVol);

# Make the distribution
s = seq(0.1,8, 0.1)

## Density function for Inverse Gamma:
dinvgamma = function(x, alpha, beta)
{
  #(beta^alpha / gamma(alpha)) * x^(-alpha - 1) * exp(-beta / x); # actual form
  alpha * log(beta) - log(gamma(alpha)) +(-alpha-1) * log(x) - (beta/x); # log form
}

## Density function for Normal-Inverse-Gamma
## with parameters mu0, lambda, alpha, beta
dninvg = function(s, sd, mu0, lambda, alpha, beta)
{
  return(dnorm(s, mean = mu0, sd = sqrt(sd/lambda)) * dinvgamma(sd, alpha, beta));
}

joint_posterior<-function(s,mean,sd, mu0, lambda, alpha, beta)
{
	return(dnorm(s, mean = mean, sd = sd)*dninvg(s, s, mu0, lambda, alpha, beta));
}

## Density function for Students-t distribution
## with parameters v,s
dstudentst = function(t, v, beta, n)
{
  s = sqrt(2*beta/n)
(gamma((v+1)/2)/(gamma(v/2))) * 1/(sqrt(pi)*s*v) * ((1+t^2/(v*s))^(-(v+1)/2));
#sqrt((n*beta) / (2*pi)) * (gamma((v+1)/2)/(gamma(v/2))) * ((1+t^2/(v*s))^(-(v-1)/2));

#   -1/2 * log(pi) - log(s) + log(gamma((v+1)/2)) - log(gamma(v/2)) 
#       -((v+1)/2) * log( 1 + t^2/(v*s));
}

#pdf('4_plots.pdf')

## Initial parameters
mu0 = 0;
n0 = 10^-6;
alpha0 = 10^-6;
beta0 = 10^-6;
# n
ngroup.cntrl = nrow(group1);  
ngroup.mild = nrow(group2);  
ngroup.dem = nrow(group3);
# alphan (shape)
alphan.cntrl = (alpha0+ngroup.cntrl/2); 
alphan.mild = (alpha0+ngroup.mild/2); 
alphan.dem = (alpha0+ngroup.dem/2);
# betan (scale)
betan.cntrl = beta0+sum((group1-mean_group1)^2)/2;
betan.mild = beta0+sum((group2-mean_group2)^2)/2;
betan.dem = beta0+sum((group3-mean_group3)^2)/2;
# rate = 1/scale
rate.cntrl = 1/betan.cntrl;
rate.mild = 1/betan.mild;
rate.dem = 1/betan.dem;
# Sample mean & var for each grp
sample_mean.cntrl     = (1/ngroup.cntrl)*(sum(group1));
sample_variance.cntrl = (1/(ngroup.cntrl-1))*(sum((group1-sample_mean.cntrl)^2));
sample_mean.mild      = (1/ngroup.mild)*(sum(group2));
sample_variance.mild  = (1/(ngroup.mild-1))*(sum((group2-sample_mean.mild)^2));
sample_mean.dem       = (1/ngroup.dem)*(sum(group3));
sample_variance.dem   = (1/(ngroup.dem-1))*(sum((group3-sample_mean.dem)^2));

## 4a
# Joint posterior
joint_posterior.cntrl = joint_posterior(s,mean_group1,sqrt(sd_group1),mu0,n0,alpha0,beta0); # cntrl
joint_posterior.mild  = joint_posterior(s,mean_group2,sqrt(sd_group2),mu0,n0,alpha0,beta0); # Mild
joint_posterior.dem   = joint_posterior(s,mean_group3,sqrt(sd_group3),mu0,n0,alpha0,beta0); # Dementia

## 4a
## Marginal posterior of \sigma is IG
samples<-seq(from=-10,to=100,length.out=10^6)
final_samples=rinvgamma(10^6, shape = alphan.cntrl, scale = betan.cntrl);
hist(final_samples,freq=FALSE,
     main=expression(bold(paste("Marginal posterior for ", sigma^2, " Group 1"))))
lines(density(final_samples), col = "red", lwd = 1);
abline(v = (sample_variance.cntrl), col='blue', lwd = 3, lty = 2)

## 4b
## Marginal posterior of \mu is Students-t
samples = seq(from=-10^4, to=10^4, length.out=10^6);

xbar_cntrl = (1/ngroup.cntrl) * sum(group1);
mun_cntrl = (n0*mu0 + xbar_cntrl*ngroup.cntrl) / (n0+ngroup.cntrl);
nn_cntrl = n0 + ngroup.cntrl;
betan.cntrl = betan.cntrl + ((n0*ngroup.cntrl)/(n0+ngroup.cntrl)) * 1/2 * (xbar_cntrl-mu0)^2;
ncp_cntrl = sqrt(nn_cntrl*alphan.cntrl / betan.cntrl) * (samples - mun_cntrl);
samples.mean1 = dstudentst(samples - mun_cntrl, 2*alphan.cntrl, 
                           betan.cntrl,nn_cntrl/50);
#samples.mean1 = dt(samples, 2*alphan.cntrl, ncp_cntrl, log=TRUE);
# samples.mean1 = exp(dst(samples, mu = mun_cntrl, sigma = sqrt(2*betan.cntrl/nn_cntrl), nu = 2*alphan.cntrl, log=TRUE));

xbar_mild = (1/ngroup.mild) * sum(group2);
mun_mild = (n0*mu0 + xbar_mild*ngroup.mild) / (n0+ngroup.mild);
nn_mild = n0 + ngroup.mild;
betan.mild = betan.mild + ((n0*ngroup.mild)/(n0+ngroup.mild)) * 1/2 * (xbar_mild-mu0)^2;
ncp_mild = sqrt(nn_mild*alphan.mild / betan.mild) * (samples - mun_mild);
samples.mean2 = dstudentst(samples - mun_mild, 2*alphan.mild, 
                           betan.mild,nn_mild/50);
#samples.mean2 = dt(samples, 2*alphan.mild, ncp_mild, log=TRUE);
#samples.mean2 = exp(dst(samples, mu = mun_mild, sigma = sqrt(2*betan.mild/nn_mild), nu = 2*alphan.mild, log=TRUE));

xbar_dem = (1/ngroup.dem) * sum(group3);
mun_dem = (n0*mu0 + xbar_dem*ngroup.dem) / (n0+ngroup.dem);
nn_dem = n0 + ngroup.dem;
betan.dem = betan.dem + ((n0*ngroup.dem)/(n0+ngroup.dem)) * 1/2 * (xbar_dem-mu0)^2;
ncp_dem = sqrt(nn_dem*alphan.dem / betan.dem) * (samples - mun_dem);
samples.mean3 = dstudentst(samples - mun_dem, 2*alphan.dem, 
                          betan.dem,nn_dem/50);
#samples.mean3 = dt(samples, 2*alphan.dem, ncp_dem, log=TRUE);
#samples.mean3 = exp(dst(samples, mu = mun_dem, sigma = sqrt(2*betan.dem/nn_dem), nu = 2*alphan.dem, log=TRUE));

## Plot them all
plot(samples, samples.mean1, type='l', col='red', lwd=1,
     ylim=c(0,2e-5),xlim=c(2800,3800),
     #ylim=c(0,max(samples.mean1, samples.mean2, samples.mean3)),xlim=c(2800,3800),
     main = expression(bold(paste("Marginal Posteriors for ",mu))),
     ylab = expression(paste("p(", mu, "|",y_ij, ")")));
 lines(samples, samples.mean2, col='blue', lwd=1);
 lines(samples, samples.mean3, col='green', lwd=1);

abline(v = (sample_mean.cntrl), col='red', lwd = 3, lty = 2)
abline(v = (sample_mean.mild), col='blue', lwd = 3, lty = 2)
abline(v = (sample_mean.dem), col='green', lwd = 3, lty = 2)
# legend('topright', c("Control", "Mild", "Dementia"),
#        col = c("red", "blue", "green"), lwd = 3, 
#        xjust = 0, yjust = 1)

## 4c
# d_ij is a normal distribution.
n = 10^6;
samples_sigma.control = 1/rgamma(n, shape = alphan.cntrl, scale = 1/betan.cntrl)
samples_sigma.mild = 1/rgamma(n, shape = alphan.mild, scale = 1/betan.mild)
samples_sigma.dem = 1/rgamma(n, shape = alphan.dem, scale = 1/betan.dem)

samples.sd = sqrt(samples_sigma.control/ngroup.cntrl + samples_sigma.mild/ngroup.mild)
samples.d12 = rnorm(n, mean = mean_group1 - mean_group2, sd = samples.sd);
samples.sd = sqrt(samples_sigma.control/ngroup.cntrl + samples_sigma.dem/ngroup.dem)
samples.d13 = rnorm(n, mean = mean_group1 - mean_group3, sd = samples.sd);
samples.sd = sqrt(samples_sigma.mild/ngroup.mild + samples_sigma.dem/ngroup.dem)
samples.d23 = rnorm(n, mean = mean_group2 - mean_group3, sd = samples.sd);

#hist(samples_d23, freq=FALSE);

### 4d
n = 100;
samples_mean.dem = rt(n, 2*alphan.dem, 
                         sqrt((nn_dem*alphan.dem) / betan.dem) * (mean_group3 - mun_dem) ); 
#samples_mean.dem = rst(n, mu = mun_cntrl, sigma = sqrt(2*betan.dem/nn_dem), nu = 2*alphan.dem)
samples_mean.mild = rt(n, 2*alphan.mild, 
                         sqrt((nn_mild*alphan.mild) / betan.mild) * (mean_group2 - mun_mild) ); 
#samples_mean.mild = rst(n, mu = mun_mild, sigma = sqrt(2*betan.mild/nn_mild), nu = 2*alphan.mild)
samples_mean.cntrl = rt(n, 2*alphan.cntrl, 
                       sqrt((nn_cntrl*alphan.cntrl) / betan.cntrl) * (mean_group1 - mun_cntrl) ); 
#samples_mean.cntrl = rst(n, mu = mun_cntrl, sigma = sqrt(2*betan.cntrl/nn_cntrl), nu = 2*alphan.cntrl)

sa.23 = samples_mean.mild - samples_mean.dem ;
sa.13 = samples_mean.cntrl - samples_mean.dem ;
sa.12 = samples_mean.cntrl - samples_mean.mild ;

p.23 = t.test(sa.23, alternative="less");
p.13 = t.test(sa.13, alternative="less");
p.12 = t.test(sa.12, alternative="less");

cat(paste("P(d_12 < 0) =", (sum(samples.d12 < 0) / length(samples.d12)), ", p-value(12) = ", p.12$p.value, "\n"));
cat(paste("P(d_13 < 0) =", (sum(samples.d13 < 0) / length(samples.d13)), ", p-value(13) = ", p.13$p.value, "\n"));
cat(paste("P(d_23 < 0) =", (sum(samples.d23 < 0) / length(samples.d23)), ", p-value(23) = ", p.23$p.value, "\n"));

#dev.off()
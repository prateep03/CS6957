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

#pdf('5_plots.pdf')

## Initial parameters
mu0 = 0;
n0 = 10^-6;
# n
ngroup.cntrl = nrow(group1);  
ngroup.mild = nrow(group2);  
ngroup.dem = nrow(group3);
# Sample mean & var for each grp
sample_mean.cntrl     = (1/ngroup.cntrl)*(sum(group1));
sample_variance.cntrl = (1/(ngroup.cntrl-1))*(sum((group1-sample_mean.cntrl)^2));
sample_mean.mild      = (1/ngroup.mild)*(sum(group2));
sample_variance.mild  = (1/(ngroup.mild-1))*(sum((group2-sample_mean.mild)^2));
sample_mean.dem       = (1/ngroup.dem)*(sum(group3));
sample_variance.dem   = (1/(ngroup.dem-1))*(sum((group3-sample_mean.dem)^2));

# 
# joint_posterior<-function(s,mean,sd)
# {
#   return (dnorm(s, mean = mean, sd = sd));
# }

## 5a
# Make the distribution
samples = seq(from=-10, to=100, length.out=10^6)
final_samples = rinvchisq(10^6, ngroup.cntrl-1, sample_variance.cntrl);
hist(final_samples, freq=FALSE,
     main=expression(bold(paste("Marginal posterior for ", sigma^2, " Group 1"))));
lines(density(final_samples), col="red", lwd=1);
abline(v = sample_variance.cntrl, col = "blue", lwd=3, lty=2)

## 5b
## Marginal posterior of \mu is Students-t
samples = seq(from=-10^5, to=10^5, length.out=10^6);

samples.mean1 = dst(samples, mu = sample_mean.cntrl, sigma = sample_variance.cntrl/ngroup.cntrl, 
                    nu = ngroup.cntrl-1);

samples.mean2 = dst(samples, mu = sample_mean.mild, sigma = sample_variance.mild/ngroup.mild, 
                    nu = ngroup.mild-1);

samples.mean3 = dst(samples, mu = sample_mean.dem, sigma = sample_variance.dem/ngroup.dem, 
                    nu = ngroup.dem-1);

## Plot them all
plot(samples, samples.mean1, type='l', col='red', lwd=1,
     #ylim=c(0,0.00020),xlim=c(2800,3800),
     ylim=c(0,max(samples.mean1, samples.mean2, samples.mean3)), xlim=c(-3e4, 3e4),
     main = expression(bold(paste("Marginal Posteriors for ",mu))),
     ylab = expression(paste("p(", mu, "|",y_ij, ")")));
lines(samples, samples.mean2, col='blue', lwd=1);
lines(samples, samples.mean3, col='green', lwd=1);

abline(v = (sample_mean.cntrl), col='red', lwd = 3, lty = 2)
abline(v = (sample_mean.mild), col='blue', lwd = 3, lty = 2)
abline(v = (sample_mean.dem), col='green', lwd = 3, lty = 2)
# legend('topleft', c("Control", "Mild", "Dementia"),
#        col = c("red", "blue", "green"), lwd = 3, 
#        xjust = 0, yjust = 1)

## 5c
n = 10^6;
samples_sigma.control = rinvchisq(n, ngroup.cntrl, sample_variance.cntrl);
samples_sigma.mild = rinvchisq(n, ngroup.mild, sample_variance.mild);
samples_sigma.dem = rinvchisq(n, ngroup.dem, sample_variance.dem);

samples_sd = sqrt(samples_sigma.control/ngroup.cntrl + samples_sigma.mild/ngroup.mild)
samples_d12 = rnorm(n, mean = mean_group1 - mean_group2, sd = samples_sd);
samples_sd = sqrt(samples_sigma.control/ngroup.cntrl + samples_sigma.dem/ngroup.dem)
samples_d13 = rnorm(n, mean = mean_group1 - mean_group3, sd = samples_sd);
samples_sd = sqrt(samples_sigma.mild/ngroup.mild + samples_sigma.dem/ngroup.dem)
samples_d23 = rnorm(n, mean = mean_group2 - mean_group3, sd = samples_sd);

cat(paste("P(d_12 < 0) =", (sum(samples_d12 < 0) / length(samples_d12)), "\n"));
cat(paste("P(d_13 < 0) =", (sum(samples_d13 < 0) / length(samples_d13)), "\n"));
cat(paste("P(d_23 < 0) =", (sum(samples_d23 < 0) / length(samples_d23)), "\n"));

#dev.off();
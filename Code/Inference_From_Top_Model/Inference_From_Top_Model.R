# Clear the workspace
rm(list=ls())

################################################################
### Script to fit Top Model and make inference for manuscript
################################################################

# Load necessary libraries
library(coda)        # For trace plots and density plots
library(truncnorm)   # For truncnorm function
library(MASS)        # For mvrnorm function
library(Matrix)      # For Kronecker function
library(mvtnorm)     # For multivariate normal distribution
library(R.utils)     # For saving files

# Define your data and parameters
n.iter <- 20000  # Number of MCMC iterations
file_name <- "Top_Model_Output"  # Base name for saving files
save_dir <- "~/Dropbox/GitHub/Mercury-Waterfowl-RiskAssessment/R_Output"  # Directory to save output

# Load data from GitHub
data <- read.csv(paste0("https://github.com/perrywilliamsunr/",
                        "Mercury-Waterfowl-RiskAssessment/raw/main/Data/data.csv"))

# Prepare data for MCMC
n <- length(unique(data$bird_id))
ind <- seq(1, n * 3, 3)
y <- c(data$y1[ind[1]:(ind[1] + 2)], data$y2[ind[1]:(ind[1] + 2)])
for(i in 2:length(ind)){
  y.tmp <- c(data$y1[ind[i]:(ind[i] + 2)], data$y2[ind[i]:(ind[i] + 2)])
  y <- c(y, y.tmp)
}
y.untransformed <- y
y <- log(y)  # Log-transform the data

# For reproducibility
set.seed(123)

# Matrix building materials
m.tmp1 <- matrix(0, 6, 6)
m.tmp1[c(1, 8, 15)] <- 1
m.tmp2 <- matrix(0, 6, 6)
m.tmp2[c(22, 29, 36)] <- 1
m.tmp3 <- matrix(0, 6, 6)
m.tmp3[c(2, 3, 7, 9, 13, 14)] <- 1
m.tmp4 <- matrix(0, 6, 6)
m.tmp4[c(23, 24, 28, 30, 34, 35)] <- 1
m.tmp5 <- matrix(0, 6, 6)
m.tmp5[c(4, 5, 6, 10, 11, 12, 16, 17, 18)] <- 1
m.tmp5[c(19, 20, 21, 25, 26, 27, 31, 32, 33)] <- 1

# Independent variables
n.mall <- 5
liver <- rep(rep(1:0, each = 3), n)
indices <- seq(1, nrow(data), by = 3)
species <- rep(data$species, each = 2)
flank.tmp <- (data$flank1 + data$flank2) / 2
flank.scaled <- (flank.tmp - 8821.888) / 9209.355  # Scaling flank
axillary.tmp <- (data$axillary1 + data$axillary2) / 2
axillary.scaled <- (axillary.tmp - mean(unique(axillary.tmp))) / sd(unique(axillary.tmp))  # Scaling axillary
breast.tmp <- (data$breast1 + data$breast2) / 2
breast.scaled <- (breast.tmp - mean(unique(breast.tmp))) / sd(unique(breast.tmp))  # Scaling breast
flank <- rep(flank.scaled, each = 2)
axillary <- rep(axillary.scaled, each = 2)
breast <- rep(breast.scaled, each = 2)

# Model selection variables
model.number <- 35
Species <- TRUE
Liver <- TRUE
Breast <- FALSE
Flank <- TRUE
Species.Liver <- TRUE
Species.Flank <- TRUE
Species.Axillary <- FALSE
Species.Breast <- FALSE

# Construct design matrix X
X <- matrix(1, length(breast), 1)
if(Liver) X <- cbind(X, liver)
if(Flank) X <- cbind(X, flank)
if(Breast) X <- cbind(X, breast)
if(Species) X <- cbind(X, species)
if(Species.Liver) X <- cbind(X, species * liver)
if(Species.Flank) X <- cbind(X, species * flank)
if(Species.Axillary) X <- cbind(X, (1 - species) * axillary)
if(Species.Breast) X <- cbind(X, species * breast)

# Set initial values, tuning parameters, and prior parameters
initial_values <- list(
  s.l = 1,
  s.b = 1,
  rho = c(0.5, 0.5, 0.5),
  beta = rep(0, ncol(X))
)
tuning_parameters <- list(
  s.l.tune = 686.3144,
  s.b.tune = 0.005948354,
  rho.tune = rep(0.2086,3)
)
prior_parameters <- list(
  #s.l.prior = c(825, 100^2),
  s.l.prior=c(0.01,0.01),
  s.b.prior = c(.01, .01),
  rho.prior = c(100, 1),
  mu.beta = matrix(0, ncol(X), 1),
  s2.beta = 1 * diag(ncol(X))
)

# Source the MCMC_Algorithm.R script from GitHub
source(paste0('https://raw.githubusercontent.com/perrywilliamsunr/',
              'Mercury-Waterfowl-RiskAssessment/main/Code/MCMC_Algorithm/',
              'MCMC_Algorithm.R'))

# Run the MCMC sampling
run_mcmc(y, X, n.iter, initial_values, tuning_parameters, prior_parameters, file_name, save_dir)

# After running this script, you should find output files in the specified directory

# Load the MCMC output file
load(file.path(save_dir, paste0(file_name, ".iter", n.iter, ".RData")))
burn <- floor(n.iter / 2)
thin <- 1
ind <- seq(burn, n.iter, thin)
length(ind)

beta.save=MCMC.Output$beta
rho.save=MCMC.Output$rho
s.l.save=MCMC.Output$s.l
s.b.save=MCMC.Output$s.b
Dbar.save=MCMC.Output$Dbar

par(mfrow=c(2,2))
for(i in 1:ncol(rho.save)){
  plot(rho.save[ind,i],type='l')
}
# plot(rho.tune.save[ind,1],type='l')
# rho.tune
q=apply(rho.save[ind,],2,quantile,c(0.025,0.975))

par(mfrow=c(3,2),mar=c(4,4,1,1))
for(i in 1:ncol(beta.save)){
  plot(beta.save[ind,i],type='l',ylim=c(-10,10))
  abline(h=0,col=2)
}
q=apply(beta.save,2,quantile,c(0.025,0.975))
q

par(mfrow=c(2,1))
plot(s.l.save[ind],type='l')

# plot(s.l.tune.save,type='l')
# s.l.tune

par(mfrow=c(2,1))
plot(s.b.save[ind],type='l')
# plot(s.b.tune.save[ind],type='l')
# s.b.tune



# ###
# ### Model checking 
# ###
# 
# Bayes.count1=matrix(0,n.iter,length(y))
# Bayes.count2=matrix(0,n.iter,1)
# 
# for(k in ind){
#   mu=matrix(X%*%beta.save[k,],n*6,1)
#   block=s.l.save[k]^2*m.tmp1+
#     s.b.save[k]^2*m.tmp2+
#     rho.save[k,1]*m.tmp3*s.l.save[k]^2+
#     rho.save[k,2]*m.tmp4*s.b.save[k]^2+
#     rho.save[k,3]*s.l.save[k]*s.b.save[k]*m.tmp5
#   Sigma=kronecker(diag(n),block)
#   Bayes.count1[k,]=ifelse((log(y.ppd.save[k,])-mu)^2/mu<
#                             ((y-mu)^2/mu),1,0)
#   if(-2*dmvnorm(log(y.ppd.save[k,]),mu,Sigma,log=TRUE)<
#      -2*dmvnorm(y,mu,Sigma,log=TRUE)){
#     Bayes.count2[k,1]=1
#   }
# }
# 
# ## (Bayes.p=apply(Bayes.count1[ind,],2,mean))
# (Bayes.p=mean(Bayes.count2[ind,]))


###
### DIC
###

calculate_means <- function(mat) {
  if (is.vector(mat)) {
    return(mean(mat))
  } else if (ncol(mat) == 1) {
    return(mean(mat))
  } else {
    return(apply(mat, 2, mean))
  }
}


post.beta.mean=calculate_means(beta.save[ind,])
post.rho.mean=apply(rho.save[ind,],2,mean)
post.sb.mean=mean(s.b.save[ind])
post.sl.mean=mean(s.l.save[ind])
post.block.mean=post.sl.mean^2*m.tmp1+
  post.sb.mean^2*m.tmp2+
  post.rho.mean[1]*m.tmp3*post.sl.mean^2+
  post.rho.mean[2]*m.tmp4*post.sb.mean^2+
  post.rho.mean[3]*post.sl.mean*post.sb.mean*m.tmp5
post.Sigma.mean=kronecker(diag(n),post.block.mean)
y.sim=rmvnorm(1,c(X%*%post.beta.mean),post.Sigma.mean)
plot(y,y.sim,xlim=c(2,10),ylim=c(2,10))
abline(0,1)


## Deviance at mean
Dhat=-2*(dmvnorm(c(y),c(X%*%post.beta.mean),post.Sigma.mean,log=TRUE))

## Mean deviance
Dbar=mean(Dbar.save[ind]) 

pD=Dbar-Dhat
DIC=Dhat+2*pD

round(Bayes.p,2)
round(DIC)
round(Dbar)
round(Dhat)
round(pD)
model.number



# ## Equivalently
# #(DIC = Dbar+pD)
# 

###
### Results
###

data=read.csv(paste0("~/Dropbox/projects/Waterfowl/Ducks/",
                     "Wood Ducks Mercury/Manuscript/Appendix/prediction_data.csv"))


###
### Emperical Means
###

## Wood duck breast
mean(unlist(data[data$species==0,4:6]))
range(unlist(data[data$species==0,4:6]))
sd(unlist(data[data$species==0,4:6]))

## Wood duck liver
mean(unlist(data[data$species==0,1:3]))
range(unlist(data[data$species==0,1:3]))
sd(unlist(data[data$species==0,1:3]))

## Wood duck Flank Feathers
mean(unlist(data[data$species==0,7:8]))
range(unlist(data[data$species==0,7:8]))
sd(unlist(data[data$species==0,7:8]))

## Mallard breast
mean(unlist(data[data$species==1,4:6]))
range(unlist(data[data$species==1,4:6]))
sd(unlist(data[data$species==1,4:6]))

## Mallard liver
mean(unlist(data[data$species==1,1:3]))
range(unlist(data[data$species==1,1:3]))
sd(unlist(data[data$species==1,1:3]))

## Mallard Flank Feathers
mean(unlist(data[data$species==1,7:8]))
range(unlist(data[data$species==1,7:8]))
sd(unlist(data[data$species==1,7:8]))


## Wood duck Breast tissue exceeding standards
sum(apply(data[data$species==0,4:6],1,max)>300)
sum(apply(data[data$species==0,4:6],1,max)>1000)

## Mallard Breast tissue exceeding standards
sum(apply(data[data$species==1,4:6],1,max)>300)
sum(apply(data[data$species==1,4:6],1,max)>1000)

## Wood duck Liver tissue exceeding standards
sum(apply(data[data$species==0,1:3],1,max)>300)
sum(apply(data[data$species==0,1:3],1,max)>1000)

## Mallard Liver tissue exceeding standards
sum(apply(data[data$species==1,1:3],1,max)>300)
sum(apply(data[data$species==1,1:3],1,max)>1000)

###############################################################
###. Marginal posterior distributions
################################################################



sd.flank

#flank breast slope in wood ducks
quantile(beta.save[ind,3],c(0.025,0.5,0.975))
plot(density(beta.save[ind,3]))

# flank breast slope in mallards
quantile(beta.save[ind,3]+beta.save[ind,6],c(0.025,0.5,0.975))

# flank liver slope in wood ducks
quantile(beta.save[ind,3] ,c(0.025,0.5,0.975))

# flank liver slope in mallards
quantile(beta.save[ind,3]+beta.save[ind,6],c(0.025,0.5,0.975))

# liver intercept in wood ducks
quantile(beta.save[ind,1]+beta.save[ind,2] ,c(0.025,0.5,0.975))

# liver intercept in mallards
quantile(beta.save[ind,1]+beta.save[ind,2]+beta.save[ind,4]+beta.save[ind,5] ,c(0.025,0.5,0.975))

# mallard breast compared to wood duck breast
quantile(beta.save[ind,4] ,c(0.025,0.5,0.975))

# mallard liver compared to wood duck liver
quantile(beta.save[ind,4]+beta.save[ind,5] ,c(0.025,0.5,0.975))



# correlation among triplicates in liver
quantile(rho.save[ind,1],c(0.025,0.5,0.975))

# correlation among triplicates in breast
quantile(rho.save[ind,2],c(0.025,0.5,0.975))

# correlation among triplicates in a bird
quantile(rho.save[ind,3],c(0.025,0.5,0.975))



###############################################################
###. Predict liver and breast tissue mercury from Flank Feathers
################################################################

head(data)
min.flank=min(data[,7:8])
max.flank=max(data[,7:8])
mean.flank=mean((data$fl1+data$fl2)/2)
sd.flank=sd((data$fl1+data$fl2)/2)
flank.hg=seq(min.flank,max.flank,1000)
flank.hg.scale=(flank.hg-mean.flank)/sd.flank

flank.ind=1:length(flank.hg.scale)
species.ind=0:1
liver.ind=0:1

y.ppd.list=list()
log.y.ppd.list=list()

s=1
i=1

# outer loop for s
for(s in 1:2){
  # inner loop for i
  for(i in flank.ind){
    y.ppd=matrix(NA,nrow=length(ind),ncol=6)
    log.y.ppd=matrix(NA,nrow=length(ind),ncol=6)
    
    # innermost loop for k
    for(k in 1:length(ind)){
      
      X.pred=matrix(c(rep(1,6),
                      rep(species.ind[s],6), # species
                      rep(0:1,each=3),       # breast or liver
                      rep(flank.hg.scale[i],6), # scaled flank measurement
                      rep(species.ind[s],6)*rep(0:1,each=3), #Species*liver
                      rep(species.ind[s],6)*rep(flank.hg.scale[i],6)),6,6) # species*flank
      mu=X.pred%*%beta.save[ind[k],]
      
      post.Sigma.pred=s.l.save[ind[k]]^2*m.tmp1+
        s.b.save[ind[k]]^2*m.tmp2+
        rho.save[ind[k],1]*m.tmp3*s.l.save[ind[k]]^2+
        rho.save[ind[k],2]*m.tmp4*s.b.save[ind[k]]^2+
        rho.save[ind[k],3]*s.l.save[ind[k]]*s.b.save[ind[k]]*m.tmp5
      
      
      log.y.ppd[k,]=rmvnorm(1,mu,post.Sigma.pred)
      y.ppd[k,]=exp(log.y.ppd[k,])
    }
    
    # store the matrix in the list with a name indicating the comnbination of s and i
    y.ppd.list[[paste("s", s, "_i", i, sep = "")]] <- y.ppd
    log.y.ppd.list[[paste("s", s, "_i", i, sep = "")]] <- log.y.ppd
    
    print(c(s,i))
  }
}

##y.ppd.list[["s1_i2"]]

quantiles=apply(log.y.ppd.list[[1]][,c(1,4)],2,quantile,c(0.025,0.5,0.975))
## c(1,4) are the first of triplicate measurements for 
##      liver and breast, respectively

df <- data.frame(
  FlankHg = flank.hg.scale[1], 
  BreastLower = quantiles[1],
  BreastMedian = quantiles[2],
  BreastUpper = quantiles[3],
  LiverLower = quantiles[4],
  LiverMedian = quantiles[5],
  LiverUpper=quantiles[6]
)


for(i in 2:length(log.y.ppd.list)){
  quantiles=apply(log.y.ppd.list[[i]][,c(1,4)],2,quantile,c(0.025,0.5,0.975))
  df[i,]=c(flank.hg.scale[i], quantiles[1], quantiles[2],
           quantiles[3], quantiles[4], quantiles[5], quantiles[6])
}
dfb=rep(flank.hg,2)
df[,1]=dfb

###
### Plot (This is an alternative plot to Fig. 2 in the manuscript)
###

##
## A
##

par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(df[1:49,1],df[1:49,2],type='l',
     ylim=c(2,12),xlim=c(0,10000),
     ylab="log([tHg]) in breast",
     xlab="",
     main="Wood duck breast",
     lty=2
)
lines(df[1:49,1],df[1:49,3],lty=1)
lines(df[1:49,1],df[1:49,4],lty=2)
points((data$fl1[1:39]+data$fl2[1:39])/2,log(data$br1[1:39]),pch=16)
points((data$fl1[1:39]+data$fl2[1:39])/2,log(data$br2[1:39]),pch=16)
points((data$fl1[1:39]+data$fl2[1:39])/2,log(data$br3[1:39]),pch=16)

##
## B
##

plot(df[50:98,1],df[50:98,2],type='l',
     ylim=c(2,12),xlim=c(0,11000),
     ylab="",
     xlab="",
     main="Mallard breast",
     lty=2
)

lines(df[50:98,1],df[50:98,3],lty=1)
lines(df[50:98,1],df[50:98,4],lty=2)
points((data$fl1[data$species==1]+data$fl2[data$species==1])/2,log(data$br1[data$species==1]),pch=16)
points((data$fl1[data$species==1]+data$fl2[data$species==1])/2,log(data$br2[data$species==1]),pch=16)
points((data$fl1[data$species==1]+data$fl2[data$species==1])/2,log(data$br3[data$species==1]),pch=16)

##
## C
##

plot(df[1:49,1],df[1:49,5],type='l',
     ylim=c(2,12),xlim=c(0,10000),
     ylab="log([tHg]) in liver",
     xlab="[tHg] in wood duck flank feathers",
     main="Wood duck liver",
     lty=2
)
lines(df[1:49,1],df[1:49,6],lty=1)
lines(df[1:49,1],df[1:49,7],lty=2)
points((data$fl1[1:39]+data$fl2[1:39])/2,log(data$lv1[1:39]),pch=16)
points((data$fl1[1:39]+data$fl2[1:39])/2,log(data$lv2[1:39]),pch=16)
points((data$fl1[1:39]+data$fl2[1:39])/2,log(data$lv3[1:39]),pch=16)

##
## D
##

plot(df[50:98,1],df[50:98,5],type='l',
     ylim=c(2,12),xlim=c(0,11000),
     ylab="",
     xlab="[tHg] in mallard flank feathers",
     main="Mallard liver",lty=2
)

lines(df[50:98,1],df[50:98,6],lty=1)
lines(df[50:98,1],df[50:98,7],lty=2)
points((data$fl1[data$species==1]+data$fl2[data$species==1])/2,log(data$lv1[data$species==1]),pch=16)
points((data$fl1[data$species==1]+data$fl2[data$species==1])/2,log(data$lv2[data$species==1]),pch=16)
points((data$fl1[data$species==1]+data$fl2[data$species==1])/2,log(data$lv3[data$species==1]),pch=16)
legend(4000,5,c("Median","95% CRI","Observed Data"), col=1, lty=c(1,2,NA),pch=c(NA,NA,16))


##############################################################################################
### Tool to predict % chance breast tissue exceeds EPA guidelines using flank feather mercury
##############################################################################################

prob.greater.than.EPA.300=matrix(,49,2)
## Wood ducks
for(i in 1:49){
  prob.greater.than.EPA.300[i,1]=sum(y.ppd.list[[i]][,1:3]>300 )/length(y.ppd.list[[i]][,1:3])
  prob.greater.than.EPA.300[i,2]=sum(y.ppd.list[[i+49]][,1:3]>300 )/length(y.ppd.list[[i+49]][,1:3])
}
par(mfrow=c(2,1))
plot(flank.hg,prob.greater.than.EPA.300[,1],type='l',
     ylim=c(0,1),
     xlab="",
     ylab="Probability exceeds USEPA threshold",
     main="Wood ducks"
)
rug(c(data$fl1[data$species==0],data$fl2[data$species==0]))

plot(flank.hg,prob.greater.than.EPA.300[,2],type='l',
     ylim=c(0,1),
     xlab="Flank feather total mercury concentration (ng/g)",
     ylab="Probability exceeds USEPA threshold",
     main="Mallards"
)
rug(c(data$fl1[data$species==1],data$fl2[data$species==1]))























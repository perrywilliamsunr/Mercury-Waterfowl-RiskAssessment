rm(list=ls())

###
### Required packages
###

library(doParallel)
library(foreach)
library(invgamma)
library(Matrix)
library(mvtnorm)
library(truncnorm)

###
### Load data
###

data=read.csv(paste0("https://github.com/perrywilliamsunr/",
                     "Wood-Duck-Mercury-Contamination/raw/",
                     "main/Data/data.csv"))

###
### K folds
###

myCluster <- makeCluster(44, # number of cores to use
                         type = "FORK") # type of cluster
detectCores()

registerDoParallel(myCluster)

K=1:length(unique(data$bird_id))


foreach(l = K, .packages=c('invgamma', 'Matrix', 'mvtnorm', 'truncnorm')) %dopar% {
  
  n.iter=20000
  
  ##
  ## Develop training data and testing data
  ##
  
  data.train=data[data$bird_id!=l,]
  
  
  n=length(unique(data.train$bird_id))
  ind=seq(1,n*3,3)  
  y=c(data.train$y1[ind[1]:(ind[1]+2)],data.train$y2[ind[1]:(ind[1]+2)])
  for(i in 2:length(ind)){
    y.tmp=c(data.train$y1[ind[i]:(ind[i]+2)],data.train$y2[ind[i]:(ind[i]+2)])
    y=c(y,y.tmp)
  }
  y.untransformed=y
  y=log(y)
  
  
  ###
  ### Matrix building materials
  ###
  
  m.tmp1=matrix(0,6,6)
  m.tmp1[c(1,8,15)]=1
  m.tmp2=matrix(0,6,6)
  m.tmp2[c(22,29,36)]=1
  m.tmp3=matrix(0,6,6)
  m.tmp3[c(2,3,7,9,13,14)]=1
  m.tmp4=matrix(0,6,6)
  m.tmp4[c(23,24,28,30,34,35)]=1
  m.tmp5=matrix(0,6,6)
  m.tmp5[c(4,5,6,10,11,12,16,17,18)]=1
  m.tmp5[c(19,20,21,25,26,27,31,32,33)]=1
  
  ###
  ### Build Design Matrix
  ###
  
  liver=rep(rep(1:0,each=3),n)
  indices <- seq(1, nrow(data.train), by = 3)
  species=rep(data.train$species,each=2)
  flank.tmp=(data.train$flank1+data.train$flank2)/2
  flank.scaled=(flank.tmp-8821.888)/9209.355 # mean(unique(flank.tmp)); sd(unique(flank.tmp))
  axillary.tmp=(data.train$axillary1+data.train$axillary2)/2
  axillary.scaled=(axillary.tmp-13246)/15903 # 13246, 15903
  breast.tmp=(data.train$breast1+data.train$breast2)/2
  breast.scaled=(breast.tmp-9938.321)/12274.35 # 9938.321; 12274.35
  flank=rep(flank.scaled,each=2)
  axillary=rep(axillary.scaled,each=2)
  breast=rep(breast.scaled,each=2)
  
  ###
  ### Which model to fit?
  ###
  
  model.number=35
  Species=TRUE
  Liver=TRUE
  Breast=FALSE
  Flank=TRUE
  Flank2=FALSE
  
  Species.Liver=TRUE
  Species.Flank=TRUE
  Species.Axillary=FALSE
  Species.Breast=FALSE
  
  
  X=matrix(1,length(breast),1)
  if(Liver==TRUE)X=cbind(X,liver)
  if(Flank==TRUE)X=cbind(X,flank)
  if(Breast==TRUE)X=cbind(X,breast)
  if(Species==TRUE)X=cbind(X,species)
  if(Species.Liver==TRUE)X=cbind(X,species*liver)
  if(Species.Flank==TRUE)X=cbind(X,species*flank)
  if(Species.Axillary==TRUE)X=cbind(X,(1-species)*axillary)
  if(Species.Breast==TRUE)X=cbind(X,species*breast)
  
  # Set initial values, tuning parameters, and prior parameters
  initial_values <- list(
    s.l = 1,
    s.b = 1,
    rho = c(0.5, 0.5, 0.5),
    beta = rep(0, ncol(X))
  )
  tuning_parameters <- list(
    s.l.tune = 0.1,
    s.b.tune = 0.1,
    rho.tune = c(0.1, 0.1, 0.1)
  )
  prior_parameters <- list(
    s.l.prior = c(825, 100^2),
    s.b.prior = c(.01, .01),
    rho.prior = c(100, 1),
    mu.beta = matrix(0, ncol(X), 1),
    s2.beta = 1 * diag(ncol(X))
  )
  
  # Source the MCMC_Algorithm.R script from GitHub
  source(paste0('https://raw.githubusercontent.com/perrywilliamsunr/',
                'Wood-Duck-Mercury-Contamination/main/Code/MCMC_Algorithm/',
                'MCMC_Algorithm.R'))
  
  # Run the MCMC sampling
  file_name <- paste0("MCMC_CV_output_", l, ".csv")
  save_dir <- "/Users/pwill/Dropbox/GitHub/Wood-Duck-Mercury-Contamination/R_Output"
  run_mcmc(y, X, n.iter, initial_values, tuning_parameters, prior_parameters, file_name, save_dir)
  
  
}


stopCluster(myCluster)

###########################################################################
### PREDICTION
###########################################################################

###
### Matrix building materials
###

m.tmp1=matrix(0,6,6)
m.tmp1[c(1,8,15)]=1
m.tmp2=matrix(0,6,6)
m.tmp2[c(22,29,36)]=1
m.tmp3=matrix(0,6,6)
m.tmp3[c(2,3,7,9,13,14)]=1
m.tmp4=matrix(0,6,6)
m.tmp4[c(23,24,28,30,34,35)]=1
m.tmp5=matrix(0,6,6)
m.tmp5[c(4,5,6,10,11,12,16,17,18)]=1
m.tmp5[c(19,20,21,25,26,27,31,32,33)]=1

l=1

classification.within=matrix(NA,44,3)
classification.within.liver=matrix(NA,44,3)
classification.under.predict=matrix(NA,44,3)
q.breast.save=matrix(NA,44,3)
q.liver.save=matrix(NA,44,3)

for(l in 1:44){
  data=read.csv(paste0("https://github.com/perrywilliamsunr/",
                       "Wood-Duck-Mercury-Contamination/raw/",
                       "main/Data/data.csv"))
  data.ho=data[data$bird_id==l,]
  
  burn=floor(n.iter/2)
  thin=1
  ind=seq(burn,n.iter,thin)

  load(paste0("~/Dropbox/projects/Waterfowl/Ducks/",
              "Wood Ducks Mercury/Manuscript/",
              "Appendix/MCMC.Output.CV.",l,".RData"))

  beta.save=MCMC.Output$beta
  rho.save=MCMC.Output$rho
  s.l.save=MCMC.Output$s.l
  s.b.save=MCMC.Output$s.b

  ###
  ### Model prediction
  ###

  species=data.ho$species[1]
  flank.tmp=unique((data.ho$flank1+data.ho$flank2)/2)
  flank.scaled=(flank.tmp-8821.888)/9209.355 # mean(unique(flank.tmp)); sd(unique(flank.tmp))

  log.y.ppd=matrix(NA,length(ind),6)
  y.ppd=matrix(NA,length(ind),6)

  X.pred=matrix(c(rep(1,6),
                  rep(0:1,each=3),       # breast =0 or liver = 1
                  rep(flank.scaled,6), # scaled flank measurement
                  #rep(flank.scaled^2,6), # scaled flank measurement
                  rep(species,6), # species
                  rep(species,6)*rep(0:1,each=3), #Species*liver
                  rep(species,6)*rep(flank.scaled,6)),6,6) # species*flank
  
  
  
  
  for(k in 1:length(ind)){

    mu=X.pred%*%beta.save[ind[k],]

    post.Sigma.pred=s.l.save[ind[k]]^2*m.tmp1+
      s.b.save[ind[k]]^2*m.tmp2+
      rho.save[ind[k],1]*m.tmp3*s.l.save[ind[k]]^2+
      rho.save[ind[k],2]*m.tmp4*s.b.save[ind[k]]^2+
      rho.save[ind[k],3]*s.l.save[ind[k]]*s.b.save[ind[k]]*m.tmp5


    log.y.ppd[k,]=rmvnorm(1,mu,post.Sigma.pred)
    y.ppd[k,]=exp(log.y.ppd[k,])
  }


  q.b=quantile(c(log.y.ppd[,1:3]),probs=c(0.025,0.5,0.975),na.rm=TRUE)
  q.l=quantile(c(log.y.ppd[,4:6]),c(0.025,0.5,0.975),na.rm=TRUE)
  q.breast.save[l,]=q.b
  q.liver.save[l,]=q.l

  print(l)
  classification.within[l,]=q.breast.save[l,1]<=log(data.ho$y2) & 
    q.breast.save[l,3]>=log(data.ho$y2)
  classification.under.predict[l,]=q.breast.save[l,1]<=log(data.ho$y2)

  classification.within.liver[l,]=q.liver.save[l,1]<=log(data.ho$y1) & 
    q.liver.save[l,3]>=log(data.ho$y1)
  
  
}

## Which birds missed prediction?
(row.id.missed=which(apply(classification.within,1,sum)<1))
(row.id.missed.liver=which(apply(classification.within.liver,1,sum)<1))
length(row.id.missed)
problem.data=data[data$bird_id%in%row.id.missed,]


###
### Plot of Quantiles vs actual values
###

par(mfrow=c(1,1))

plot(x=1:44,y=(q.breast.save[,1]),
     ylim=c(-5,20),
     xlab="Bird ID",
     ylab="log([tHg])",
     pch=16,cex=1.5)
#points(x=1:44,y=q.breast.save[,2],pch=16)
points(x=1:44,y=(q.breast.save[,3]),pch=16,cex=1.5)


points(data$bird_id,log((data$y2)),col=2,pch=16,cex=.75)
abline(v=row.id.missed,lty=2)

# ## Add liver data
# points(data$bird_id,log((data$y1)),col=4,pch=16)

## Add flank data used to make predictions
flank.tmp=(data$flank1+data$flank2)/2
flank.scaled=(flank.tmp-8821.888)/9209.355 # mean(unique(flank.tmp)); sd(unique(flank.tmp))
points(data$bird_id,flank.scaled,col=7,pch=16)

###
### What is unique about the problem data?
###

# difference in liver and breast tissue measurements?

plot(1:length(data$flank1),flank.tmp,col=data$bird_id,pch=16)
abline(v=row.id.missed*3-2)
abline(h=0)


# difference in liver and breast tissue measurements?
plot(data$y1-data$y2,col=data$bird_id,pch=16)
abline(v=row.id.missed*3-2)
abline(h=0)

# difference in breast and flank measurements?
plot(data$y2-data$flank1,col=data$bird_id,pch=16)
abline(v=row.id.missed*3-2)
abline(h=0)

# difference in breast and mean flank measurements?
plot(data$y2-(data$flank1+data$flank2)/2,col=data$bird_id,pch=16)
abline(v=row.id.missed*3-2)
abline(h=0)

# difference in liver and flank measurements?
plot(data$y1-(data$flank1+data$flank2)/2,col=data$bird_id,pch=16)
abline(v=row.id.missed*3-2)
abline(h=0)

##########################
### Plot of prediction
##########################

head(data)

x.l=-.9#*9209.355+8821.888
x.u=4#*9209.355+8821.888
y=rep(0,44)
y[row.id.missed]=1
y.b.mat=matrix(data$y2,3,)
y.b.m=colMeans(y.b.mat)
flank=unique(flank.scaled)
flank.pred=seq(x.l,x.u,.01)


##################################################
### Start of Figure
##################################################
alpha=0.07

par(mfrow=c(2,2),mar=c(4,4,1,1))

###
### Wood duck breast tissue
###

plot(NA,NA,
     ylim=c(0,15),xlim=c(x.l,x.u),
     ylab="log([tHg]) in breast",
     xlab="",
     xaxt='n',
     main="Wood Duck Breast")
# wood duck
for(i in 1:nrow(beta.save)){
  mu.pred=beta.save[i,1] + 
    beta.save[i,2]*0 + # 0 = breast, 1=liver
    beta.save[i,3] * flank.pred + 
    beta.save[i,4]*0 + # 0 = wood duck, 1 = mallard
    beta.save[i,5]*0*0 + # species * liver
    beta.save[i,6]*0*flank.pred  # species*flank
  lines(flank.pred,mu.pred,col=rgb(205/255, 175/255, 100/255, alpha = .1))  
}
points(flank[1:39],log(y.b.m[1:39]),
       col=rgb(0.545, 0.271, 0.075),
       pch=16)
points(flank[row.id.missed[row.id.missed<40]],
       log(y.b.m[row.id.missed[row.id.missed<40]]),
       pch=16,
       col=rgb(0, 0.749, 1),
)
legend(0,4,legend=c("Within 95% CRI", "Outside 95% CRI"), 
       col=c(rgb(0.545, 0.271, 0.075), rgb(0, 0.749, 1)),
       pch=16,title="Leave-one-out cross-validation")



###
### Mallard breast tissue
###

plot(NA,NA,
     ylim=c(0,15),xlim=c(x.l,x.u),
     ylab="",
     xlab="",
     xaxt='n',
     main="Mallard Breast")


for(i in ind){
  mu.pred=beta.save[i,1] + 
    beta.save[i,2]*0 + # 0 = breast, 1=liver
    beta.save[i,3] * flank.pred + 
    beta.save[i,4]*1 + # 0 = wood duck, 1 = mallard
    beta.save[i,5]*0*0 + # species * liver
    beta.save[i,6]*1*flank.pred  # species*flank
  lines(flank.pred,mu.pred,col=rgb(0, 0.502, 0, alpha = alpha))  
}

points(flank[40:44],log(y.b.m[40:44]),
       col=rgb(0.545, 0.271, 0.075),
       pch=16)
points(flank[row.id.missed[row.id.missed>39]],
             log(y.b.m[row.id.missed[row.id.missed>39]]),
             pch=16,
             col=rgb(0, 0.749, 1)
)

###
### Wood duck liver tissue
###

plot(NA,NA,
     ylim=c(0,15),xlim=c(x.l,x.u),
     ylab="log([tHg]) in liver",
     main="Wood Duck Liver",
     xlab="Flank feather [tHg] (z-scaled)")
# wood duck
for(i in 1:nrow(beta.save)){
  mu.pred=beta.save[i,1] + 
    beta.save[i,2]*1 + # 0 = breast, 1=liver
    beta.save[i,3] * flank.pred + 
    beta.save[i,4]*0 + # 0 = wood duck, 1 = mallard
    beta.save[i,5]*0*1 + # species * liver
    beta.save[i,6]*0*flank.pred  # species*flank
  lines(flank.pred,mu.pred,col=rgb(205/255, 175/255, 100/255, alpha = .1))  
}

y.l.mat=matrix(data$y1,3,)
y.l.m=colMeans(y.l.mat)


points(flank[1:39],log(y.l.m[1:39]),
       col=rgb(0.545, 0.271, 0.075),
       pch=16)
points(flank[row.id.missed.liver[row.id.missed<40]],
       log(y.l.m[row.id.missed.liver[row.id.missed.liver<40]]),
       pch=16,
       col=rgb(0, 0.749, 1),
)


###
### Mallard liver tissue
###

plot(NA,NA,
     ylim=c(0,15),xlim=c(x.l,x.u),
     ylab="",
     xlab="Flank feather [tHg] (z-scaled)",
     main="Mallard Liver")
# mallard
for(i in ind){
  mu.pred=beta.save[i,1] + 
    beta.save[i,2]*1 + # 0 = breast, 1=liver
    beta.save[i,3] * flank.pred + 
    beta.save[i,4]*1 + # 0 = wood duck, 1 = mallard
    beta.save[i,5]*1*1 + # species * liver
    beta.save[i,6]*1*flank.pred  # species*flank
  lines(flank.pred,mu.pred,col=rgb(0, 0.502, 0, alpha = alpha))  
}

points(flank[40:44],log(y.l.m[40:44]),
       col=rgb(0.545, 0.271, 0.075),
       pch=16)
points(flank[row.id.missed.liver[row.id.missed.liver>39]],
       log(y.l.m[row.id.missed.liver[row.id.missed.liver>39]]),
       pch=16,
       col=rgb(0, 0.749, 1)
)


##################################################
### End of Figure
##################################################


##################################################
### Beta Inference
##################################################

beta0=beta.save[ind,1]
beta.liver=beta.save[ind,2]
beta.flank=beta.save[ind,3]
beta.species=beta.save[ind,4]
beta.species.liver=beta.save[ind,5]
beta.species.flank=beta.save[ind,6]

## increase in wood duck breast tissue for each increase in flank tissue sd
quantile(beta.flank,c(0.025,0.5,0.975))

## increase in mallard breast tissue for each increase in flank tissue sd
quantile(beta.flank+beta.species.flank,c(0.025,0.5,0.975))

## Difference between liver [thg] and breast [tHg] in wood ducks
quantile(beta.liver,c(0.025,0.5,0.975))

## Difference between liver [thg] and breast [tHg] in mallards
quantile(beta.liver+beta.species.liver,c(0.025,0.5,0.975))

## Which species has higher [tHg]?
quantile(beta.species,c(0.025,0.5,0.975))

## slope difference for flank for mallards
quantile(beta.species.flank,c(0.025,0.5,0.975))


















###
### Quantiles of the full analysis
###

setwd("~/Dropbox/projects/Waterfowl/Ducks/Wood Ducks Mercury/Manuscript/Appendix")
load("MCMC.Output.Real.35.RData")

beta.save=MCMC.Output$beta
rho.save=MCMC.Output$rho
s.l.save=MCMC.Output$s.l
s.b.save=MCMC.Output$s.b

###
### Model prediction
###

species=0
flank=(unique(data$flank1)+unique(data$flank2))/2
flank.mean=mean(flank)
flank.sd=sd(flank)
flank.tmp=mean(c(data$flank1,data$flank2))
flank=(flank.tmp-flank.mean)/flank.sd

log.y.ppd=matrix(NA,length(ind),6)
y.ppd=matrix(NA,length(ind),6)

X.pred=matrix(c(rep(1,6),
                rep(species,6), # species
                rep(0:1,each=3),       # breast =0 or liver = 1
                rep(flank,6), # scaled flank measurement
                rep(species,6)*rep(0:1,each=3), #Species*liver
                rep(species,6)*rep(flank,6)),6,6) # species*flank

for(k in 1:length(ind)){


  mu=X.pred%*%beta.save[ind[k],]

  post.Sigma.pred=s.l.save[ind[k]]^2*m.tmp1+
    s.b.save[ind[k]]^2*m.tmp2+
    rho.save[ind[k],1]*m.tmp3*s.l.save[ind[k]]^2+
    rho.save[ind[k],2]*m.tmp4*s.b.save[ind[k]]^2+
    rho.save[ind[k],3]*s.l.save[ind[k]]*s.b.save[ind[k]]*m.tmp5


  log.y.ppd[k,]=rmvnorm(1,mu,post.Sigma.pred)
  y.ppd[k,]=exp(log.y.ppd[k,])
}


classification.within[l,]=q.breast.save[l,1]<=data.ho$y2 &q.breast.save[l,3]>=data.ho$y2
classification.under.predict[l,]=q[1,]<=c(data.ho$y2,data.ho$y1)&q[2,]



## Average successful prediction rate within 95% CI
mean(classification.within)

## Average successful prediction rate lower bound 95% CI
mean(classification.under.predict)

## Average successful prediction rate of breast tissue
mean(classification.within[,4:6])
mean(classification.under.predict[,4:6])

## Average successful prediction rate of liver tissue
mean(classification.within[,1:3])
mean(classification.under.predict[,1:3])



##
## Plot problem data vs good data
##

par(mfrow=c(1,1))
m.flank=scale((data$flank1+data$flank2)/2)
y=data$y2
plot(m.flank,log(y),pch=16)
m2.flank=(((problem.data$flank1+problem.data$flank2)/2)-8821.888)/9138.784
y2=log(problem.data$y2)
points(m2.flank,y2,col=2,pch=16)

m.flank=scale((data$flank1+data$flank2)/2)
text(m2.flank,y,data$bird_id)


## Properties of missed predictions
par(mfrow=c(1,1))
plot(scale((data$flank1+data$flank2)/2),log(data$y1))
points(scale((problem.data$flank1+problem.data$flank2)/2),log(problem.data$y1),col=2)

par(mfrow=c(2,2))
image(y.m)
image(classification.within)
image(classification.under.predict)
y.miss=y.m[!classification.within]
image(y.miss)


###
###
###

l=row.id.missed[2]
load(paste0("~/Dropbox/projects/Waterfowl/Ducks/",
            "Wood Ducks Mercury/Manuscript/",
            "Appendix/MCMC.Output.CV.",l,".RData"))



data.ho=data[data$bird_id==l,]
m.flank.data.ho=mean(unlist(data.ho[,4:5]))
flank=(m.flank.data.ho-MCMC.Output$flank[1])/MCMC.Output$flank[2]
species.ho=0

x=c(1,0,0,flank,0,0)
dim(MCMC.Output$beta)
l.mu.pred=MCMC.Output$beta%*%x
plot(density(l.mu.pred),xlim=c(5,8))
abline(v=log(1394.581),col=2)


plot(data$flank2,data$y2)
for(i in 1:12){
  l=row.id.missed[i]
  load(paste0("~/Dropbox/projects/Waterfowl/Ducks/",
              "Wood Ducks Mercury/Manuscript/",
              "Appendix/MCMC.Output.CV.",l,".RData"))



  data.ho=data[data$bird_id==l,]

  points(data.ho$flank2,data.ho$y2,col=2,pch=16)
}

mu.pred=exp(l.y.pred)
hist(mu.pred)
hist(l.mu.pred)
data.ho

plot(scale(data$flank1),log(data$y1),xlim=c(-3,3),ylim=c(0,10))
points(flank,log(1422.837),col=2,pch=16)




apply(MCMC.Output$y.ppd,2,quantile,c(0.025,0.975))


f


















###
### FIGURE 1
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


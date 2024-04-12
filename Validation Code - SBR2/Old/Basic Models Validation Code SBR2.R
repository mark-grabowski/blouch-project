#Basic Model Validation Code
#SBR2 - 022924
rm(list=ls())
########################################################################################################
#Complete code to make simulated data, setup, run analyses
#Also need to load blouch.prep.R code to aid in simulations (calculate Vs, etc.), and set.converge.regimes.R
#Prior and Posterior Plots Code.Rmd - to make plots
########################################################################################################
#Basic Setup
library(geiger)
library(MASS)
library(phytools)
library(ggplot2)
library(rstan)
########################################################################################################
#For execution on a local, multicore CPU with excess RAM we recommend calling
options(mc.cores = parallel::detectCores())
#options(mc.cores = 8)
rstan_options(auto_write = TRUE)

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
arch <- ifelse(R.version$arch == "aarch64", "arm64", "x86_64")
cat(paste("\nCXX14FLAGS += -O3 -mtune=native -arch", arch, "-ftemplate-depth-256"),
    file = M, sep = "\n", append = FALSE)
########################################################################################################
#source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/R Setup Code/blouch.prep.SBR2.R")
tree.10K<-read.tree("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/blouch project - for ms - SBR2/Sharable Data/10KPrimateTree.tre")
source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/R Packages - SBR2 - cloned/testing/R/set.converge.regimes.R")
source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/R Packages - SBR2 - cloned/testing/R/blouch.prep.R")
########################################################################################################
#Create phylogeny - for all models
########################################################################################################
N<-50 #Number of species
set.seed(10) #Set seed to get same random species each time

phy <- ape::keep.tip(tree.10K,sample(tree.10K$tip.label)[1:N]) 
phy<-ape::multi2di(phy)

l.tree<-max(ape::branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree 

tip.label<-phy$tip.label

########################################################################################################


########################################################################################################
#Direct Effect Model - with Measurement Error
########################################################################################################
#Setup parameters
Z<-1 #Number of traits
hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.01 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));

vX0<-0
vY0 <- 0

Sxx<-10 #Look at effects

V<-calc_direct_V(phy,sigma2_y,a)
X<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

optima<-2 #Intecept
beta<-0.25 #Slope

mu<-optima+X*beta #Simulate mu for Y
#Simulate direct effect Y trait
Y<-MASS::mvrnorm(n=1,mu,V)

df<-data.frame(Y=Y,X=X)
names(df)<-c("Y","X")

ggplot2::ggplot(data=df,ggplot2::aes(x=X,y=Y))+
  ggplot2::geom_point()
summary(lm(Y~X,df))

X_error<-rep(0.01,N)
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-X+rnorm(N,0,0.01)

############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-5
optima.prior<-c(2,0.25) #Informed by linear model
beta.prior<-c(0.25,0.25) #Informed by linear model

#View priors in Prior and Posterior Plots Code.Rmd

############################################################################################################
#Make trdata file
trait.data<-data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error))
trdata<-treeplyr::make.treedata(phy,trait.data)

############################################################################################################
#Blouch prep code - direct effect model - blouch.direct.prep()
dat<-blouch.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",hl.prior,vy.prior,optima.prior,beta.prior)
############################################################################################################
#Run Blouch Direct Effect model with measurement error
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_direct.stan") #Macbook Pro
  
stan_model <- stan_model("blouchOU_direct.stan")
  
fit.npi.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter = 2000)
  
print(fit.npi.direct,pars = c("hl","vy","optima","beta"))
post<-extract(fit.npi.direct)

########################################################################################################
#Direct effect model with multiple X traits
########################################################################################################
#Simulate multiple X traits - non-correlated
Z<-2
hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.1 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));

vX0<-0
vY0 <- 0
Sxx<-10 #Look at effe
ts<-ts_fxn(phy)
ta<-ts[[1]]
tij<-ts[[2]]

V<-calc_direct_V(phy,sigma2_y,a)
#vcv<-matrix(c(1,0.75,0.75,1),2,2) #Correlation between traits
vcv<-matrix(c(1,0,0,1),2,2) #No correlation between traits

Xs<-sim.corrs(phy,vcv) #Simulated correlated BM Xs
optima<-2
beta<-c(0.35,0.25) #Slope
mu<-optima+Xs%*%beta #Simulate mu for Y

#Simulate direct effect Y trait
Y<-mvrnorm(n=1,mu,V)
df<-data.frame(Y=Y,X=Xs)

ggplot(data=df,aes(x=X.1,y=X.2))+
  geom_point()
summary(lm(X.2~X.1,df))

ggplot(data=df,aes(x=X.1,y=Y))+
  geom_point()
summary(lm(Y~X.1,df))

ggplot(data=df,aes(x=X.2,y=Y))+
  geom_point()
summary(lm(Y~X.2,df))

########################################################################################################
#Simulate errors with multiple traits - original Hansen setup
X_error<-matrix(0.01,nrow=N,ncol=Z)
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-apply(Xs,2,function(X){X+rnorm(N,0,0.01)})


############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-5
optima.prior<-c(2,0.2) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model

#View priors in Prior and Posterior Plots Code.Rmd

############################################################################################################
#Make trdata file
trait.data<-data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error))
trdata<-treeplyr::make.treedata(phy,trait.data)
names(trdata$dat)[3]<-"X_with_error_1"
names(trdata$dat)[4]<-"X_with_error_2"
names(trdata$dat)[5]<-"X_error_1"
names(trdata$dat)[6]<-"X_error_2"
############################################################################################################
#Blouch prep code - direct effect model - blouch.direct.prep()
dat<-blouch.direct.prep(trdata,"Y_with_error","Y_error",c("X_with_error_1","X_with_error_2"),c("X_error_1","X_error_2"),hl.prior,vy.prior,optima.prior,beta.prior)
############################################################################################################
#Run Blouch Direct Effect model with measurement error
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_direct.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_direct.stan")

fit.npi.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter = 2000)

print(fit.npi.direct,pars = c("hl","vy","optima","beta"))
post<-extract(fit.npi.direct)


########################################################################################################
#Adaptive Model - with Measurement Error
########################################################################################################
#Adaptive Model
#Setup parameters
Z_adaptive<-1 #Number of traits
hl<-0.1
a<-log(2)/hl
vy<-0.1 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));
vX0<-0
vY0 <- 0
Sxx<-10 #Look at effects

X<-fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
sigma2_x<-matrix(1,1,1)

optima<-2 #Intecept
beta<-0.25 #Slope

dmX<-calc_adaptive_dmX(phy,a,X) #Calculate the design matrix
mu<-optima+dmX%*%beta #Simulate mu for Y

V<-calc_adaptive_V(phy,a, sigma2_y, beta,  sigma2_x, Z_adaptive)

#Simulate direct effect Y trait
Y<-mvrnorm(n=1,mu,V)

df<-data.frame(Y=Y,X=X)
names(df)<-c("Y","X")

ggplot(data=df,aes(x=X,y=Y))+
  geom_point()

summary(lm(Y~X,df))

########################################################################################################
#Simulate errors
X_error<-rep(0.01,N)
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-X+rnorm(N,0,0.01)

########################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-5
optima.prior<-c(2,0.2) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
#View priors in Prior and Posterior Plots Code.Rmd

############################################################################################################
#Make trdata file
trait.data<-data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error))
trdata<-treeplyr::make.treedata(phy,trait.data)

############################################################################################################
#Blouch prep code
dat<-blouch.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",hl.prior,vy.prior,optima.prior,beta.prior)
############################################################################################################
#Run Blouch Adaptive model with measurement error
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_adapt.stan") #Macbook Pro

stan_model <- rstan::stan_model("blouchOU_adapt.stan")

fit.npi.adaptive<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000)

print(fit.npi.adaptive,pars = c("hl","vy","optima","beta","beta_e"))
post<-extract(fit.npi.adaptive)

############################################################################################################
#Simulate multiple adaptive traits - not correlated
Z<-2
hl<-0.1
a<-log(2)/hl
vy<-0.1 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));
vX0<-0
vY0 <- 0
optima<-2 #Intecept

vcv<-matrix(c(1,0,0,1),2,2) #No correlation between traits
Xs<-sim.corrs(phy,vcv) #Simulated correlated BM Xs
sigma2_x<-ratematrix(phy,Xs) #Calculate evolutionary v/cv matrix

beta<-c(0.35,0.25) #Slope
dmX<-calc_adaptive_dmX(phy,a,Xs) #Calculate the design matrix
mu<-optima+dmX%*%beta #Simulate mu for Y

V<-calc_adaptive_V(phy,a, sigma2_y, beta,  sigma2_x, Z_adaptive)
Y<-mvrnorm(n=1,mu,V)

df<-data.frame(Y=Y,X=Xs)

ggplot(data=df,aes(x=X.1,y=X.2))+
  geom_point()
summary(lm(X.2~X.1,df))

ggplot(data=df,aes(x=X.1,y=Y))+
  geom_point()

summary(lm(Y~X.1,df))

ggplot(data=df,aes(x=X.2,y=Y))+
  geom_point()

summary(lm(Y~X.2,df))

########################################################################################################
#Simulate errors with multiple traits - original Hansen setup
X_error<-matrix(0.01,nrow=N,ncol=Z)
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-apply(Xs,2,function(X){X+rnorm(N,0,0.01)})

############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-5
optima.prior<-c(2,0.2) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model

#View priors in Prior and Posterior Plots Code.Rmd

############################################################################################################
#Make trdata file
trait.data<-data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error))
trdata<-treeplyr::make.treedata(phy,trait.data)
names(trdata$dat)[3]<-"X_with_error_1"
names(trdata$dat)[4]<-"X_with_error_2"
names(trdata$dat)[5]<-"X_error_1"
names(trdata$dat)[6]<-"X_error_2"
############################################################################################################
#Blouch prep code
dat<-blouch.adapt.prep(trdata,"Y_with_error","Y_error",c("X_with_error_1","X_with_error_2"),c("X_error_1","X_error_2"),hl.prior,vy.prior,optima.prior,beta.prior)
############################################################################################################
#Run Blouch Adaptive model with measurement error
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_adapt.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_adapt.stan")

fit.npi.adaptive<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter = 2000)

print(fit.npi.adaptive,pars = c("hl","vy","optima","beta","beta_e"))
post<-extract(fit.npi.adaptive)

############################################################################################################
#Combination Direct effect and Adaptive Model
############################################################################################################
#Direct effect + Adaptive Model
############################################################################################################
hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.1 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));

vX0<-0
vY0 <- 0

ts<-ts_fxn(phy)
ta<-ts[[1]]
tij<-ts[[2]]
T_term<-ts[[3]]
tja<-ts[[4]]
Z_direct<-1
Z_adaptive<-1
Z<-Z_direct+Z_adaptive
sigma2_x<-matrix(1,1,1)

Xd<-rnorm(N,0,1)
names(Xd)<-phy$tip.label
phenogram(phy,Xd,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
Xa<-fastBM(phy,a=vX0,sig2=sigma2_x[1,1],internal=FALSE) #Simulate X BM variable on tree, with scaling 10
phenogram(phy,Xa,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
#sigma2_x<-ratematrix(phy,Xa) #Calculate evolutionary v/cv matrix
Xs<-cbind(Xd,Xa)

alpha<-2 #Intecept
beta<-c(0.35,0.25) #Slopes
dmX<-calc_mixed_dmX(phy,a,Xs,Z_direct,Z_adaptive)
mu<-alpha+dmX%*%beta #Simulate mu for Y

V<-calc_adaptive_V(phy,a, sigma2_y, beta[(Z_direct+1):(Z_adaptive+Z_direct)],sigma2_x,Z_adaptive)
Y<-MASS::mvrnorm(n=1,mu,V)

df<-data.frame(Y=Y,X=Xs)

ggplot(data=df,aes(x=X.Xd,y=X.Xa))+
  geom_point()
summary(lm(X.Xa~X.Xd,df))

ggplot(data=df,aes(x=X.Xd,y=Y))+
  geom_point()

summary(lm(Y~X.Xd,df))

ggplot(data=df,aes(x=X.Xa,y=Y))+
  geom_point()

summary(lm(Y~X.Xa,df))

########################################################################################################
#Simulate errors - for use with blouchOU_reg_direct_adaptive_ME
Z_X_error<-2 #Number of X traits with error
X_error<-matrix(0.01,nrow=N,ncol=Z_X_error)
X_error<-data.frame(X_error)
names(X_error)<-c("Xd_error","Xa_error")
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-apply(Xs,2,function(X){X+rnorm(N,0,0.01)})

############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-5
optima.prior<-c(2,0.2) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model

#View priors in Prior and Posterior Plots Code.Rmd
############################################################################################################
#Make trdata file
trait.data<-data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error))
trdata<-treeplyr::make.treedata(phy,trait.data)
############################################################################################################
#Blouch prep code
dat<-blouch.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Run combination of direct effect and adaptive predictors with ME
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_direct_adapt.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_direct_adapt.stan")

fit.npi.mixed<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.npi.mixed,pars = c("hl","vy","optima","beta","beta_e"))
post<-extract(fit.npi.mixed)
########################################################################################################
#Multi-optima Models
############################################################################################################

#Set regimes - manually - 2 regimes
#Locate nodes
plot(phy,no.margin=TRUE,edge.width=2,cex=0.7)
ape::nodelabels(frame="none",adj=c(1.1,-0.4))
ape::tiplabels()

shifts<-c(84) #Location of nodes with regime shifts
trdata<-data.frame(phy$tip.label)
trdata<-treeplyr::make.treedata(phy,trdata)
trdata<-set.converge.regimes(trdata,shifts)

#Check if manual setting code worked
shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
print(edge.regimes)

reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(2)
print(reg.colors)
plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

############################################################################################################
#Format tree with lineage data
n<-length(trdata$phy$tip.label)
regimes_internal <-trdata$phy$node.label
regimes_tip <- trdata$dat$regimes
regimes <- concat.factor(regimes_tip, regimes_internal)
anc_maps<-"regimes"
lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

############################################################################################################
#Now we will simulate Y data based on our generative model
############################################################################################################
#Set true values for parameters
hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.01 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));
optima<-c(0.5,0.25) #Optima for two regimes

dmX<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
mu<-dmX%*%optima #Simulate mu for Y
V<-calc_direct_V(phy, sigma2_y, a)
Y<-MASS::mvrnorm(n=1,mu,V)

########################################################################################################
#Simulate errors
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)

phytools::phenogram(phy,Y_with_error,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model

############################################################################################################
#Make trdata file
trait.data<-data.frame(cbind(Y_with_error,Y_error))
trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error)))

############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior)

########################################################################################################
#Run Multi-optima Model
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg.stan")

fit.reg<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg,pars = c("hl","vy","optima"))
post<-extract(fit.reg)

############################################################################################################
#Multilevel Multi-Optima Model
############################################################################################################
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
sigma.prior<-1
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.mlm.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior,sigma.prior)

########################################################################################################
#Run Multi-optima Model
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_mlm_vi.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_mlm_vi.stan")

fit.mlm.reg<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.mlm.reg,pars = c("hl","vy","optima","sigma"))
post<-extract(fit.mlm.reg)

############################################################################################################
#Multilevel Multi-Optima Model - Non-Centered
############################################################################################################
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-5
optima.prior<-c(0,1) #Informed by linear model
sigma.prior<-1
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.mlm.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior,sigma.prior)

########################################################################################################
#Run Multi-optima Model
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_mlm_vi_nc.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_mlm_vi_nc.stan")

fit.mlm.nc.reg<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.mlm.nc.reg,pars = c("hl","vy","optima","sigma"))
post<-extract(fit.mlm.nc.reg)


############################################################################################################
############################################################################################################
#Multi-Optima + Direct Effect Model
############################################################################################################
#Set regimes - manually - 2 regimes
#Locate nodes
plot(phy,no.margin=TRUE,edge.width=2,cex=0.7)
ape::nodelabels(frame="none",adj=c(1.1,-0.4))
ape::tiplabels()

#Paint Regimes on Tree
#source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Simulation Code/Functions/set.converge.regimes.R") #Macbook Pro
shifts<-c(84) #Location of nodes with regime shifts
trdata<-data.frame(phy$tip.label)
trdata<-treeplyr::make.treedata(phy,trdata)
trdata<-set.converge.regimes(trdata,shifts)


#Check if manual setting code worked
shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
print(edge.regimes)

reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(2)
print(reg.colors)
plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

#Phylogeny info
n<-length(trdata$phy$tip.label)
mrca1 <- ape::mrca(trdata$phy)
times <- ape::node.depth.edgelength(trdata$phy)

regimes_internal <-trdata$phy$node.label
regimes_tip <- trdata$dat$regimes
regimes <- concat.factor(regimes_tip, regimes_internal)
anc_maps<-"regimes"
lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

###########################################################################
#Simulate Data
###########################################################################
hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.01 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));

vX0<-0
vY0 <- 0
Sxx<-10

Z_direct<-1

V<-calc_direct_V(phy,sigma2_y,a)
X<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
#X<-rnorm(N,0,1)
names(X)<-phy$tip.label

phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

dmX<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
dmX<-cbind(dmX,X)

beta<-c(2,1,0.25) #Two Optima/One Slope
mu<-dmX%*%beta #Simulate mu for Y

V<-calc_direct_V(phy,sigma2_y,a)
Y<-MASS::mvrnorm(n=1,mu,V)

#Plot data
df<-data.frame(Y=Y,X=X)

ggplot2::ggplot(data=df,ggplot2::aes(x=X,y=Y))+
  ggplot2::geom_point()
summary(lm(Y~X,df))

#################################################################################################################Simulate errors
Z_X_error<-1 #Number of X traits with error
X_error<-rep(0.01,N)
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-X+rnorm(N,0,0.01)

############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
beta.prior<-c(0.25,1.25) #Informed by linear model

############################################################################################################
#Make trdata file
trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))

############################################################################################################
dat<-blouch.reg.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)

############################################################################################################
#Run Multi-Optima Direct Effect Model
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_direct.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_direct.stan")

fit.reg.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.direct,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.direct)

############################################################################################################
#Multilevel Multi-Optima Direct Effect Model
############################################################################################################
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
sigma.prior<-1

############################################################################################################
#Make trdata file
trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))

dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)

############################################################################################################
#Run Multilevel Multi-Optima Direct Effect Model
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_direct_mlm_vi.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_direct_mlm_vi.stan")

fit.reg.direct.mlm<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.direct.mlm,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.direct.mlm)

############################################################################################################
#Run Multilevel Multi-Optima Direct Effect Model - Non-centered
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_direct_mlm_vi_nc.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_direct_mlm_vi_nc.stan")

fit.reg.direct.mlm.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.direct.mlm.nc,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.direct.mlm.nc)


########################################################################################################################################################################################################################
#Multi-optima Adaptive Model
############################################################################################################
#Set regimes - manually - 2 regimes
#Locate nodes
plot(phy,no.margin=TRUE,edge.width=2,cex=0.7)
ape::nodelabels(frame="none",adj=c(1.1,-0.4))
ape::tiplabels()

#Paint Regimes on Tree
shifts<-c(84) #Location of nodes with regime shifts
trdata<-data.frame(phy$tip.label)
trdata<-treeplyr::make.treedata(phy,trdata)
trdata<-set.converge.regimes(trdata,shifts)

#Check if manual setting code worked
shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
print(edge.regimes)

reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(2)
print(reg.colors)
plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

#Phylogeny info
n<-length(trdata$phy$tip.label)

regimes_internal <-trdata$phy$node.label
regimes_tip <- trdata$dat$regimes
regimes <- concat.factor(regimes_tip, regimes_internal)
anc_maps<-"regimes"
lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

############################################################################################################
########################################################################################################
ts<-ts_fxn(phy)
ta<-ts[[1]]
tij<-ts[[2]]
T_term<-ts[[3]]
tja<-ts[[4]]

#########################
hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.01 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));

vX0<-0
vY0 <- 0

Z_adaptive<-1
sigma2_x<-matrix(10,1,1)

X<-phytools::fastBM(phy,a=vX0,sig2=sigma2_x[1,1],internal=FALSE) #Simulate X BM variable on tree, with scaling 10
phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
#sigma2_x<-ratematrix(phy,Xa) #Calculate evolutionary v/cv matrix

dmX<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
dmX<-cbind(dmX,calc_adaptive_dmX(phy,a,X))

beta<-c(2,1,0.25) #Two Optima/Two Slopes
mu<-dmX%*%beta #Simulate mu for Y

V<-calc_adaptive_V(phy,a, sigma2_y, beta[3],sigma2_x,Z_adaptive)

Y<-mvrnorm(n=1,mu,V)

################################################################################################################
#Plot data
df<-data.frame(Y=Y,X=X)

ggplot2::ggplot(data=df,ggplot2::aes(x=X,y=Y))+
  ggplot2::geom_point()
summary(lm(Y~X,df))
################################################################################################################
#Simulate errors
Z_X_error<-1 #Number of X traits with error
X_error<-rep(0.01,N)
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-X+rnorm(N,0,0.01)

############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-5
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Make trdata file
trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))
############################################################################################################
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)

############################################################################################################
#Run Multi-Optima Adaptive Model
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_adapt.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_adapt.stan")

fit.reg.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.adapt,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.adapt)

############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-5
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-1

############################################################################################################
#Make trdata file
trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))
############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)

############################################################################################################
#Run Multi-Optima Adaptive Model
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_adapt_mlm_vi.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_adapt_mlm_vi.stan")

fit.reg.adapt.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.adapt.mlm.vi,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.adapt.mlm.vi)

############################################################################################################
#Run Multi-Optima Adaptive Model - Non-centered
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_adapt_mlm_vi_nc.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_adapt_mlm_vi_nc.stan")

fit.reg.adapt.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.adapt.mlm.vi.nc,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.adapt.mlm.vi.nc)


############################################################################################################
#Multi-optima Direct Effect and Adaptive Model
############################################################################################################
N<-50 #Number of species
set.seed(10) #Set seed to get same random species each time

phy <- ape::keep.tip(tree.10K,sample(tree.10K$tip.label)[1:N]) 
phy<-ape::multi2di(phy)

l.tree<-max(ape::branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree 

tip.label<-phy$tip.label

#Set regimes - manually - 2 regimes
#Locate nodes
plot(phy,no.margin=TRUE,edge.width=2,cex=0.7)
ape::nodelabels(frame="none",adj=c(1.1,-0.4))
ape::tiplabels()

#Paint Regimes on Tree
shifts<-c(84) #Location of nodes with regime shifts
trdata<-data.frame(phy$tip.label)
trdata<-treeplyr::make.treedata(phy,trdata)
trdata<-set.converge.regimes(trdata,shifts)

#Check if manual setting code worked
shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
print(edge.regimes)

reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(2)
print(reg.colors)
plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

#Phylogeny info
n<-length(trdata$phy$tip.label)

regimes_internal <-trdata$phy$node.label
regimes_tip <- trdata$dat$regimes
regimes <- concat.factor(regimes_tip, regimes_internal)
anc_maps<-"regimes"
lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch
############################################################################################################
hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.01 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));

vX0<-0
vY0 <- 0

Z_direct<-1
Z_adaptive<-1
Z<-Z_direct+Z_adaptive
sigma2_x<-matrix(1,1,1) #Variance of BM Process

Xd<-rnorm(N,0,1)
names(Xd)<-phy$tip.label
phytools::phenogram(phy,Xd,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
Xa<-phytools::fastBM(phy,a=vX0,sig2=sigma2_x[1,1],internal=FALSE) #Simulate X BM variable on tree, with scaling 10
phytools::phenogram(phy,Xa,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
#sigma2_x<-ratematrix(phy,Xa) #Calculate evolutionary v/cv matrix
Xs<-cbind(Xd,Xa)

dmX<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
dmX<-cbind(dmX,calc_mixed_dmX(phy,a,Xs,Z_direct,Z_adaptive))

beta<-c(2,1,0.35,0.25) #Two Optima/Two Slopes
mu<-dmX%*%beta #Simulate mu for Y

V<-calc_adaptive_V(phy,a, sigma2_y,   beta[length(beta)],  sigma2_x,Z_adaptive)
Y<-MASS::mvrnorm(n=1,mu,V)

################################################################################################################
#Plot data
df<-data.frame(Y=Y,Xd=Xs[,1],Xa=Xs[,2])

ggplot2::ggplot(data=df,ggplot2::aes(x=Xd,y=Y))+
  ggplot2::geom_point()

ggplot2::ggplot(data=df,ggplot2::aes(x=Xa,y=Y))+
  ggplot2::geom_point()
summary(lm(Y~Xs,df))

########################################################################################################
#Simulate errors - for use with blouchOU_reg_direct_adaptive_ME
Z_X_error<-2 #Number of X traits with error
X_error<-matrix(0.01,nrow=N,ncol=Z_X_error)
X_error<-data.frame(X_error)
names(X_error)<-c("Xd_error","Xa_error")
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-apply(Xs,2,function(X){X+rnorm(N,0,0.01)})

############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-5
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Make trdata file
trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))

############################################################################################################
dat<-blouch.reg.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)

############################################################################################################
#Run Multi-Optima Direct Effect and Adaptive Model
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_direct_adapt.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_direct_adapt.stan")

fit.reg.direct.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.direct.adapt,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.direct.adapt)

############################################################################################################


#Next we run the Multilevel version of the same model
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-5
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-5
############################################################################################################
#Make trdata file
trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))


############################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)

############################################################################################################
#Run Multi-Optima Direct Effect and Adaptive Model
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_direct_adapt_mlm_vi.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_direct_adapt_mlm_vi.stan")

fit.reg.direct.adapt.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.direct.adapt.mlm.vi,pars = c("hl","vy","optima","beta","sigma"))
post<-extract(fit.reg.direct.adapt.mlm.vi)

############################################################################################################
#Run Multi-Optima Direct Effect and Adaptive Model - Non-centered
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_direct_adapt_mlm_vi_nc.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_direct_adapt_mlm_vi_nc.stan")

fit.reg.direct.adapt.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.direct.adapt.mlm.vi,pars = c("hl","vy","optima","beta","sigma"))
post<-extract(fit.reg.direct.adapt.mlm.vi)



#Varying Effects Models Validation Code
#SBR2 - 022924
rm(list=ls())
########################################################################################################
#Complete code to make simulated data, setup, run analyses
#Also need to load blouch.prep.R code to aid in simulations (calculate Vs, etc.), and set.converge.regimes.R
#Use Prior and Posterior Plots Code.Rmd - to make plots
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
###########################################################################################################
#Multi-Optima Direct Effect - Varying Effects Model
########################################################################################################
set.seed(10)

N<-50 #Number of species
#set.seed(1) #Set seed to get same random species each time

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

reg_tips<-trdata$dat$regimes
reg_tips<-as.numeric(as.factor(reg_tips))

reg.colors<-ggsci::pal_npg(palette=c("nrc"),alpha=1)(2)
plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

#Phylogeny info
n<-length(trdata$phy$tip.label)

regimes_internal <-trdata$phy$node.label
regimes_tip <- trdata$dat$regimes
regimes <- concat.factor(regimes_tip, regimes_internal)

anc_maps<-"regimes"
lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

####################################################################################################
hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.01 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));

vX0<-0
vY0 <- 0
Sxx<-10 #Look at effects

X<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
Z_direct<-1
names(X)<-phy$tip.label
phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

dmX<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
mu<-matrix(NA,N,1)
for(i in 1:N){
  mu[i]<-dmX[i,]%*%optima+beta[reg_tips[i]]%*%X[i];
}

V<-calc_direct_V(phy,sigma2_y,a)
Y<-MASS::mvrnorm(n=1,mu,V)

###############################################################################################################
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
beta.prior<-c(0,0.25)
############################################################################################################
#Make trdata file
trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))

############################################################################################################
#Test Blouch prep code - Regimes + Direct Efffect model
#source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/R Setup Code/blouch.prep.R")
dat<-blouch.reg.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)


############################################################################################################
#Run Multi-Optima Direct Effect Model with Varying Effects
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_direct_ve.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_direct_ve.stan")

fit.reg.direct.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.direct.ve,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.direct.ve)

############################################################################################################
#Multilevel Multioptima Direct Effect Model with Varying Slopes
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-5
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
sigma.prior<-c(0,1)
############################################################################################################
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)

############################################################################################################
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_direct_mlm_ve.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_direct_mlm_ve.stan")

fit.reg.direct.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.direct.mlm.ve,pars = c("hl","vy","optima","beta","sigma"))
post<-extract(fit.reg.direct.mlm.ve)

############################################################################################################
#Multilevel Multi-optima Direct Effect Model with Varying Slopes - Non-centered
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-5
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
sigma.prior<-c(0,1)
############################################################################################################
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)

############################################################################################################
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_direct_mlm_ve_nc.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_direct_mlm_ve_nc.stan")

fit.reg.direct.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.direct.mlm.ve.nc,pars = c("hl","vy","optima","beta","sigma"))
post<-extract(fit.reg.direct.mlm.ve.nc)
############################################################################################################
#Multi-optima Adaptive Model with Varying Slopes
############################################################################################################
#Two regimes with adaptive trait and multiple slopes per optima but single alpha parameter
set.seed(10)
N<-50 #Number of species
#set.seed(1) #Set seed to get same random species each time

phy <- ape::keep.tip(tree.10K,sample(tree.10K$tip.label)[1:N]) 
phy<-ape::multi2di(phy)

l.tree<-max(branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree 

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

reg_tips<-trdata$dat$regimes
reg_tips<-as.numeric(as.factor(reg_tips))

print(reg.colors)
plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

#Phylogeny info
n<-length(trdata$phy$tip.label)

regimes_internal <-trdata$phy$node.label
regimes_tip <- trdata$dat$regimes
regimes <- concat.factor(regimes_tip, regimes_internal)
anc_maps<-"regimes"
lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

##############################################
hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.01 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));

vX0<-0
vY0 <- 0
Sxx<-10 #Look at effects


X<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
sigma2_x<-matrix(1,1,1)
Z_adaptive<-1
names(X)<-phy$tip.label
phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

optima_matrix<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
pred_X<-calc_adaptive_dmX(phy,a,X)
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes

mu<-matrix(NA,N,1)
for(i in 1:N){
  mu[i] = optima_matrix[i,]%*%optima+beta[reg_tips[i]]%*%pred_X[i]
}

n_reg<-length(unique(regimes))
V<-calc_adaptive_V(phy,a, sigma2_y, beta,  sigma2_x, Z_adaptive)
Y<-MASS::mvrnorm(n=1,mu,V)
##################################################################################################################
#Simulate errors 
Z_X_error<-1 #Number of X traits with error
X_error<-rep(0.01,N)
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-X+rnorm(N,0,0.01)
############################################################################################################
#Make trdata file
trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))

############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-5
optima.prior<-c(1.5,0.5)
beta.prior<-c(0,0.25)
############################################################################################################
#Test Blouch prep code - Regimes + Direct Efffect model
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)

############################################################################################################
#Plot of data
df<-data.frame(Y=dat$Y_obs,X=dat$X_obs,Regimes=regimes_tip)


#slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
slope.plot.1<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X,color=Regimes))+
  
  geom_abline(intercept=optima[1],slope=beta[1],alpha=0.5,linetype=2)+
  geom_abline(intercept=optima[2],slope=beta[2],alpha=0.5,linetype=2)+
  
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ylab("Y") + xlab("Adaptive trait")+
  ggsci::scale_color_npg()

slope.plot.1
############################################################################################################
#Run Stan code
############################################################################################################
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_adapt_ve.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_adapt_ve.stan")

fit.reg.adapt.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.adapt.ve,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.adapt.ve)

############################################################################################################
#Multilevel Multi-optima Adaptive Model with Varying Slopes
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-5
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
sigma.prior<-c(0,1)

############################################################################################################
#Test Blouch prep code - Regimes + Direct Efffect model
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)

############################################################################################################
#Run Stan code
############################################################################################################
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_adapt_mlm_ve.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_adapt_mlm_ve.stan")

fit.reg.adapt.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.adapt.mlm.ve,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.adapt.mlm.ve)


############################################################################################################
#Multilevel Multi-optima Adaptive Model with Varying Slopes - Non-centered
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-5
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
sigma.prior<-c(0,1)

############################################################################################################
#Test Blouch prep code - Regimes + Direct Efffect model
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)

############################################################################################################
#Run Stan code
############################################################################################################
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_adapt_mlm_ve_nc.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_adapt_mlm_ve_nc.stan")

fit.reg.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.adapt.mlm.ve.nc,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.adapt.mlm.ve.nc)

########################################################################################################
## Multi-optima Direct Effect and Adaptive models with Varying Effects - Multiple predictors
#Two regimes with 1 direct and 1 adaptive trait and two slopes per regime
########################################################################################################
set.seed(10)

N<-50 #Number of species
#set.seed(1) #Set seed to get same random species each time

phy <- ape::keep.tip(tree.10K,sample(tree.10K$tip.label)[1:N]) 
phy<-ape::multi2di(phy)

l.tree<-max(ape::branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree 

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

reg_tips<-trdata$dat$regimes
reg_tips<-as.numeric(as.factor(reg_tips))

reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(2)
print(reg.colors)
plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

#Phylogeny info
n<-length(trdata$phy$tip.label)

regimes_internal <-trdata$phy$node.label
regimes_tip <- trdata$dat$regimes
regimes <- concat.factor(regimes_tip, regimes_internal)
anc_maps<-"regimes"
lineages <- lapply(1:n, function(e) lineage.constructor(phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch
################################################################################################################
#Simulate X and Y data using generative model
################################################################################################################
hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.01 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));

vX0<-0
vY0 <- 0
sigma2_x<-matrix(1,1,1) #Variance of BM Process

Xa<-phytools::fastBM(phy,a=vX0,sig2=sigma2_x[1,1],internal=FALSE) #Simulate X BM variable on tree, with scaling 10
Xd<-phytools::fastBM(phy,a=vX0,sig2=sigma2_x[1,1],internal=FALSE) #Simulate X BM variable on tree, with scaling 10

names(Xa)<-phy$tip.label
names(Xd)<-phy$tip.label
phytools::phenogram(phy,Xd,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
phytools::phenogram(phy,Xa,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

Xs<-cbind(Xd,Xa)
sigma2_x<-matrix(1,1,1) #Variance of BM Process

Z_adaptive<-1
Z_direct<-1

optima_matrix<-weight.matrix(phy, a, lineages) #Slouch approach
pred_X<-calc_mixed_dmX(phy,a,Xs,Z_direct,Z_adaptive)
optima<-c(2,1)
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Two traits on columns, two regimes on vertical

mu<-matrix(NA,N,1)
for(i in 1:N){
  mu[i] = optima_matrix[i,]%*%optima+pred_X[i,]%*%t(beta[reg_tips[i],])
}

n_reg<-length(unique(regimes))
V<-calc_adaptive_V(phy,a, sigma2_y, beta[,2],  sigma2_x, Z_adaptive)
Y<-MASS::mvrnorm(n=1,mu,V)
###############################################################################################################
summary(lm(Y~Xs))
###############################################################################################################
Z_X_error<-2 #Number of X traits with error
X_error<-matrix(0.01,nrow=N,ncol=Z_X_error)
X_error<-data.frame(X_error)
names(X_error)<-c("Xd_error","Xa_error")
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-Xs+rnorm(N,0,0.01)
X_with_error<-data.frame(X_with_error)
names(X_with_error)<-c("Xd","Xa")
############################################################################################################
#Make trdata file
trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))
#names(trdata$dat)[6:7]<-c("Xd_error","Xa_error")
#We will use the helper function blouch.reg.adapt.prep() to setup the dat file for Stan. Here "Z_adaptive" is the number of predictors, and "regimes" is the name of the regime column in trdata$dat.
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-5
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
########################################################################################################
#Test Blouch prep code - Regimes + Direct Effect model
dat<-blouch.reg.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Plot of data
df<-data.frame(Y=dat$Y_obs,X=dat$X_obs,Regimes=regimes_tip)

#slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
slope.plot.1<-ggplot2::ggplot()+  
  ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=Xd,color=Regimes))+
  
  ggplot2::geom_abline(intercept=optima[1],slope=beta[1,1],alpha=0.5,linetype=2)+
  ggplot2::geom_abline(intercept=optima[2],slope=beta[2,1],alpha=0.5,linetype=2)+
  
  
  ggplot2::theme_bw()+
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ggplot2::ylab("Y") + ggplot2::xlab("Direct effect trait")+
  ggsci::scale_color_npg()

slope.plot.1

#slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
slope.plot.2<-ggplot2::ggplot()+  
  ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=Xa,color=Regimes))+
  
  ggplot2::geom_abline(intercept=optima[1],slope=beta[1,2],alpha=0.5,linetype=2)+
  ggplot2::geom_abline(intercept=optima[2],slope=beta[2,2],alpha=0.5,linetype=2)+
  
  ggplot2::theme_bw()+
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ggplot2::ylab("Y") + ggplot2::xlab("Adaptive trait")+
  ggsci::scale_color_npg()

slope.plot.2
########################################################################################################
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_direct_adapt_ve.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_direct_adapt_ve.stan")

fit.reg.direct.adapt.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.direct.adapt.ve,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.direct.adapt.ve)

########################################################################################################
#Multilevel Multi-Optima Direct Effect and Adaptive Model with Varying Effects
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-5
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
sigma.prior<-c(0,1)
########################################################################################################
#Test Blouch prep code - Regimes + Direct Effect model
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)

########################################################################################################
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_direct_adapt_mlm_ve.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_direct_adapt_mlm_ve.stan")

fit.reg.direct.adapt.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.direct.adapt.mlm.ve,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.direct.adapt.mlm.ve)

########################################################################################################
########################################################################################################
#Multilevel Multi-Optima Direct Effect and Adaptive Model with Varying Effects - Non-centered
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-5
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
sigma.prior<-c(0,1)
########################################################################################################
#Test Blouch prep code - Regimes + Direct Effect model
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)

########################################################################################################
setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/") #Macbook Pro
stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/Stan Models Milestones/Finished Versions - SBR2/blouchOU_reg_direct_adapt_mlm_ve_nc.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_direct_adapt_mlm_ve_nc.stan")

fit.reg.direct.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.reg.direct.adapt.mlm.ve.nc,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.direct.adapt.mlm.ve.nc)

########################################################################################################


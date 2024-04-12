#Direct effect model simulations
#Generative Code
#Script to simulate X an Y OU Data - for Milstone 1 Code
#04/26/2023

rm(list=ls())

calc_direct_V<-function(phy, sigma2_y, a){
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]] #Same as Cophenetic Distance matrix
  Vt<-sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) * exp(-a * tij)) #From Hansen (1997)
  #ta - time from root to tips, tij  - total time separating spcies
  #Variance of y/()
  return(Vt)
}
ts_fxn<-function(phy){ #Calculate t
  n<-length(phy$tip.label)
  mrca1 <- ape::mrca(phy) #Node numbers for MRCA of tips
  times <- ape::node.depth.edgelength(phy) #Time from root to tips of each node, starting with the tips
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(phy$tip.label, phy$tip.label)) #Matrix with time from root to MRCA of pairs of tips - pulls out values of times that correspond with node numbers - integers
  T.term <- times[1:n] #Times from root for tips
  tia <- times[1:n] - ta #Times from root to tips - times from root to MRCA = times from MRCA to tips
  tja <- t(tia) #Transpose of the times from MRCA to tips matrix
  #tij <- tja + tia #Sum of matrix and its transpose - total time separating species
  tij<-cophenetic(phy)
  #return(list(ta,tia,tja,tij,T.term))
  return(list(ta,tij,T.term,tja))
}

########################################################################################################
#Basic Setup
library(geiger)
library(MASS)
library(phytools)
library(ggplot2)
#Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997) Original Stan

tree.10K<-read.tree("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Original Submission/Blouch Testing/Phylogeny/10KPrimateTree.tre")
#tree.10K<-read.tree("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/Original Submission/Blouch Testing/Phylogeny/10KPrimateTree.tre")
N<-50 #Number of species
set.seed(10) #Set seed to get same random species each time

phy <- keep.tip(tree.10K,sample(tree.10K$tip.label)[1:N]) 
phy<-multi2di(phy)

l.tree<-max(branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree 

tip.label<-phy$tip.label
#Dmat<-cophenetic(phy) #Time separating tips, same as tij matrix in Slouch/Blouch code
#Dmat<-Dmat[tip.label,tip.label]/max(Dmat)




########################################################################################################
#Direct Effect Model
#Setup parameters
Z<-1 #Number of traits
hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.01 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));

vX0<-0
vY0 <- 0

Sxx<-10 #Look at effects
#var_y_anc<-1
ts<-ts_fxn(phy)
ta<-ts[[1]]
tij<-ts[[2]]

V<-calc_direct_V(phy,sigma2_y,a)
X<-fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

alpha<-2 #Intecept
beta<-0.25 #Slope

mu<-alpha+X*beta #Simulate mu for Y
#Simulate direct effect Y trait
Y<-mvrnorm(n=1,mu,V)

df<-data.frame(Y=Y,X=X)
names(df)<-c("Y","X")

ggplot(data=df,aes(x=X,y=Y))+
  geom_point()
summary(lm(Y~X,df))

########################################################################################################
#Direct effect model w/o ME following Hansen et al. (2008)
dat<-list(N=N,Z=Z,Y_obs=Y,X_obs=matrix(X,nrow=N,ncol=Z),ta=ta,tij=tij)


########################################################################################################
#Direct effect model w/o ME - Statistical Rethinking Approach to V
dat<-list(N=N,Z=Z,Y=Y,X=matrix(X,nrow=N,ncol=Z),Dmat=Dmat)


#Slouch testing
library(slouch)
slouch.results<-slouch.fit(phy = phy,
                           species = phy$tip.label,
                           response = dat$Y,
                           direct.cov = X)

summary(slouch.results)


########################################################################################################
#Simulate errors - original Hansen setup
X_error<-matrix(0.01,nrow=N,ncol=Z)
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-X+rnorm(N,0,0.01)

#Direct effect model w/ Statistical Rethinking ME Correction
dat<-list(N=N,Z=Z,Y_obs=Y_with_error,X_obs=matrix(X_with_error,nrow=N,ncol=Z),Y_error=Y_error,X_error=matrix(X_error,nrow=N,ncol=Z),ta=ta,tij=tij)

########################################################################################################
#Simulate multiple X traits - correlated
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
alpha<-2
beta<-c(0.35,0.1) #Slope
mu<-alpha+Xs%*%beta #Simulate mu for Y

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

#Without ME
dat<-list(N=N,Z=Z,Y_obs=Y,X_obs=matrix(Xs,nrow=N,ncol=Z),ta=ta,tij=tij)

########################################################################################################
#Simulate errors with multiple traits - original Hansen setup
X_error<-matrix(0.01,nrow=N,ncol=Z)
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-apply(Xs,2,function(X){X+rnorm(N,0,0.01)})
#X_with_error<-X+rnorm(N,0,0.1)

#Direct effect model w/ Statistical Rethinking ME Correction
dat<-list(N=N,Z=Z,Y_obs=Y_with_error,X_obs=matrix(X_with_error,nrow=N,ncol=Z),Y_error=Y_error,X_error=matrix(X_error,nrow=N,ncol=Z),ta=ta,tij=tij)

########################################################################################################
#Simulate multiple X traits - correlated
Z<-3
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
vcv<-matrix(c(1,0,0,0,1,0,0,0,1),3,3) #Correlation between traits
Xs<-sim.corrs(phy,vcv) #Simulated correlated BM Xs
alpha<-2
beta<-c(0.5,0.35,0.1) #Slope
mu<-alpha+Xs%*%beta #Simulate mu for Y

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

ggplot(data=df,aes(x=X.3,y=Y))+
  geom_point()
summary(lm(Y~X.3,df))

#Without ME
dat<-list(N=N,Z=Z,Y_obs=Y,X_obs=matrix(Xs,nrow=N,ncol=Z),ta=ta,tij=tij)

########################################################################################################
#Simulate errors with multiple traits - original Hansen setup
X_error<-matrix(0.01,nrow=N,ncol=Z)
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-apply(Xs,2,function(X){X+rnorm(N,0,0.01)})
#X_with_error<-X+rnorm(N,0,0.1)

#Direct effect model w/ Statistical Rethinking ME Correction
dat<-list(N=N,Z=Z,Y_obs=Y_with_error,X_obs=matrix(X_with_error,nrow=N,ncol=Z),Y_error=Y_error,X_error=matrix(X_error,nrow=N,ncol=Z),ta=ta,tij=tij)




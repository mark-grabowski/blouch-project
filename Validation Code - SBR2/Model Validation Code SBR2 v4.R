#Blouch Model Validation Code - 2024
#SBR2 04/09/24
rm(list=ls())
########################################################################################################
#Complete code to make simulated data, setup, run model validation for SI
#Also need to load blouch.prep.R code to aid in simulations (calculate Vs, etc.), and set.converge.regimes.R
#Using versions in blouch package
#Code tests each function external from the package before it is placed within the package
#Also makes all plots as PDFs - for SI
########################################################################################################
#Basic Setup
#library(geiger)
#library(MASS)
#library(phytools)
#library(ggplot2)
#library(rstan)
#library(rethinking)
#library(blouch)
########################################################################################################
#For execution on a local, multicore CPU with excess RAM we recommend calling
options(mc.cores = parallel::detectCores())
#options(mc.cores = 8)
rstan::rstan_options(auto_write = TRUE)

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
arch <- ifelse(R.version$arch == "aarch64", "arm64", "x86_64")
cat(paste("\nCXX14FLAGS += -O3 -mtune=native -arch", arch, "-ftemplate-depth-256"),
    file = M, sep = "\n", append = FALSE)
########################################################################################################
#source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/bloutch-testing/R Setup Code/blouch.prep.SBR2.R")
tree.10K<-ape::read.tree("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/My Drive/Current/Blouch/Blouch project/blouch project - for ms - SBR2/Sharable Data/10KPrimateTree.tre")
source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/R/set.converge.regimes.R")
source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/R/blouch.prep.R")
source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/R/simulate.data.R")
source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/R/simulation.helpers.R")
source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/R/create.plots.R")
############################################################################################################
set.seed(10) #Set seed to get same random species each time
N<-50
phy <- ape::keep.tip(tree.10K,sample(tree.10K$tip.label)[1:N])
phy<-ape::multi2di(phy)
l.tree<-max(ape::branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree
########################################################################################################
#Direct Effect Model
########################################################################################################
#N<-50 #Number of species
Z<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-2
beta<-0.25
trdata<-sim.direct.data(phy,N,Z,hl,vy,Sxx,optima,beta)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(2,0.25) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
direct.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch prep code - direct effect model - blouch.direct.prep()
dat<-blouch.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Direct effect model - Running from R Package Version
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_direct.stan")
stan_model <- rstan::stan_model("blouchOU_direct.stan")
fit.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2, cores=2,iter =2000)
print(fit.direct,pars = c("hl","vy","optima","beta"))

post<-rstan::extract(fit.direct)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-direct.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot, ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Direct Effect Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig1.pdf", plot = fig, width=7, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
#N<-50 #Number of species
Z<-1
hl<-0.25
vy<-0.01
Sxx<-10
optima<-2
beta<-0.25
trdata<-sim.direct.data(phy,N,Z,hl,vy,Sxx,optima,beta)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(2,0.25) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
direct.prior.plot.code(trdata,optima.prior,beta.prior)

############################################################################################################
#Blouch prep code - direct effect model - blouch.direct.prep()
dat<-blouch.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Direct effect model - Running from R Package Version
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_direct.stan")
stan_model <- rstan::stan_model("blouchOU_direct.stan")
fit.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2, cores=2,iter =2000)
post<-rstan::extract(fit.direct)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-direct.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot, ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Direct Effect Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig2.pdf", plot = fig, width=7, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z<-1
hl<-0.75
vy<-0.01
Sxx<-10
optima<-2
beta<-0.25
trdata<-sim.direct.data(phy,N,Z,hl,vy,Sxx,optima,beta)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(2,0.25) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
direct.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch prep code - direct effect model - blouch.direct.prep()
dat<-blouch.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Direct effect model - Running from R Package Version
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_direct.stan")
stan_model <- rstan::stan_model("blouchOU_direct.stan")
fit.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2, cores=2,iter =2000)
post<-rstan::extract(fit.direct)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-direct.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot, ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Direct Effect Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig3.pdf", plot = fig, width=7, height=7 )
############################################################################################################
########################################################################################################
#Adaptation Models
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
#N<-50 #Number of species
Z<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-2
beta<-0.25
trdata<-sim.adaptive.data(phy,N,Z,hl,vy,Sxx,optima,beta)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(2,0.25) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch prep code
dat<-blouch.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Adaptive model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_adapt.stan")
stan_model <- rstan::stan_model("blouchOU_adapt.stan")
fit.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2, cores=2,iter =2000)
print(fit.adapt,pars = c("hl","vy","optima","beta"))

post<-rstan::extract(fit.adapt)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot, ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Adaptation Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig4.pdf", plot = fig, width=7, height=7 )

########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z<-1
hl<-0.25
vy<-0.01
Sxx<-10
optima<-2
beta<-0.25
trdata<-sim.adaptive.data(phy,N,Z,hl,vy,Sxx,optima,beta)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(2,0.25) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
direct.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch prep code
dat<-blouch.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Adaptive model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_adapt.stan")
stan_model <- rstan::stan_model("blouchOU_adapt.stan")
fit.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2, cores=2,iter =2000)
print(fit.adapt,pars = c("hl","vy","optima","beta"))

post<-rstan::extract(fit.adapt)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot, ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Adaptation Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig5.pdf", plot = fig, width=7, height=7 )

########################################################################################################
#set.seed(10) #Set seed to get same random species each time
#N<-50 #Number of species
Z<-1
hl<-0.75
vy<-0.01
Sxx<-10
optima<-2
beta<-0.25
trdata<-sim.adaptive.data(phy,N,Z,hl,vy,Sxx,optima,beta)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(2,0.25) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch prep code
dat<-blouch.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Adaptive model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_adapt.stan")
stan_model <- rstan::stan_model("blouchOU_adapt.stan")
fit.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2, cores=2,iter =2000)
print(fit.adapt,pars = c("hl","vy","optima","beta"))

post<-rstan::extract(fit.adapt)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot, ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Adaptation Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig6.pdf", plot = fig, width=7, height=7 )
############################################################################################################
########################################################################################################
#Direct Effect + Adaptation Model
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
#N<-50 #Number of species
Z_direct<-1
Z_adaptive<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-2
beta<-c(0.35,0.25)
trdata<-sim.direct.adaptive.data(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(2,0.2) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#View priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch prep code
dat<-blouch.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,hl.prior,vy.prior,optima.prior,beta.prior)
############################################################################################################
#Direct Effect and Adaptation Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_direct_adapt.stan")
stan_model <- rstan::stan_model("blouchOU_direct_adapt.stan")
fit.direct.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.direct.adapt)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-direct.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Direct Effect + Adaptation Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig7.pdf", plot = fig, width=10.5, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
#N<-50 #Number of species
Z_direct<-1
Z_adaptive<-1
hl<-0.25
vy<-0.01
Sxx<-10
optima<-2
beta<-c(0.35,0.25)
trdata<-sim.direct.adaptive.data(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(2,0.2) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#View priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch prep code
dat<-blouch.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,hl.prior,vy.prior,optima.prior,beta.prior)
############################################################################################################
#Direct Effect and Adaptation Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_direct_adapt.stan")
stan_model <- rstan::stan_model("blouchOU_direct_adapt.stan")
fit.direct.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.direct.adapt)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-direct.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Direct Effect + Adaptation Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig8.pdf", plot = fig, width=10.5, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
Z_adaptive<-1
hl<-0.75
vy<-0.01
Sxx<-10
optima<-2
beta<-c(0.35,0.25)
trdata<-sim.direct.adaptive.data(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(2,0.2) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#View priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch prep code
dat<-blouch.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,hl.prior,vy.prior,optima.prior,beta.prior)
############################################################################################################
#Direct Effect and Adaptation Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_direct_adapt.stan")
stan_model <- rstan::stan_model("blouchOU_direct_adapt.stan")
fit.direct.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.direct.adapt)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-direct.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Direct Effect + Adaptation Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig9.pdf", plot = fig, width=10.5, height=7 )

############################################################################################################
########################################################################################################
#Multi-Optima Model
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
#N<-50 #Number of species
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(0.5,0.25)
shifts<-c(84)
############################################################################################################
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
trdata<-sim.reg.data(phy,N,hl,vy,Sxx,optima,shifts)
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
optima.prior.plot.code(trdata,optima.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior)
############################################################################################################
#Multi-Optima Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg.stan")
stan_model <- rstan::stan_model("blouchOU_reg.stan")
fit.reg<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-optima.prior.post.plot.code(trdata,optima.prior,post,optima)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig10.pdf", plot = fig, width=7, height=7 )

########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
hl<-0.25
vy<-0.01
Sxx<-10
optima<-c(0.5,0.25)
shifts<-c(84)
############################################################################################################
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
trdata<-sim.reg.data(phy,N,hl,vy,Sxx,optima,shifts)
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
optima.prior.plot.code(trdata,optima.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior)
############################################################################################################
#Multi-Optima Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg.stan")
stan_model <- rstan::stan_model("blouchOU_reg.stan")
fit.reg<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-optima.prior.post.plot.code(trdata,optima.prior,post,optima)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Model, hl=",hl,sep=""))

ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig11.pdf", plot = fig, width=7, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
hl<-0.75
vy<-0.01
Sxx<-10
optima<-c(0.5,0.25)
shifts<-c(84)
############################################################################################################
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
trdata<-sim.reg.data(phy,N,hl,vy,Sxx,optima,shifts)
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
optima.prior.plot.code(trdata,optima.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior)
############################################################################################################
#Multi-Optima Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg.stan")
stan_model <- rstan::stan_model("blouchOU_reg.stan")
fit.reg<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-optima.prior.post.plot.code(trdata,optima.prior,post,optima)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig12.pdf", plot = fig, width=7, height=7 )
########################################################################################################
########################################################################################################
#Multilevel Multi-Optima Model - Varying Intercepts
############################################################################################################
#set.seed(10) #Set seed to get same random species each time
#N<-50 #Number of species
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(0.5,0.25)
shifts<-c(84)
############################################################################################################
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
sigma.prior<-c(0,1)
trdata<-sim.reg.data(phy,N,hl,vy,Sxx,optima,shifts)
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
optima.prior.plot.code(trdata,optima.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.mlm.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior,sigma.prior)
############################################################################################################
#Multilevel Multi-Optima Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_mlm_vi.stan")
stan_model <- rstan::stan_model("blouchOU_reg_mlm_vi.stan")
fit.reg.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.mlm.vi)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-optima.prior.post.plot.code(trdata,optima.prior,post,optima)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig13.pdf", plot = fig, width=7, height=7 )

############################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
hl<-0.25
vy<-0.01
Sxx<-10
optima<-c(0.5,0.25)
shifts<-c(84)
############################################################################################################
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
sigma.prior<-c(0,1)
trdata<-sim.reg.data(phy,N,hl,vy,Sxx,optima,shifts)
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
optima.prior.plot.code(trdata,optima.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.mlm.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior,sigma.prior)
############################################################################################################
#Multilevel Multi-Optima Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_mlm_vi.stan")
stan_model <- rstan::stan_model("blouchOU_reg_mlm_vi.stan")
fit.reg.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.mlm.vi)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-optima.prior.post.plot.code(trdata,optima.prior,post,optima)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig14.pdf", plot = fig, width=7, height=7 )
############################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
hl<-0.75
vy<-0.01
Sxx<-10
optima<-c(0.5,0.25)
shifts<-c(84)
############################################################################################################
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
sigma.prior<-c(0,1)
trdata<-sim.reg.data(phy,N,hl,vy,Sxx,optima,shifts)
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
optima.prior.plot.code(trdata,optima.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.mlm.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior,sigma.prior)
############################################################################################################
#Multilevel Multi-Optima Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_mlm_vi.stan")
stan_model <- rstan::stan_model("blouchOU_reg_mlm_vi.stan")
fit.reg.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.mlm.vi)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-optima.prior.post.plot.code(trdata,optima.prior,post,optima)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig15.pdf", plot = fig, width=7, height=7 )
########################################################################################################
############################################################################################################
#Multilevel Multi-Optima Model- Varying Intercepts - Non-Centered
############################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(0.5,0.25)
shifts<-c(84)
############################################################################################################
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
sigma.prior<-c(0,1)
trdata<-sim.reg.data(phy,N,hl,vy,Sxx,optima,shifts)
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
optima.prior.plot.code(trdata,optima.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.mlm.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior,sigma.prior)
############################################################################################################
#Multilevel Multi-Optima Model - Non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_mlm_vi_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_mlm_vi_nc.stan")
fit.reg.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.mlm.vi.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-optima.prior.post.plot.code(trdata,optima.prior,post,optima)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Model\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig16.pdf", plot = fig, width=7, height=7 )

############################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
hl<-0.25
vy<-0.01
Sxx<-10
optima<-c(0.5,0.25)
shifts<-c(84)
############################################################################################################
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
sigma.prior<-c(0,1)
trdata<-sim.reg.data(phy,N,hl,vy,Sxx,optima,shifts)
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
optima.prior.plot.code(trdata,optima.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.mlm.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior,sigma.prior)
############################################################################################################
#Multilevel Multi-Optima Model - Non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_mlm_vi_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_mlm_vi_nc.stan")
fit.reg.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.mlm.vi.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-optima.prior.post.plot.code(trdata,optima.prior,post,optima)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Model\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig17.pdf", plot = fig, width=7, height=7 )

############################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
hl<-0.75
vy<-0.01
Sxx<-10
optima<-c(0.5,0.25)
shifts<-c(84)
############################################################################################################
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
sigma.prior<-c(0,1)
trdata<-sim.reg.data(phy,N,hl,vy,Sxx,optima,shifts)
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
optima.prior.plot.code(trdata,optima.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.mlm.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior,sigma.prior)
############################################################################################################
#Multilevel Multi-Optima Model - Non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_mlm_vi_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_mlm_vi_nc.stan")
fit.reg.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.mlm.vi.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-optima.prior.post.plot.code(trdata,optima.prior,post,optima)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Model\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig18.pdf", plot = fig, width=7, height=7 )

############################################################################################################
########################################################################################################
#Multi-Optima Direct Effect Model
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
#N<-50 #Number of species
Z_direct<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.25)
shifts<-c(84)
trdata<-sim.reg.direct.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
############################################################################################################
#Multi-Optima Direct Effect Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct.stan")
fit.reg.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-reg.direct.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig19.pdf", plot = fig, width=7, height=7 )

########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
hl<-0.25
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.direct.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
############################################################################################################
#Multi-Optima Direct Effect Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct.stan")
fit.reg.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-reg.direct.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig20.pdf", plot = fig, width=7, height=7 )

########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
hl<-0.75
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.direct.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
############################################################################################################
#Multi-Optima Direct Effect Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct.stan")
fit.reg.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-reg.direct.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig21.pdf", plot = fig, width=7, height=7 )

############################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
#N<-50 #Number of species
Z_direct<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.direct.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
############################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_vi.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_vi.stan")
fit.reg.direct.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.mlm.vi)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-reg.direct.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect Model - Varying Intercepts, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig22.pdf", plot = fig, width=7, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
hl<-0.25
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.direct.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
############################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_vi.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_vi.stan")
fit.reg.direct.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.mlm.vi)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-reg.direct.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect Model - Varying Intercepts, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig23.pdf", plot = fig, width=7, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
hl<-0.75
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.direct.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
############################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_vi.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_vi.stan")
fit.reg.direct.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.mlm.vi)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-reg.direct.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect Model - Varying Intercepts, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig24.pdf", plot = fig, width=7, height=7 )

############################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts - Non-centered
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.direct.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts - Non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_vi_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_vi_nc.stan")
fit.reg.direct.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.mlm.vi.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-reg.direct.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect Model - Varying Intercepts\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig25.pdf", plot = fig, width=7, height=7 )

########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
hl<-0.25
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.direct.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts - Non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_vi_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_vi_nc.stan")
fit.reg.direct.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.mlm.vi.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-reg.direct.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect Model - Varying Intercepts\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig26.pdf", plot = fig, width=7, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
hl<-0.75
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.direct.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts - Non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_vi_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_vi_nc.stan")
fit.reg.direct.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.mlm.vi.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-reg.direct.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect Model - Varying Intercepts\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig27.pdf", plot = fig, width=7, height=7 )
############################################################################################################
########################################################################################################
#Multi-Optima Direct Effect Model - Varying Effects
########################################################################################################
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10 #Look at effects
shifts<-c(84)
trdata<-sim.reg.direct.ve.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
dat<-blouch.reg.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_ve.stan")
fit.reg.direct.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig28.pdf", plot = fig, width=7, height=7 )
########################################################################################################
hl<-0.25 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10 #Look at effects
shifts<-c(84)
trdata<-sim.reg.direct.ve.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
dat<-blouch.reg.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_ve.stan")
fit.reg.direct.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig29.pdf", plot = fig, width=7, height=7 )
########################################################################################################
hl<-0.75 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10 #Look at effects
shifts<-c(84)
trdata<-sim.reg.direct.ve.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
dat<-blouch.reg.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_ve.stan")
fit.reg.direct.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig30.pdf", plot = fig, width=7, height=7 )
########################################################################################################
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects
########################################################################################################
N<-50
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10 #Look at effects
shifts<-c(84)
trdata<-sim.reg.direct.ve.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_ve.stan")
fit.reg.direct.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.mlm.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig31.pdf", plot = fig, width=7, height=7 )
############################################################################################################
N<-50
hl<-0.25 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10 #Look at effects
shifts<-c(84)
trdata<-sim.reg.direct.ve.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_ve.stan")
fit.reg.direct.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.mlm.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig32.pdf", plot = fig, width=7, height=7 )
############################################################################################################
N<-50
hl<-0.75 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10 #Look at effects
shifts<-c(84)
trdata<-sim.reg.direct.ve.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_ve.stan")
fit.reg.direct.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.mlm.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig33.pdf", plot = fig, width=7, height=7 )
########################################################################################################
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects - Non-centered
########################################################################################################
N<-50
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10 #Look at effects
shifts<-c(84)
trdata<-sim.reg.direct.ve.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects - Non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_ve_nc.stan")
fit.reg.direct.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.mlm.ve.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig34.pdf", plot = fig, width=7, height=7 )
########################################################################################################
########################################################################################################
N<-50
hl<-0.25 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10 #Look at effects
shifts<-c(84)
trdata<-sim.reg.direct.ve.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_ve_nc.stan")
fit.reg.direct.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.mlm.ve.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig35.pdf", plot = fig, width=7, height=7 )
########################################################################################################
########################################################################################################
N<-50
hl<-0.75 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10 #Look at effects
shifts<-c(84)
trdata<-sim.reg.direct.ve.data(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25)
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_ve_nc.stan")
fit.reg.direct.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.mlm.ve.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig36.pdf", plot = fig, width=7, height=7 )

########################################################################################################
########################################################################################################
#Multi-Optima Adaptation Model
########################################################################################################
set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_adapt<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.adapt.data(phy,N,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.25) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Adaptation Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt.stan")
fit.reg.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Adaptation Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig37.pdf", plot = fig, width=7, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_adapt<-1
hl<-0.25
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.adapt.data(phy,N,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Adaptation Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt.stan")
fit.reg.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Adaptation Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig38.pdf", plot = fig, width=7, height=7 )

########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_adapt<-1
hl<-0.75
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.adapt.data(phy,N,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Adaptation Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt.stan")
fit.reg.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Adaptation Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig39.pdf", plot = fig, width=7, height=7 )

########################################################################################################
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Intercepts
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_adapt<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.adapt.data(phy,N,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_vi.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_vi.stan")
fit.reg.adapt.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.mlm.vi)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptation Model - Varying Intercepts, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig40.pdf", plot = fig, width=7, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_adapt<-1
hl<-0.25
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.adapt.data(phy,N,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_vi.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_vi.stan")
fit.reg.adapt.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.mlm.vi)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptation Model - Varying Intercepts, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig41.pdf", plot = fig, width=7, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_adapt<-1
hl<-0.75
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.adapt.data(phy,N,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_vi.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_vi.stan")
fit.reg.adapt.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.mlm.vi)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptation Model - Varying Intercepts, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig42.pdf", plot = fig, width=7, height=7 )

########################################################################################################
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Intercepts - Non-centered
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_adapt<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.adapt.data(phy,N,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_vi_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_vi_nc.stan")
fit.reg.adapt.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.mlm.vi.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptation Model - Varying Intercepts\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig43.pdf", plot = fig, width=7, height=7 )

########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_adapt<-1
hl<-0.25
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.adapt.data(phy,N,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_vi_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_vi_nc.stan")
fit.reg.adapt.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.mlm.vi.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptation Model - Varying Intercepts\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig44.pdf", plot = fig, width=7, height=7 )

########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_adapt<-1
hl<-0.75
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-0.25
shifts<-c(84)
trdata<-sim.reg.adapt.data(phy,N,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_vi_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_vi_nc.stan")
fit.reg.adapt.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.mlm.vi.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptation Model - Varying Intercepts\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig45.pdf", plot = fig, width=7, height=7 )

########################################################################################################
########################################################################################################
## Multi-Optima Adaptation Model - Varying Effects
########################################################################################################
N<-50
Z_adaptive<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.adapt.ve.data(phy,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
########################################################################################################
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_ve.stan")
fit.reg.adapt.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig46.pdf", plot = fig, width=7, height=7 )
########################################################################################################
N<-50
Z_adaptive<-1
hl<-0.25 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.adapt.ve.data(phy,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
########################################################################################################
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_ve.stan")
fit.reg.adapt.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig47.pdf", plot = fig, width=7, height=7 )
########################################################################################################
N<-50
Z_adaptive<-1
hl<-0.75 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.adapt.ve.data(phy,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
########################################################################################################
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_ve.stan")
fit.reg.adapt.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig48.pdf", plot = fig, width=7, height=7 )
########################################################################################################
############################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Effects
############################################################################################################
N<-50
Z_adaptive<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.adapt.ve.data(phy,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve.stan")
fit.reg.adapt.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.mlm.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig49.pdf", plot = fig, width=7, height=7 )
############################################################################################################
N<-50
Z_adaptive<-1
hl<-0.25 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.adapt.ve.data(phy,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve.stan")
fit.reg.adapt.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.mlm.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig50.pdf", plot = fig, width=7, height=7 )
############################################################################################################
N<-50
Z_adaptive<-1
hl<-0.75 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.adapt.ve.data(phy,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve.stan")
fit.reg.adapt.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.mlm.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig51.pdf", plot = fig, width=7, height=7 )
############################################################################################################
############################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Effects - Non-centered
############################################################################################################

N<-50
Z_adaptive<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.adapt.ve.data(phy,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve_nc.stan")
fit.reg.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.mlm.ve.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig52.pdf", plot = fig, width=7, height=7 )
############################################################################################################
N<-50
Z_adaptive<-1
hl<-0.25 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.adapt.ve.data(phy,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve_nc.stan")
fit.reg.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.mlm.ve.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig53.pdf", plot = fig, width=7, height=7 )
############################################################################################################
N<-50
Z_adaptive<-1
hl<-0.75 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.adapt.ve.data(phy,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve_nc.stan")
fit.reg.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.adapt.mlm.ve.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig54.pdf", plot = fig, width=7, height=7 )

############################################################################################################
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
Z_adapt<-1
hl<-0.1
vy<-0.01
Sxx<-1
optima<-c(2,1)
beta<-c(0.35,0.25)
#beta<-c(2,1,0.35,0.25) #Two Optima/Two Slopes
shifts<-c(84)
trdata<-sim.reg.direct.adapt.data(phy,N,Z_direct,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
dat<-blouch.reg.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt.stan")
fit.reg.direct.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect + Adaptation Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig55.pdf", plot = fig, width=10.5, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
Z_adapt<-1
hl<-0.25
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.35,0.25)
#beta<-c(2,1,0.35,0.25) #Two Optima/Two Slopes
shifts<-c(84)
trdata<-sim.reg.direct.adapt.data(phy,N,Z_direct,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
dat<-blouch.reg.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt.stan")
fit.reg.direct.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect + Adaptation Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig56.pdf", plot = fig, width=10.5, height=7 )
########################################################################################################
#set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
Z_adapt<-1
hl<-0.75
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.35,0.25)
#beta<-c(2,1,0.35,0.25) #Two Optima/Two Slopes
shifts<-c(84)
trdata<-sim.reg.direct.adapt.data(phy,N,Z_direct,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
dat<-blouch.reg.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt.stan")
fit.reg.direct.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect + Adaptation Model, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig57.pdf", plot = fig, width=10.5, height=7 )

############################################################################################################
############################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts
############################################################################################################
N<-50 #Number of species
Z_direct<-1
Z_adapt<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.35,0.25) #Two Optima/Two Slopes
shifts<-c(84)
trdata<-sim.reg.direct.adapt.data(phy,N,Z_direct,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_vi.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_vi.stan")
fit.reg.direct.adapt.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.mlm.vi)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig58.pdf", plot = fig, width=10.5, height=7 )
############################################################################################################
N<-50 #Number of species
Z_direct<-1
Z_adapt<-1
hl<-0.25
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.35,0.25) #Two Optima/Two Slopes
shifts<-c(84)
trdata<-sim.reg.direct.adapt.data(phy,N,Z_direct,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_vi.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_vi.stan")
fit.reg.direct.adapt.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.mlm.vi)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig59.pdf", plot = fig, width=10.5, height=7 )
############################################################################################################
N<-50 #Number of species
Z_direct<-1
Z_adapt<-1
hl<-0.75
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.35,0.25) #Two Optima/Two Slopes
shifts<-c(84)
trdata<-sim.reg.direct.adapt.data(phy,N,Z_direct,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_vi.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_vi.stan")
fit.reg.direct.adapt.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.mlm.vi)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig60.pdf", plot = fig, width=10.5, height=7 )

############################################################################################################
############################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts - Non-cetntered
############################################################################################################
N<-50 #Number of species
Z_direct<-1
Z_adapt<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.35,0.25) #Two Optima/Two Slopes
shifts<-c(84)
trdata<-sim.reg.direct.adapt.data(phy,N,Z_direct,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_vi_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_vi_nc.stan")
fit.reg.direct.adapt.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.mlm.vi.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig61.pdf", plot = fig, width=10.5, height=7 )
############################################################################################################
N<-50 #Number of species
Z_direct<-1
Z_adapt<-1
hl<-0.25
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.35,0.25) #Two Optima/Two Slopes
shifts<-c(84)
trdata<-sim.reg.direct.adapt.data(phy,N,Z_direct,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_vi_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_vi_nc.stan")
fit.reg.direct.adapt.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.mlm.vi.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig62.pdf", plot = fig, width=10.5, height=7 )
############################################################################################################
N<-50 #Number of species
Z_direct<-1
Z_adapt<-1
hl<-0.75
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.35,0.25) #Two Optima/Two Slopes
shifts<-c(84)
trdata<-sim.reg.direct.adapt.data(phy,N,Z_direct,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_vi_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_vi_nc.stan")
fit.reg.direct.adapt.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.mlm.vi.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig63.pdf", plot = fig, width=10.5, height=7 )


############################################################################################################
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects
########################################################################################################
N<-50
Z_direct<-1
Z_adaptive<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1) #Regime 1, 2, 3, 4
#optima<-c(2,1.5,1,0.5) #Regime 1, 2, 3, 4
#beta<-c(0.35,0.25) #Two Optima/Two Slopes
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Direct effect - column 1, Adaptive - column 2, slopes for regimes on rows
Sxx<-10
shifts<-c(84)
#shifts<-c(94,54,72) #Location of nodes with regime shifts
trdata<-sim.reg.direct.adapt.ve.data(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
########################################################################################################
dat<-blouch.reg.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_ve.stan")
fit.reg.direct.adapt.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect + Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig64.pdf", plot = fig, width=7, height=7 )
############################################################################################################
N<-50
Z_direct<-1
Z_adaptive<-1
hl<-0.25 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1.5)
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Direct effect - column 1, Adaptive - column 2, slopes for regimes on rows
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.direct.adapt.ve.data(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
########################################################################################################
dat<-blouch.reg.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_ve.stan")
fit.reg.direct.adapt.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect + Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig65.pdf", plot = fig, width=7, height=7 )
############################################################################################################
N<-50
Z_direct<-1
Z_adaptive<-1
hl<-0.75 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1.5)
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Direct effect - column 1, Adaptive - column 2, slopes for regimes on rows
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.direct.adapt.ve.data(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
#reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior,optima,beta)
########################################################################################################
dat<-blouch.reg.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_ve.stan")
fit.reg.direct.adapt.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect + Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig66.pdf", plot = fig, width=7, height=7 )
############################################################################################################
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Effects
########################################################################################################
Z_direct<-1
Z_adaptive<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1.5)
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Direct effect - column 1, Adaptive - column 2, slopes for regimes on rows
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.direct.adapt.ve.data(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
########################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_ve.stan")
fit.reg.direct.adapt.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.mlm.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect + Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig67.pdf", plot = fig, width=7, height=7 )
############################################################################################################
Z_direct<-1
Z_adaptive<-1
hl<-0.25 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1.5)
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Direct effect - column 1, Adaptive - column 2, slopes for regimes on rows
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.direct.adapt.ve.data(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
########################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_ve.stan")
fit.reg.direct.adapt.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.mlm.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect + Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig68.pdf", plot = fig, width=7, height=7 )
########################################################################################################
Z_direct<-1
Z_adaptive<-1
hl<-0.75 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1.5)
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Direct effect - column 1, Adaptive - column 2, slopes for regimes on rows
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.direct.adapt.ve.data(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
########################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_ve.stan")
fit.reg.direct.adapt.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.mlm.ve)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect + Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig69.pdf", plot = fig, width=7, height=7 )
############################################################################################################
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Effects - Non-centered
########################################################################################################
Z_direct<-1
Z_adaptive<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1.5)
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Direct effect - column 1, Adaptive - column 2, slopes for regimes on rows
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.direct.adapt.ve.data(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
########################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_ve_nc.stan")
fit.reg.direct.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.mlm.ve.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect + Adaptive Model - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig70.pdf", plot = fig, width=7, height=7 )
############################################################################################################
Z_direct<-1
Z_adaptive<-1
hl<-0.25 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1.5)
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Direct effect - column 1, Adaptive - column 2, slopes for regimes on rows
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.direct.adapt.ve.data(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
########################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_ve_nc.stan")
fit.reg.direct.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.mlm.ve.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect + Adaptive Model - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig71.pdf", plot = fig, width=7, height=7 )
########################################################################################################
Z_direct<-1
Z_adaptive<-1
hl<-0.75 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1.5)
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Direct effect - column 1, Adaptive - column 2, slopes for regimes on rows
Sxx<-10
shifts<-c(84)
trdata<-sim.reg.direct.adapt.ve.data(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(1.5,0.5) #Informed by linear model
beta.prior<-c(0,0.25) #Informed by linear model
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
#reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior,optima,beta)
sigma.prior.plot.code(sigma.prior)
########################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_ve_nc.stan")
fit.reg.direct.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.reg.direct.adapt.mlm.ve.nc)
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.direct.adapt.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,optima,beta)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot[[1]],regression.plot[[2]],ncol=3,nrow=2, labels = c("A)","B)","C)","D)","E)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect + Adaptive Model - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Validation/Valid_Fig72.pdf", plot = fig, width=7, height=7 )
############################################################################################################
############################################################################################################
#Multi-Optima Direct Effect Model - Varying Effects - Multiple Traits
############################################################################################################

############################################################################################################
#Multi-Optima Adaptation Model - Varying Effects - Multiple Traits
############################################################################################################














#MLM
########################################################################################################
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects
########################################################################################################
set.seed(10)
N<-100
Z_direct<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
#optima<-c(5,4.5,4,3.5,3,2.5,2,1.5,1,0.5)
#beta<-c(1,0.75,0.5,0.25,0.1,-0.1,-0.25,-0.5,-0.75,-1)
Sxx<-10 #Look at effects
optima.bar<-0
optima.sd<-2
beta.bar<-0.5
beta.sd<-1
rho<- -0.8

#Create new phylogeny
set.seed(10)
phy<-ape::rtree(n=N)
phy<-ape::multi2di(phy)
l.tree<-max(ape::branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree

#Set regimes - manually - 2 regimes
#Locate nodes
plot(phy,no.margin=TRUE,edge.width=2,cex=0.7)
ape::nodelabels(frame="none",adj=c(1.1,-0.4))
#ape::tiplabels()

#shifts<-c(338,276,213,214,244,258,279,296,320) #10 regimes, 20 species per regime
shifts<-c(189,164,169,175,111,130,127,121,102) #10 regimes, 10 species per regime

mlm.data<-sim.reg.direct.mlm.ve.data(phy,N,Z_direct,hl,vy,Sxx,optima.bar,optima.sd,beta.bar,beta.sd,rho,shifts)

trdata<-mlm.data[[1]]
vary.effects<-mlm.data[[2]]
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.25)
vy.prior<-20
optima.prior<-c(0,2.5) #Informed by linear model
beta.prior<-c(0,1.25)
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_ve.stan")
fit.reg.direct.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 1,cores=1,iter =2000)

post<-rstan::extract(fit.reg.direct.mlm.ve)
print(fit.reg.direct.mlm.ve,pars = c("hl","vy","optima","optima_bar","beta","beta_bar","sigma","Rho[1,2]"))

############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-reg.direct.mlm.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,vary.effects)
sigma.prior.post.plot<-sigma.prior.post.plot.code(sigma.prior,post)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect Model- Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Figures/Fig_Direct_mlm_ve.pdf", plot = fig, width=7, height=7 )
########################################################################################################
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects - Non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_ve_nc.stan")
fit.reg.direct.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 1,cores=1,iter =2000)

post<-rstan::extract(fit.reg.direct.mlm.ve.nc)
print(fit.reg.direct.mlm.ve.nc,pars = c("hl","vy","optima","optima_bar","beta","beta_bar","sigma","Z","L_Rho[1,2]"))

############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-reg.direct.mlm.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,vary.effects)
sigma.prior.post.plot<-sigma.prior.post.plot.code(sigma.prior,post)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Direct Effect Model- Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Figures/Fig_Direct_mlm_ve_nc.pdf", plot = fig, width=7, height=7 )

############################################################################################################
#Non MLM version of Same Model - Multi-Optima Direct Effect Model - Varying Effects - Non-centered
############################################################################################################
dat<-blouch.reg.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_direct_ve.stan")
fit.reg.direct.ve<- rstan::sampling(object = stan_model,data = dat,chains = 1,cores=1,iter =2000)

post<-rstan::extract(fit.reg.direct.ve)
print(fit.reg.direct.ve,pars = c("hl","vy","optima","beta"))

############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.direct.plot.code(hl.prior,vy.prior,post)
regression.plot<-reg.direct.mlm.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,vary.effects)
sigma.prior.post.plot<-sigma.prior.post.plot.code(sigma.prior,post)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multi-Optima Direct Effect Model- Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Figures/Fig_Direct_ve.pdf", plot = fig, width=7, height=7 )

############################################################################################################
#Compare models using PSIS
loo_mlm_ve <- loo::loo(fit.reg.direct.mlm.ve, save_psis = TRUE)
print(loo_mlm_ve)
plot(loo_mlm_ve) #4X6
plot(loo_mlm_ve,label_points=TRUE) #Label outliers

loo_mlm_ve_nc <- loo::loo(fit.reg.direct.mlm.ve.nc, save_psis = TRUE)
print(loo_mlm_ve_nc)
plot(loo_mlm_ve_nc) #4X6
plot(loo_mlm_ve_nc,label_points=TRUE) #Label outliers

loo_ve <- loo::loo(fit.reg.direct.ve, save_psis = TRUE)
print(loo_ve)
plot(loo_ve) #4X6
plot(loo_ve,label_points=TRUE) #Label outliers

loo::loo_compare(loo_mlm_ve, loo_mlm_ve_nc, loo_ve)

1.3 + c(-1,1)*0.5*1.86 #95% comparability interval of the difference, z-score of 1.86










############################################################################################################
########################################################################################################
#Multilevel Multi-Optima Adaptive Model - Varying Effects - non-centered - 200 species, 10 regimes
########################################################################################################
set.seed(10)
N<-200
Z_adapt<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
Sxx<-10 #Look at effects
optima.bar<-0
optima.sd<-2
beta.bar<-1
beta.sd<-1
rho<- -0.8
shifts<-c(338,276,213,214,244,258,279,296,320) #10 regimes, 20 species per regime
#shifts<-c(189,164,169,175,111,130,127,121,102) #10 regimes

set.seed(10)
phy<-ape::rtree(n=N)
phy<-ape::multi2di(phy)
l.tree<-max(ape::branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree

#Set regimes - manually - 2 regimes
#Locate nodes
plot(phy,no.margin=TRUE,edge.width=2,cex=0.7)
ape::nodelabels(frame="none",adj=c(1.1,-0.4))
#ape::tiplabels()


mlm.data<-sim.reg.adapt.mlm.ve.data(phy,N,Z_adapt,hl,vy,Sxx,optima.bar,optima.sd,beta.bar,beta.sd,rho,shifts)

trdata<-mlm.data[[1]]
vary.effects<-mlm.data[[2]]
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,3) #Informed by linear model
beta.prior<-c(0,1)
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adapt=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects - 
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve_nc.stan")
fit.reg.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 1,cores=1,iter =2000)

post<-rstan::extract(fit.reg.adapt.mlm.ve.nc)
print(fit.reg.adapt.mlm.ve.nc,pars = c("hl","vy","optima","optima_bar","beta","beta_bar","sigma","Rho"))

############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.mlm.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,vary.effects)
sigma.prior.post.plot<-sigma.prior.post.plot.code(sigma.prior,post)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Figures/Fig_Adapt_ve_nc.pdf", plot = fig, width=7, height=7 )

########################################################################################################
#Multi-Optima Direct Effect Model - Varying Effects - 
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adapt=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################

setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_ve.stan")
fit.reg.adapt.ve<- rstan::sampling(object = stan_model,data = dat,chains = 1,cores=1,iter =2000)

post<-rstan::extract(fit.reg.adapt.ve)
print(fit.reg.adapt.ve,pars = c("hl","vy","optima","beta"))

############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.mlm.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,vary.effects)
sigma.prior.post.plot<-sigma.prior.post.plot.code(sigma.prior,post)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Figures/Fig_Adapt_ve_nc.pdf", plot = fig, width=7, height=7 )


############################################################################################################
#Compare models using PSIS
loo_mlm_ve_nc <- loo::loo(fit.reg.adapt.mlm.ve.nc, save_psis = TRUE)
print(loo_mlm_ve_nc)
plot(loo_mlm_ve_nc) #4X6
plot(loo_mlm_ve_nc,label_points=TRUE) #Label outliers

loo_ve <- loo::loo(fit.reg.adapt.ve, save_psis = TRUE)
print(loo_ve)
plot(loo_ve) #4X6
plot(loo_ve,label_points=TRUE) #Label outliers

loo::loo_compare(loo_mlm_ve_nc, loo_ve)


1.2 + c(-1,1)*2*1.86 #95% comparability interval of the difference, z-score of 1.86








########################################################################################################
########################################################################################################
#Multilevel Multi-Optima Adaptive Model - Varying Effects - non-centered - 100 species, 5 regimes
########################################################################################################
N<-100
Z_adaptive<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
Sxx<-10 #Look at effects
optima.bar<-0
optima.sigma<-2
beta.bar<-0
beta.sigma<-1
rho<- -0.25

set.seed(10)
phy<-ape::rtree(n=N)
phy<-ape::multi2di(phy)
l.tree<-max(ape::branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree
plot(phy)
ape::nodelabels()

#shifts<-c(338,276,213,214,244,258,279,296,320) #10 regimes, 20 species per regime
#shifts<-c(189,164,169,175,111,130,127,121,102) #10 regimes, 10 species per regime
shifts<-c(102,111,186,157) #5 Regimes 

#mlm.data<-sim.reg.adapt.mlm.ve.data(phy,N.regimes,avg.N.sp.regimes,Z_adapt,hl,vy,Sxx,optima.bar,optima.sd,beta.bar,beta.sd,rho,shifts)
mlm.data<-sim.reg.adapt.mlm.ve.data(phy,N,Z_adaptive,hl,vy,Sxx,optima.bar,optima.sigma,beta.bar,beta.sigma,rho,shifts)

trdata<-mlm.data[[1]]
vary.effects<-mlm.data[[2]]
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
beta.prior<-c(0,1)
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adapt=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptive Model - Varying Effects - Non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve_nc.stan")
fit.reg.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 1,cores=1,iter =2000)

post<-rstan::extract(fit.reg.adapt.mlm.ve.nc)
print(fit.reg.adapt.mlm.ve.nc,pars = c("hl","vy","optima","optima_bar","beta","beta_bar","sigma","Rho[1,2]"))

############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.mlm.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,vary.effects)
sigma.prior.post.plot<-sigma.prior.post.plot.code(sigma.prior,post)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Figures/Fig_Adapt_mlm_ve_nc.pdf", plot = fig, width=7, height=7 )

############################################################################################################
#Non MLM version of Same Model - Multi-Optima Adaptation Model - Varying Effects
############################################################################################################
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adapt=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Adaptive Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_ve.stan")
fit.reg.adapt.ve<- rstan::sampling(object = stan_model,data = dat,chains = 1,cores=1,iter =2000)

post<-rstan::extract(fit.reg.adapt.ve)
print(fit.reg.adapt.ve,pars = c("hl","vy","optima","beta"))

############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.mlm.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,vary.effects)
sigma.prior.post.plot<-sigma.prior.post.plot.code(sigma.prior,post)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Figures/Fig_Adapt_ve.pdf", plot = fig, width=7, height=7 )

############################################################################################################
#Compare models using PSIS
loo_mlm_ve_nc <- loo::loo(fit.reg.adapt.mlm.ve.nc, save_psis = TRUE)
print(loo_mlm_ve_nc)
plot(loo_mlm_ve_nc) #4X6
plot(loo_mlm_ve_nc,label_points=TRUE) #Label outliers

loo_ve <- loo::loo(fit.reg.adapt.ve, save_psis = TRUE)
print(loo_ve)
plot(loo_ve) #4X6
plot(loo_ve,label_points=TRUE) #Label outliers

loo::loo_compare(loo_mlm_ve_nc, loo_ve)


0.7 + c(-1,1)*0.5*1.86 #95% comparability interval of the difference, z-score of 1.86



#START
############################################################################################################
########################################################################################################
#Multilevel Multi-Optima Adaptive Model - Varying Effects - non-centered - 10 regimes, 100 tips
########################################################################################################
N<-100
N.regimes<-10
Z_adaptive<-1
hl<-0.1
vy<-0.01 #0.25,0.5 - testing options
Sxx<-10 #Look at effects
optima.bar<-0
optima.sigma<-1
beta.bar<-0.5
beta.sigma<-0.5
rho<- -0.8

set.seed(10)
phy<-ape::rtree(n=N)
phy<-ape::multi2di(phy)
l.tree<-max(ape::branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree
plot(phy)
ape::nodelabels()

#shifts<-c(338,276,213,214,244,258,279,296,320) #10 regimes, 20 species per regime
shifts<-c(189,164,169,175,111,130,127,121,102) #10 regimes, 10 species per regime
#shifts<-c(94,61,82,65) #5 Regimes, 10 species per regime

#mlm.data<-sim.reg.adapt.mlm.ve.data(phy,N.regimes,avg.N.sp.regimes,Z_adapt,hl,vy,Sxx,optima.bar,optima.sd,beta.bar,beta.sd,rho,shifts)
mlm.data<-sim.reg.adapt.mlm.ve.data(phy,N,Z_adaptive,hl,vy,Sxx,optima.bar,optima.sigma,beta.bar,beta.sigma,rho,shifts)

trdata<-mlm.data[[1]]
vary.effects<-mlm.data[[2]]
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
beta.prior<-c(0.5,0.5)
sigma.prior<-c(0,1)
############################################################################################################
#Explore Priors
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
#sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adapt=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptive Model - Varying Effects - Non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve_nc.stan")
fit.reg.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 1,cores=1,iter =2000)

post<-rstan::extract(fit.reg.adapt.mlm.ve.nc)
print(fit.reg.adapt.mlm.ve.nc,pars = c("hl","vy","optima","optima_bar","beta","beta_bar","sigma","Rho[1,2]"))
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.mlm.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,vary.effects)
#sigma.prior.post.plot<-sigma.prior.post.plot.code(sigma.prior,post)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Figures/Fig_Adapt_mlm_ve_nc.pdf", plot = fig, width=7, height=7 )
########################################################################################################
#Multilevel Multi-Optima Adaptive Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve_nc.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve_nc.stan")
fit.reg.adapt.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 1,cores=1,iter =2000)

post<-rstan::extract(fit.reg.adapt.mlm.ve)
print(fit.reg.adapt.mlm.ve,pars = c("hl","vy","optima","optima_bar","beta","beta_bar","sigma","Rho[1,2]"))
############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.mlm.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,vary.effects)
#sigma.prior.post.plot<-sigma.prior.post.plot.code(sigma.prior,post)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Figures/Fig_Adapt_mlm_ve.pdf", plot = fig, width=7, height=7 )

############################################################################################################
#Non MLM version of Same Model - Multi-Optima Adaptation Model - Varying Effects
############################################################################################################
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adapt=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Adaptive Model - Varying Effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan")
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_ve.stan")
stan_model <- rstan::stan_model("blouchOU_reg_adapt_ve.stan")
fit.reg.adapt.ve<- rstan::sampling(object = stan_model,data = dat,chains = 1,cores=1,iter =2000)

post<-rstan::extract(fit.reg.adapt.ve)
print(fit.reg.adapt.ve,pars = c("hl","vy","optima","beta"))

############################################################################################################
hl.plot<-hl.prior.post.plot.code(hl.prior,post,hl)
vy.plot<-vy.prior.post.plot.code(vy.prior,post,vy)
covariance.plot<-covariance.prior.post.adapt.plot.code(hl.prior,vy.prior,beta.prior,post)
regression.plot<-reg.adapt.mlm.ve.prior.post.plot.code(trdata,optima.prior,beta.prior,post,vary.effects)
#sigma.prior.post.plot<-sigma.prior.post.plot.code(sigma.prior,post)
############################################################################################################
fig<-ggpubr::ggarrange(hl.plot, vy.plot,covariance.plot,regression.plot,ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")
fig<-ggpubr::annotate_figure(fig,top=paste("Multilevel Multi-Optima Adaptive Model - Varying Effects\n - Non-centered, hl=",hl,sep=""))
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Figures/Fig_Adapt_ve.pdf", plot = fig, width=7, height=7 )############################################################################################################
############################################################################################################
#Compare models using PSIS
loo_mlm_ve_nc <- loo::loo(fit.reg.adapt.mlm.ve.nc, save_psis = TRUE)
print(loo_mlm_ve_nc)
plot(loo_mlm_ve_nc) #4X6
plot(loo_mlm_ve_nc,label_points=TRUE) #Label outliers

loo_mlm_ve <- loo::loo(fit.reg.adapt.mlm.ve, save_psis = TRUE)
print(loo_mlm_ve)
plot(loo_mlm_ve) #4X6
plot(loo_mlm_ve_nc,label_points=TRUE) #Label outliers


loo_ve <- loo::loo(fit.reg.adapt.ve, save_psis = TRUE)
print(loo_ve)
plot(loo_ve) #4X6
plot(loo_ve,label_points=TRUE) #Label outliers

loo::loo_compare(loo_mlm_ve_nc, loo_mlm_ve, loo_ve)
############################################################################################################

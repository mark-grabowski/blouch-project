#Blouch Model Checking Code - 2024
#SBR2 04/09/24
rm(list=ls())
########################################################################################################
#Prior and posterior predictive checks for all models
#Complete code to make simulated data, setup, run model checking for SI
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
library(rstan)
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
############################################################################################################
########################################################################################################
#Direct Effect Model
########################################################################################################
set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-2
beta<-0.25
trdata<-direct.effect.sim(tree.10K,N,Z,hl,vy,Sxx,optima,beta)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(2,0.25) #Informed by linear model
beta.prior<-c(0.25,0.25) #Informed by linear model
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
##reg.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch prep code - direct effect model - blouch.direct.prep()
dat<-blouch.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Prior Predictive Check - Direct effect model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_direct_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_direct_priorpc.stan")
#fit.direct.prior<- rstan::sampling(object = blouch:::stanmodels$blouchOU_direct_priorpc,data = dat,chains = 2,cores=2,iter =2000,algorithm=c("Fixed_param"))
fit.direct.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000,algorithm=c("Fixed_param"))
post<-rstan::extract(fit.direct.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Direct Effect Model - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_1.pdf",
                plot = plots, width=10, height=2,)

plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Posterior Predictive Check - Direct effect model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_direct_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_direct_postpc.stan")
fit.direct.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.direct.postpc,pars = c("hl","vy","optima","beta"))
post<-rstan::extract(fit.direct.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Direct Effect Model - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_2.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
############################################################################################################
########################################################################################################
#Adaptation Models
########################################################################################################
set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-2
beta<-0.25
trdata<-adaptation.sim(tree.10K,N,Z,hl,vy,Sxx,optima,beta)
############################################################################################################
#Set Priors
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(2,0.25) #Informed by linear model
beta.prior<-c(0.25,0.25) #Informed by linear model
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.adapt.plot.code(hl.prior,vy.prior,beta.prior)
#reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch prep code
dat<-blouch.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Prior Predictive Check - Adaptive model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_adapt_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_adapt_priorpc.stan")

fit.adaptive.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000,algorithm=c("Fixed_param"))
post<-rstan::extract(fit.adaptive.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Adaptive Model - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_3.pdf",
                plot = plots, width=10, height=2)
plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Posterior Predictive Check - Adaptive model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_adapt_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_adapt_postpc.stan")

fit.adaptive.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2, cores=2,iter =2000)
print(fit.adaptive.postpc,pars = c("hl","vy","optima","beta","beta_e"))
post<-rstan::extract(fit.adaptive.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Adaptive Model - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_4.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
############################################################################################################
########################################################################################################
#Direct Effect + Adaptation Model
########################################################################################################
set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
Z_adaptive<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-2
beta<-c(0.35,0.25)
trdata<-direct.effect.adaptation.sim(tree.10K,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta)
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
#reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch prep code
dat<-blouch.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,hl.prior,vy.prior,optima.prior,beta.prior)
############################################################################################################
#Direct Effect and Adaptation Model - Prior Predictive Checks
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_direct_adapt_priorpc.stan")
stan_model <- rstan::stan_model("blouchOU_direct_adapt_priorpc.stan")
fit.mixed.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.mixed.priorpc)
num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Direct Effect + Adaptation Model - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_5.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Direct Effect and Adaptation Model - Posterior Predictive Checks
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_direct_adapt_postpc.stan")
stan_model <- rstan::stan_model("blouchOU_direct_adapt_postpc.stan")
fit.mixed.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

post<-rstan::extract(fit.mixed.postpc)
num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Direct Effect + Adaptation Model - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_6.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
############################################################################################################
########################################################################################################
#Multi-Optima Model - Varying Intercepts
########################################################################################################
set.seed(10) #Set seed to get same random species each time
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
trdata<-reg.sim(tree.10K,N,hl,vy,Sxx,optima,shifts)
############################################################################################################
hl.prior.plot.code(hl.prior)
vy.prior.plot.code(vy.prior)
covariance.prior.direct.plot.code(hl.prior,vy.prior)
optima.prior.plot.code(trdata,optima.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.prep(trdata,"Y_with_error","Y_error","regimes",hl.prior,vy.prior,optima.prior)
############################################################################################################
#Run Multi-Optima Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_priorpc.stan")
fit.reg.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.priorpc)
num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Model - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_7.pdf",
                plot = plots, width=10, height=2)

############################################################################################################
#Run Multi-Optima Model - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_postpc.stan")
fit.reg.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.postpc,pars = c("hl","vy","optima"))

post<-rstan::extract(fit.reg.postpc)
num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Model - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_8.pdf",
                plot = plots, width=10, height=2)
############################################################################################################
#Multilevel Multi-Optima Model - Varying Intercepts
############################################################################################################
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
sigma.prior<-c(0,1)
trdata<-reg.sim(tree.10K,N,hl,vy,Sxx,optima,shifts)
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
#Run Multi-Optima Model
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_mlm_vi_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_mlm_vi_priorpc.stan")
fit.reg.mlm.vi.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.mlm.vi.priorpc)
num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Model - Varying Intercepts - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_9.pdf",
                plot = plots, width=10, height=2)
############################################################################################################
#Run Multi-Optima Model - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_mlm_vi_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_mlm_vi_postpc.stan")
fit.reg.mlm.vi.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.mlm.vi.postpc,pars = c("hl","vy","optima","optima_bar","sigma"))

post<-rstan::extract(fit.reg.mlm.vi.postpc)
num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Model - Varying Intercepts - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_10.pdf",
                plot = plots, width=10, height=2)

############################################################################################################
#Multilevel Multi-Optima Model- Varying Intercepts - Non-Centered
############################################################################################################
hl.prior<-c(log(0.25),0.75)
vy.prior<-20
optima.prior<-c(0,1) #Informed by linear model
sigma.prior<-c(0,1)
trdata<-reg.sim(tree.10K,N,hl,vy,Sxx,optima,shifts)
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
#Multi-Optima Model - Non-centered - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_mlm_vi_nc_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_mlm_vi_nc_priorpc.stan")
fit.reg.mlm.vi.nc.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.mlm.vi.nc.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Model - Varying Intercepts -Non-centered - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_11.pdf",
                plot = plots, width=10, height=2)
############################################################################################################
#Multi-Optima Model - Non-centered - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_mlm_vi_nc_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_mlm_vi_nc_postpc.stan")
fit.reg.mlm.vi.nc.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.mlm.vi.nc.postpc,pars = c("hl","vy","optima","optima_bar","sigma"))

post<-rstan::extract(fit.reg.mlm.vi.nc.postpc)
num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Model - Varying Intercepts - Non-centered - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_12.pdf",
                plot = plots, width=10, height=2)

############################################################################################################
########################################################################################################
#Multi-Optima Direct Effect Model - Varying Intercepts
########################################################################################################
set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.25)
shifts<-c(84)
trdata<-reg.direct.sim(tree.10K,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
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
#reg.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect Model - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_priorpc.stan")
fit.reg.direct.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.direct.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Direct Effect Model - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_13.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Multi-Optima Direct Effect Model - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_postpc.stan")

fit.reg.direct.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.direct.postpc,pars = c("hl","vy","optima","beta"))
post<-rstan::extract(fit.reg.direct.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Direct Effect Model - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_14.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)
############################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts
########################################################################################################
set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.25)
shifts<-c(84)
trdata<-reg.direct.sim(tree.10K,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
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
#reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts- Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_vi_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_vi_priorpc.stan")
fit.reg.direct.mlm.vi.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.direct.mlm.vi.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect Model - Varying Intercepts - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_15.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_vi_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_vi_postpc.stan")

fit.reg.direct.mlm.vi.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.direct.mlm.vi.postpc,pars = c("hl","vy","optima","optima_bar","beta","sigma"))
post<-rstan::extract(fit.reg.direct.mlm.vi.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect Model - Varying Intercepts - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_16.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)

############################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts - Non-centered
########################################################################################################
set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.25)
shifts<-c(84)
trdata<-reg.direct.sim(tree.10K,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
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
#reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts- Non-centered - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_vi_nc_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_vi_nc_priorpc.stan")
fit.reg.direct.mlm.vi.nc.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.direct.mlm.vi.nc.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect Model - Varying Intercepts - Non-centered - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_17.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Intercepts - Non-centered - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_vi_nc_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_vi_nc_postpc.stan")

fit.reg.direct.mlm.vi.nc.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.direct.mlm.vi.nc.postpc,pars = c("hl","vy","optima","optima_bar","beta","sigma"))
post<-rstan::extract(fit.reg.direct.mlm.vi.nc.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect Model - Varying Intercepts - Non-centered - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_18.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)
########################################################################################################
########################################################################################################
#Multi-Optima Adaptation Model with Varying Intercepts
########################################################################################################
set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_adapt<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.25)
shifts<-c(84)
trdata<-reg.adapt.sim(tree.10K,N,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
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
#reg.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
#Blouch Prep Code
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Adaptation Model - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_adapt_priorpc.stan")
fit.reg.adapt.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.adapt.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Adaptation Model - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_19.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Multi-Optima Adaptation Model - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_adapt_postpc.stan")

fit.reg.adapt.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.adapt.postpc,pars = c("hl","vy","optima","beta"))
post<-rstan::extract(fit.reg.adapt.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Adaptation Model - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_20.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)

########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Intercepts 
########################################################################################################
set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_adapt<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.25)
shifts<-c(84)
trdata<-reg.adapt.sim(tree.10K,N,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
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
#reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)

############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Intercepts -Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_vi_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_vi_priorpc.stan")
fit.reg.adapt.mlm.vi.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.adapt.mlm.vi.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Adaptation Model - Varying Intercepts -Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_21.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Intercepts - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_vi_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_vi_postpc.stan")

fit.reg.adapt.mlm.vi.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.adapt.mlm.vi.postpc,pars = c("hl","vy","optima","optima_bar","beta","sigma"))
post<-rstan::extract(fit.reg.adapt.mlm.vi.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Adaptation Model - Varying Intercepts - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_22.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)

########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Intercepts - Non-centered - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_vi_nc_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_vi_nc_priorpc.stan")
fit.reg.adapt.mlm.vi.nc.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.adapt.mlm.vi.nc.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Adaptation Model - Varying Intercepts - Non-centered - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_23.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Multilevel Multi-Optima Adaptation Model - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_vi_nc_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_vi_nc_postpc.stan")

fit.reg.adapt.mlm.vi.nc.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.adapt.mlm.vi.nc.postpc,pars = c("hl","vy","optima","optima_bar","beta","sigma"))
post<-rstan::extract(fit.reg.adapt.mlm.vi.nc.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Adaptation Model - Varying Intercepts - Non-centered - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_24.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)

############################################################################################################
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model with Varying Intercepts
########################################################################################################
set.seed(10) #Set seed to get same random species each time
N<-50 #Number of species
Z_direct<-1
Z_adapt<-1
hl<-0.1
vy<-0.01
Sxx<-10
optima<-c(2,1)
beta<-c(0.35,0.25) #Two Optima/Two Slopes
shifts<-c(84)
trdata<-reg.direct.adapt.sim(tree.10K,N,Z_direct,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
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
#reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
dat<-blouch.reg.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_priorpc.stan")
fit.reg.direct.adapt.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.direct.adapt.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Direct Effect + Adaptation Model - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_25.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Multi-Optima Direct Effect Adaptation Model - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_postpc.stan")
fit.reg.direct.adapt.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.direct.adapt.postpc,pars = c("hl","vy","optima","beta"))

post<-rstan::extract(fit.reg.direct.adapt.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Direct Effect + Adaptation Model - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_26.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)
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
trdata<-reg.direct.adapt.sim(tree.10K,N,Z_direct,Z_adapt,hl,vy,Sxx,optima,beta,shifts)
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
#reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_vi_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_vi_priorpc.stan")
fit.reg.direct.adapt.mlm.vi.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.direct.adapt.mlm.vi.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_27.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_vi_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_vi_postpc.stan")
fit.reg.direct.adapt.mlm.vi.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.direct.adapt.mlm.vi.postpc,pars = c("hl","vy","optima","optima_bar","beta","sigma"))

post<-rstan::extract(fit.reg.direct.adapt.mlm.vi.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_28.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)

########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts - Non-Centered - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_vi_nc_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_vi_nc_priorpc.stan")
fit.reg.direct.adapt.mlm.vi.nc.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.direct.adapt.mlm.vi.nc.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts - Non-centered - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_29.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts - Non-Centered - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_vi_nc_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_vi_nc_postpc.stan")
fit.reg.direct.adapt.mlm.vi.nc.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.direct.adapt.mlm.vi.nc.postpc,pars = c("hl","vy","optima","optima_bar","beta","sigma"))

post<-rstan::extract(fit.reg.direct.adapt.mlm.vi.nc.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Intercepts - Non-centered - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_30.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)

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
trdata<-reg.direct.ve.sim(tree.10K,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
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
#reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
############################################################################################################
dat<-blouch.reg.direct.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect Model - Varying Effects - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_ve_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_ve_priorpc.stan")
fit.reg.direct.ve.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.direct.ve.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Direct Effect Model - Varying Efects - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_31.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Multi-Optima Direct Effect Model - Varying Effects - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_ve_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_ve_postpc.stan")
fit.reg.direct.ve.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.direct.ve.postpc,pars = c("hl","vy","optima","beta"))

post<-rstan::extract(fit.reg.direct.ve.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Direct Effect Model - Varying Efects - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_32.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects
N<-50
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10 #Look at effects
shifts<-c(84)
trdata<-reg.direct.ve.sim(tree.10K,N,Z_direct,hl,vy,Sxx,optima,beta,shifts)
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
#reg.direct.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.direct.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_direct=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_ve_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_ve_priorpc.stan")
fit.reg.direct.mlm.ve.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.direct.mlm.ve.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect Model - Varying Efects - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_33.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_ve_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_ve_postpc.stan")
fit.reg.direct.mlm.ve.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.direct.mlm.ve.postpc,pars = c("hl","vy","optima","optima_bar","beta","beta_bar","sigma","Rho"))

post<-rstan::extract(fit.reg.direct.mlm.ve.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect Model - Varying Efects - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_34.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects - Non-centered - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_ve_nc_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_ve_nc_priorpc.stan")
fit.reg.direct.mlm.ve.nc.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.direct.mlm.ve.nc.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect Model - Varying Efects - Non-centered - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_35.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Multilevel Multi-Optima Direct Effect Model - Varying Effects - Non-centered - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_mlm_ve_nc_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_mlm_ve_nc_postpc.stan")
fit.reg.direct.mlm.ve.nc.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.direct.mlm.ve.nc.postpc,pars = c("hl","vy","optima","optima_bar","beta","beta_bar","sigma","L_Rho"))

post<-rstan::extract(fit.reg.direct.mlm.ve.nc.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect Model - Varying Efects - Non-centered Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_36.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)

############################################################################################################
########################################################################################################
## Multi-Optima Adaptation Effect Model - Varying Effects
########################################################################################################
N<-50
Z_adaptive<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10
shifts<-c(84)
trdata<-reg.adapt.ve.sim(tree.10K,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
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
#reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
########################################################################################################
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Adaptation Model - Varying Effects - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_ve_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_adapt_ve_priorpc.stan")
fit.reg.adapt.ve.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.adapt.ve.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Adaptation Model - Varying Efects - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_37.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Multi-Optima Adaptation Model - Varying Effects - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_ve_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_adapt_ve_postpc.stan")
fit.reg.adapt.ve.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.adapt.ve.postpc,pars = c("hl","vy","optima","beta"))

post<-rstan::extract(fit.reg.adapt.ve.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Adaptation Model - Varying Efects - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_38.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)
############################################################################################################
############################################################################################################
#Multilevel Multi-Optima Adaptation Model - Varying Effects - WORKING
############################################################################################################
N<-50
Z_adaptive<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1)
beta<-c(0.25,0.15) #Two Optima/Two Slopes
Sxx<-10
shifts<-c(84)
trdata<-reg.adapt.ve.sim(tree.10K,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
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
#reg.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
############################################################################################################
dat<-blouch.reg.adapt.mlm.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error",Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)

########################################################################################################
#Multi-Optima Adaptation Model - Varying Effects - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve_priorpc.stan")
fit.reg.adapt.mlm.ve.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.adapt.mlm.ve.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Adaptation Model - Varying Efects - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_39.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Multi-Optima Adaptation Model - Varying Effects - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve_postpc.stan")
fit.reg.adapt.mlm.ve.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.adapt.mlm.ve.postpc,pars = c("hl","vy","optima","beta"))

post<-rstan::extract(fit.reg.adapt.mlm.ve.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Adaptation Model - Varying Efects - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_40.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)


########################################################################################################
#Multi-Optima Adaptation Model - Varying Effects - Non-centered - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve_nc_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve_nc_priorpc.stan")
fit.reg.adapt.mlm.ve.nc.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.adapt.mlm.ve.nc.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Adaptation Model - Varying Efects - Non-centered - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_41.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Multi-Optima Adaptation Model - Varying Effects - Non-centered - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_adapt_mlm_ve_nc_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_adapt_mlm_ve_nc_postpc.stan")
fit.reg.adapt.mlm.ve.nc.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.adapt.mlm.ve.nc.postpc,pars = c("hl","vy","optima","beta"))

post<-rstan::extract(fit.reg.adapt.mlm.ve.nc.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Adaptation Model - Varying Efects - Non-centered - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_42.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)
############################################################################################################
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects
########################################################################################################
N<-50
Z_direct<-1
Z_adaptive<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1.5)
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Direct effect - column 1, Adaptive - column 2, slopes for regimes on rows
Sxx<-10
shifts<-c(84)
trdata<-reg.direct.adapt.ve.sim(tree.10K,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
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
#reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
########################################################################################################
dat<-blouch.reg.direct.adapt.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_ve_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_ve_priorpc.stan")
fit.reg.direct.adapt.ve.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.direct.adapt.ve.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Direct Effect + Adaptation Model - Varying Efects - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_43.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Multi-Optima Direct Effect + Adaptation Model - Varying Effects - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_ve_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_ve_postpc.stan")
fit.reg.adapt.ve.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.adapt.ve.postpc,pars = c("hl","vy","optima","beta"))

post<-rstan::extract(fit.reg.adapt.ve.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multi-Optima Direct Effect + Adaptation Model - Varying Efects - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_44.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Effects
########################################################################################################
N<-50
Z_direct<-1
Z_adaptive<-1
hl<-0.1 #0.1, 0.25, 0.75 - testing options
vy<-0.01 #0.25,0.5 - testing options
optima<-c(2,1.5)
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Direct effect - column 1, Adaptive - column 2, slopes for regimes on rows
Sxx<-10
shifts<-c(84)
trdata<-reg.direct.adapt.ve.sim(tree.10K,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts)
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
#reg.direct.adapt.prior.plot.code(trdata,optima.prior,beta.prior)
sigma.prior.plot.code(sigma.prior)
########################################################################################################
dat<-blouch.reg.direct.adapt.mlm.prep(trdata,"Y_with_error","Y_error",c("Xd","Xa"),c("Xd_error","Xa_error"),Z_direct=1,Z_adaptive=1,"regimes",hl.prior,vy.prior,optima.prior,beta.prior,sigma.prior)
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Effects - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_ve_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_ve_priorpc.stan")
fit.reg.direct.adapt.ve.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.direct.adapt.ve.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Efects - Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_45.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Effects - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_ve_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_ve_postpc.stan")
fit.reg.adapt.ve.mlm.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.adapt.ve.mlm.postpc,pars = c("hl","vy","optima","optima_bar","beta","beta_bar","sigma","Rho"))

post<-rstan::extract(fit.reg.adapt.ve.mlm.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Efects - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_46.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)

########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Effects - Non-centered - Prior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_ve_nc_priorpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_ve_nc_priorpc.stan")
fit.reg.direct.adapt.mlm.ve.nc.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2, iter =2000,algorithm=c("Fixed_param"))

post<-rstan::extract(fit.reg.direct.adapt.mlm.ve.nc.priorpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Efects - Non-centered -Prior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_47.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Effects - Non-centered - Posterior Predictive Check
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/") #Macbook Pro
rstan::stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/R Packages - SBR2 - cloned/testing_v3/blouch/inst/stan/blouchOU_reg_direct_adapt_mlm_ve_nc_postpc.stan") #Macbook Pro
stan_model <- rstan::stan_model("blouchOU_reg_direct_adapt_mlm_ve_nc_postpc.stan")
fit.reg.adapt.mlm.ve.nc.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.adapt.mlm.ve.nc.postpc,pars = c("hl","vy","optima","optima_bar","beta","beta_bar","sigma","L_Rho"))

post<-rstan::extract(fit.reg.adapt.mlm.ve.nc.postpc)

num.plots<-5
plots<-ysim.ppc.plot.code(dat,post,1:num.plots)
plots<-gridExtra::grid.arrange(grobs=plots, ncol=5, nrow=1,common.legend=TRUE,legend="top")
tgrob <- ggpubr::text_grob("Multilevel Multi-Optima Direct Effect + Adaptation Model - Varying Efects - Non-centered - Posterior PC",size = 10)
plots<-ggpubr::annotate_figure(plots,top=tgrob)
ggplot2::ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch Project - SBR2/Blouch ms/R2/Latex Documents/Figures/Model Checking/Fig_PPC_48.pdf",
                plot = plots, width=10, height=2)

plot(post$Y_sim_obs[4,],dat$Y_obs)

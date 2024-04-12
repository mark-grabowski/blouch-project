#Blouch Model Checking Code - 2023



library(rstan)
library(rethinking)
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
#Prior and posterior predictive checks for all models

########################################################################################################
#Milestone 1 - Prior Predictive Code
#Using direct effect model with ME
#Also works fo multivariariate Xs
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_direct_priorpcheck.stan") #Macbook Pro
stan_model <- stan_model("blouchOU_direct_priorpcheck.stan")
fit.npi.direct.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000, cores=2,algorithm=c("Fixed_param"))
#print(fit.npi.direct,pars = c("hl","vy","alpha","beta"))
post<-extract(fit.npi.direct.priorpc)
plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Milestone 1 - Posterior Predictive Code
#Using direct effect model with ME
#Also works fo multivariariate Xs
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_direct_postpcheck.stan") #Macbook Pro
stan_model <- stan_model("blouchOU_direct_postpcheck.stan")
fit.npi.direct.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000, cores=2)
#print(fit.npi.direct.postpc,pars = c("hl","vy","alpha","beta"))
post<-extract(fit.npi.direct.postpc)
plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Milestone 2 - Prior Predictive Code
#Using adaptive model with ME 
#Also works for multivariariate Xs
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_adapt_priorpc.stan") #Macbook Pro
stan_model <- stan_model("blouchOU_adapt_priorpc.stan")

fit.npi.adaptive.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000, cores=2,algorithm=c("Fixed_param"))
#print(fit.npi.adaptive,pars = c("hl","vy","alpha","beta","beta_e"))
post<-extract(fit.npi.adaptive.priorpc)
plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Milestone 2 - Posterior Predictive Code
#Using adaptive model with ME 
#Also works for multivariariate Xs
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_adapt_postpc.stan") #Macbook Pro
stan_model <- stan_model("blouchOU_adapt_postpc.stan")

fit.npi.adaptive.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000, cores=2)
post<-extract(fit.npi.adaptive.postpc)
plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Milestone 3 - Prior Predictive Code
#Combination of direct effect and adaptive predictors with ME
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_direct_adapt_priorpc.stan")
stan_model <- stan_model("blouchOU_direct_adapt_priorpc.stan")
fit.npi.mixed.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2,algorithm=c("Fixed_param"))

#print(fit.npi.mixed,pars = c("hl","vy","alpha","beta","beta_e"))
post<-extract(fit.npi.mixed.priorpc)
plot(post$Y_sim_obs[2,],dat$Y_obs)
########################################################################################################
#Milestone 3 - Posterior Predictive Code
#Combination of direct effect and adaptive predictors with ME
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_direct_adapt_postpc.stan")
stan_model <- stan_model("blouchOU_direct_adapt_postpc.stan")
fit.npi.mixed.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

#print(fit.npi.mixed,pars = c("hl","vy","alpha","beta","beta_e"))
post<-extract(fit.npi.mixed.postpc)
plot(post$Y_sim_obs[2,],dat$Y_obs)

########################################################################################################
#Milestone 4 - Prior Predictive Check
#Regime model - regimes - painted or SIMMAP
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_reg_priorpcheck.stan")
stan_model <- stan_model("blouchOU_reg_priorpcheck.stan")

fit.npi.regimes.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2, algorithm=c("Fixed_param"))

#print(fit.npi.regimes.priorpc,pars = c("hl","vy","optima"))
post<-extract(fit.npi.regimes.priorpc)
plot(post$Y_sim_obs[1,],dat$Y_obs)

########################################################################################################
#Milestone 4 - Posterior Predictive Check
#Regime model - regimes - painted or SIMMAP
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_reg_postpc.stan")
stan_model <- stan_model("blouchOU_reg_postpc.stan")

fit.npi.regimes.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

print(fit.npi.regimes.postpc,pars = c("hl","vy","optima"))
post<-extract(fit.npi.regimes.postpc)
plot(post$Y_sim_obs[4,],dat$Y_obs)

########################################################################################################
#Milestone 5 - Prior Posterior Check
#Regimes with direct effect model and measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_reg_direct_priorpc.stan")
stan_model <- stan_model("blouchOU_reg_direct_priorpc.stan")

fit.npi.regimes.direct.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000, algorithm=c("Fixed_param"))
#print(fit.npi.regimes.pred,pars = c("hl","vy","optima","beta"))
post<-extract(fit.npi.regimes.direct.priorpc)
plot(post$Y_sim_obs[4,],dat$Y_obs)

########################################################################################################
#Milestone 5 - Posterior Posterior Check
#Regimes with direct effect model and measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_reg_direct_postpc.stan")
stan_model <- stan_model("blouchOU_reg_direct_postpc.stan")

fit.npi.regimes.direct.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000)
#print(fit.npi.regimes.pred,pars = c("hl","vy","optima","beta"))
post<-extract(fit.npi.regimes.direct.postpc)
plot(post$Y_sim_obs[4,],dat$Y_obs)

########################################################################################################
#Milestone 6 - Prior Predictive Check
#Regimes with adaptation model with measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_reg_adapt_priorpc.stan")
stan_model <- stan_model("blouchOU_reg_adapt_priorpc.stan")
fit.npi.regimes.adapt.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2, algorithm=c("Fixed_param"))
#print(fit.npi.regimes.adapt.priorpc,pars = c("hl","vy","optima","beta","beta_e"))
post<-extract(fit.npi.regimes.adapt.priorpc)

plot(post$Y_sim_obs[4,],dat$Y_obs)

########################################################################################################
#Milestone 6 - Posterior Predictive Check
#Regimes with adaptation model with measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_reg_adapt_postpc.stan")
stan_model <- stan_model("blouchOU_reg_adapt_postpc.stan")
fit.npi.regimes.adapt.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)
#print(fit.npi.regimes.adapt.priorpc,pars = c("hl","vy","optima","beta","beta_e"))
post<-extract(fit.npi.regimes.adapt.postpc)

plot(post$Y_sim_obs[4,],dat$Y_obs)

#STOPPED
########################################################################################################
#Milestone 7 - Prior Predictive Check
#Regimes with direct and adaptation model and measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_reg_direct_adaptive_ME.stan")
stan_model <- stan_model("blouchOU_reg_direct_adaptive_ME.stan")

fit.npi.regimes.pred<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

print(fit.npi.regimes.pred,pars = c("hl","vy","optima","beta","beta_e"))
post<-extract(fit.npi.regimes.pred)
########################################################################################################
#Milestone 7 - Posterior Predictive Check
#Regimes with direct and adaptation model and measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_reg_direct_adaptive_ME.stan")
stan_model <- stan_model("blouchOU_reg_direct_adaptive_ME.stan")

fit.npi.regimes.pred<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

print(fit.npi.regimes.pred,pars = c("hl","vy","optima","beta","beta_e"))
post<-extract(fit.npi.regimes.pred)

########################################################################################################
#Milestone 8
#Multilevel model - multilevel optima
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg_mli.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_mli.stan")
#First run
fit.mli.regimes<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)
print(fit.mli.regimes,pars = c("hl","vy","optima","optima_bar","sigma"))
plot(precis(fit.mli.regimes,depth=2,pars = c("hl","vy","optima","optima_bar","sigma")))

post<-extract(fit.mli.regimes)
########################################################################################################
#Milestone 8
#Multilevel model - multilevel optima, non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg_mli_nc.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_mli_nc.stan")
#First run
fit.mli.nc.regimes<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2) #,control=list(adapt_delta=0.95))
print(fit.mli.nc.regimes,pars = c("hl","vy","optima","optima_bar","sigma"))
plot(precis(fit.mli.nc.regimes,depth=2,pars = c("hl","vy","optima","optima_bar","sigma")))
post<-extract(fit.mli.nc.regimes)
########################################################################################################
########################################################################################################
#Milestone 9
#Multilevel model - multilevel optima with direct effects predictor
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_mli_direct_ME.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_mli_direct_ME.stan")
#First run
fit.mli.regi.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)
print(fit.mli.regi.direct,pars = c("hl","vy","optima","optima_bar","beta","sigma"))
plot(precis(fit.mli.regi.direct,depth=2,pars = c("hl","vy","optima","optima_bar","sigma")))

post<-extract(fit.mli.regi.direct)

########################################################################################################
#Milestone 10
#Multilevel model - multilevel optima with adaptive predictor and measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_mli_adaptive_ME.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_mli_adaptive_ME.stan")
#First run
fit.mli.adapt.regimes<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)
print(fit.mli.adapt.regimes,pars = c("hl","vy","optima","optima_bar","beta","sigma"))
plot(precis(fit.mli.adapt.regimes,depth=2,pars = c("hl","vy","optima","optima_bar","beta","sigma")))

post<-extract(fit.mli.adapt.regimes)
########################################################################################################
#Milestone 11
#Multilevel model - multilevel optima with direct effect and adaptive predictor and measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_mli_directadaptive.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_mli_directadaptive.stan")
#First run
fit.mli.directadapt.regimes<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)
print(fit.mli.directadapt.regimes,pars = c("hl","vy","optima","optima_bar","beta","sigma"))
plot(precis(fit.mli.directadapt.regimes,depth=2,pars = c("hl","vy","optima","optima_bar","beta","sigma")))

post<-extract(fit.mli.directadapt.regimes)
########################################################################################################
#Milestone 12
#Combination of regime model with direct effect model with measurement error and varying slopes
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_direct_ME_Vs.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_direct_ME_VarSlopes.stan")
fit.reg.direct.Vs<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.direct.Vs,pars = c("hl","vy","optima_beta","beta"))
plot(precis(fit.reg.direct.Vs,depth=2,pars = c("hl","vy","optima_beta","beta")))
post<-extract(fit.reg.direct.Vs)
########################################################################################################

#Milestone 12
#Combination of regime model with multiple traits direct effect model with measurement error and varying slopes
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multidirect_ME_VarSlopes.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_multidirect_ME_VarSlopes.stan")
fit.reg.multidirect.VarSlopes<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.multidirect.VarSlopes,pars = c("hl","vy","optima_beta","beta"))
plot(precis(fit.reg.multidirect.VarSlopes,depth=3,pars = c("hl","vy","optima_beta","beta")))
post<-extract(fit.reg.multidirect.VarSlopes)

########################################################################################################
#Milestone 13
#Combination of regime model with direct effect model with measurement error and correlated varying effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_direct_ME_VarEff.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_direct_ME_VarEff.stan")
fit.reg.direct.VarEff<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.multidirect.VarEff,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta"))
plot(precis(fit.reg.multidirect.VarEff,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta")))
post<-extract(fit.reg.direct.VarEff)
########################################################################################################
#Milestone 13
#Combination of regime model with multiple traits direct effect model with measurement error and  varying effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multidirect_ME_VarEff.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_multidirect_ME_VarEff.stan")
fit.reg.multidirect.VarEff<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.multidirect.VarEff,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta"))
plot(precis(fit.reg.multidirect.VarEff,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta")))
post<-extract(fit.reg.multidirect.VarEff)
########################################################################################################
#Milestone 14
#Combination of regime model with direct effect model with mesurement error and correlated varying effects - non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_direct_ME_VarEff_nc.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_direct_ME_VarEff_nc.stan")
fit.reg.direct.VarEff<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.direct.VarEff,pars = c("hl","vy","optima","beta","sigma","optima_bar","beta_bar"))
plot(precis(fit.reg.direct.VarEff,depth=3,pars = c("hl","vy","optima","beta","Rho","sigma","optima_bar","beta_bar")))
post<-extract(fit.reg.direct.VarEff)
########################################################################################################
#Milestone 14
#Combination of regime model with multitrait direct effect model with measurement error and correlated varying effects - non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multidirect_ME_VarEff_nc.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_multidirect_ME_VarEff_nc.stan")
fit.reg.multidirect.VarEff.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.multidirect.VarEff.nc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta"))
plot(precis(fit.reg.multidirect.VarEff.nc,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta")))
post<-extract(fit.reg.multidirect.VarEff.nc)
########################################################################################################
#Milestone 15
#Combination of regime model with adaptive model with measurement error and varying slopes
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multiadaptive_ME_VarSlopes.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multiadaptive_ME_VarSlopes.stan")
stan_model <- stan_model("blouchOU_reg_multiadaptive_ME_VarSlopes.stan")
fit.reg.adapt.VarSlopes<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.adapt.VarSlopes,pars = c("hl","vy","optima","beta","beta_e"))
plot(precis(fit.reg.adapt.VarSlopes,depth=3,pars = c("hl","vy","optima","beta","beta_e")))
post<-extract(fit.reg.adapt.VarSlopes)

#Cmdstanr - setup
path<-"/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multiadaptive_ME_VarSlopes.stan"
file <- file.path(path)
mod <- cmdstan_model(file)

fit <- mod$sample(
  data = dat, 
  #seed = 10, 
  chains = 2, 
  parallel_chains = 2,
  refresh = 500 # print update every 500 iters
)
fit$summary(variables = c("hl","vy","optima","beta","beta_e"))

########################################################################################################
########################################################################################################
#Based on Milestone 16 - mlm with varying effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/Model Checking/blouchOU_reg_adapt_mlm_ve_priorpcheck.stan")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/chatgpt_blouchOU_reg_adapt_mlm_ve.stan")

stan_model <- stan_model("blouchOU_reg_adapt_mlm_ve_priorpcheck.stan")
fit.reg.adapt.mlm.ve.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000, algorithm=c("Fixed_param"))
post<-extract(fit.reg.adapt.mlm.ve.priorpc)
plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Based on Milestone 16 - mlm with varying effects
#Posterior predictive checks
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/Model Checking/blouchOU_reg_adapt_mlm_ve_postpcheck.stan")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/chatgpt_blouchOU_reg_adapt_mlm_ve.stan")

stan_model <- stan_model("blouchOU_reg_adapt_mlm_ve_postpcheck.stan")
fit.reg.adapt.mlm.ve.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
post<-extract(fit.reg.adapt.mlm.ve.postpc)
plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################

########################################################################################################
#Milestone 17
#Regime model with multiadaptive model with measurement error and varying effects - non-centered version
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multiadaptive_ME_VarEff_nc.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multiadaptive_ME_VarSlopes.stan")
stan_model <- stan_model("blouchOU_reg_multiadaptive_ME_VarEff_nc.stan")
fit.reg.adapt.VarEff.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =1000)
print(fit.reg.adapt.VarEff.nc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e"))
plot(precis(fit.reg.adapt.VarEff.nc,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e")))
post<-extract(fit.reg.adapt.VarEff.nc)
########################################################################################################
#Milestone 18
#Regime model with multi-direct and multi-adaptive model with measurement error and varying effects - non-centered version
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multidirectadaptive_ME_VarEff.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multiadaptive_ME_VarSlopes.stan")
stan_model <- stan_model("blouchOU_reg_multidirectadaptive_ME_VarEff.stan")
fit.reg.multidirectadapt.VarEff<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =1000)
print(fit.reg.multidirectadapt.VarEff,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e"))
plot(precis(fit.reg.multidirectadapt.VarEff,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e")))
post<-extract(fit.reg.multidirectadapt.VarEff)

########################################################################################################
#Milestone 19
#Regime model with multi-direct and multi-adaptive model with measurement error and varying effects - non-centered version
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multidirectadaptive_ME_VarEff_nc.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multiadaptive_ME_VarSlopes.stan")
stan_model <- stan_model("blouchOU_reg_multidirectadaptive_ME_VarEff_nc.stan")
fit.reg.multidirectadapt.VarEff.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =1000)
print(fit.reg.multidirectadapt.VarEff.nc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e"))
plot(precis(fit.reg.multidirectadapt.VarEff.nc,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e")))
post<-extract(fit.reg.multidirectadapt.VarEff.nc)
########################################################################################################

########################################################################################################
#Milestone 1
#Using direct effect model with ME - Statistical Rethinking Version
#Also works fo multivariariate Xs
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_direct_ME_SR.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_direct_ME_SR.stan")

fit.npi.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000, cores=2)

print(fit.npi.direct,pars = c("hl","vy","alpha","beta"))
post<-extract(fit.npi.direct)

########################################################################################################
#Based on Milestone 16 - mlm with varying effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/Model Checking/blouchOU_reg_adapt_mlm_ve_priorpcheck.stan")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/chatgpt_blouchOU_reg_adapt_mlm_ve.stan")

stan_model <- stan_model("blouchOU_reg_adapt_mlm_ve_priorpcheck.stan")
fit.reg.adapt.mlm.ve.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000, algorithm=c("Fixed_param"))
print(fit.reg.adapt.mlm.ve.priorpc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e"))
plot(precis(fit.reg.adapt.mlm.ve.priorpc,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e")))
post<-extract(fit.reg.adapt.mlm.ve.priorpc)
plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################
#Posterior predictive checks
#Based on Milestone 16 - mlm with varying effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/Model Checking/blouchOU_reg_adapt_mlm_ve_postpcheck.stan")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/chatgpt_blouchOU_reg_adapt_mlm_ve.stan")

stan_model <- stan_model("blouchOU_reg_adapt_mlm_ve_postpcheck.stan")
fit.reg.adapt.mlm.ve.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
print(fit.reg.adapt.mlm.ve.postpc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e","Y_sim_obs"))
plot(precis(fit.reg.adapt.mlm.ve.postpc,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e","Y_sim_obs"))))
post<-extract(fit.reg.adapt.mlm.ve.postpc)
plot(post$Y_sim_obs[1,],dat$Y_obs)
########################################################################################################

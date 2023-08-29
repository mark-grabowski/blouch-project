#Blouch Model Validation Code
#Run after Simulate OU XY Data.R
#Using recoded versions of Blouch


#library(cmdstanr)
#check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)

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
#Milestone 1
#Using direct effect model w/o ME from Hansen (1997) for V/CV matrix
#Works for multiple traits
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_direct.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_direct.stan")

fit.npi.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter = 2000)

print(fit.npi.direct,pars = c("hl","vy","alpha","beta"))
post<-extract(fit.npi.direct)

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
#Milestone 2
#Using adaptive model w/o ME
#Also works for multivariate Xs
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_adaptive.stan") #Macbook Pro
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_adaptive.stan")
stan_model <- stan_model("blouchOU_adaptive.stan")

fit.npi.adaptive<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000)

print(fit.npi.adaptive,pars = c("hl","vy","alpha","beta","beta_e"))
post<-extract(fit.npi.adaptive)

########################################################################################################
#Milestone 2
#Using adaptive model with ME - Statistical Rethinking Version
#Also works fo multivariariate Xs
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_adaptive_ME_SR.stan") #Macbook Pro
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_adaptive_ME_SR.stan")
stan_model <- stan_model("blouchOU_adaptive_ME_SR.stan")

fit.npi.adaptive<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000)

print(fit.npi.adaptive,pars = c("hl","vy","alpha","beta","beta_e"))
post<-extract(fit.npi.adaptive)

########################################################################################################
#Milestone 3
#Combination of direct effect and adaptive predictors w/o ME
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_direct_adaptive.stan")

stan_model <- stan_model("blouchOU_direct_adaptive.stan")

fit.npi.mixed<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000)

print(fit.npi.mixed,pars = c("hl","vy","alpha","beta","beta_e"))
post<-extract(fit.npi.mixed)

########################################################################################################
#Milestone 3
#Combination of direct effect and adaptive predictors with ME
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_direct_adaptive_ME.stan")

stan_model <- stan_model("blouchOU_direct_adaptive_ME.stan")

fit.npi.mixed<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.npi.mixed,pars = c("hl","vy","alpha","beta","beta_e"))
post<-extract(fit.npi.mixed)

########################################################################################################
#Milestone 4
#Regime model - regimes - painted or SIMMAP
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached//blouchOU_reg.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")
stan_model <- stan_model("blouchOU_reg.stan")

fit.npi.regimes<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

print(fit.npi.regimes,pars = c("hl","vy","optima"))
post<-extract(fit.npi.regimes)
########################################################################################################
#Milestone 5
#Regimes with direct effect model without measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg_direct.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")
stan_model <- stan_model("blouchOU_reg_direct.stan")

fit.npi.regimes.pred<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000)

print(fit.npi.regimes.pred,pars = c("hl","vy","optima","beta"))
post<-extract(fit.npi.regimes.pred)
########################################################################################################
#Milestone 5
#Regimes with direct effect model and measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_reg_direct_ME.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")
stan_model <- stan_model("blouchOU_reg_direct_ME.stan")

fit.npi.regimes.pred<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

print(fit.npi.regimes.pred,pars = c("hl","vy","optima","beta"))
post<-extract(fit.npi.regimes.pred)
########################################################################################################
#Milestone 6
#Regimes with adaptation model without measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg_adaptive.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")
stan_model <- stan_model("blouchOU_reg_adaptive.stan")

fit.npi.regimes.pred<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000)

print(fit.npi.regimes.pred,pars = c("hl","vy","optima","beta","beta_e"))
post<-extract(fit.npi.regimes.pred)

########################################################################################################
#Milestone 6
#Regimes with adaptation model with measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg_adaptive_ME.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")
stan_model <- stan_model("blouchOU_reg_adaptive_ME.stan")

fit.npi.regimes.pred<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

print(fit.npi.regimes.pred,pars = c("hl","vy","optima","beta","beta_e"))
post<-extract(fit.npi.regimes.pred)

########################################################################################################
#Milestone 7
#Regimes with direct and adaptation model without measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_reg_direct_adaptive.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")
stan_model <- stan_model("blouchOU_reg_direct_adaptive.stan")

fit.npi.regimes.pred<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,cores=2)

print(fit.npi.regimes.pred,pars = c("hl","vy","optima","beta","beta_e"))
post<-extract(fit.npi.regimes.pred)
########################################################################################################
#Milestone 7
#Regimes with direct and adaptation model and measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_reg_direct_adaptive_ME.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")
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
#Milestone 16-S66-S68
#Regime model with adaptive model with measurement error and varying effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multiadaptive_ME_VarEff.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_multiadaptive_ME_VarSlopes.stan")
stan_model <- stan_model("blouchOU_reg_multiadaptive_ME_VarEff.stan")
fit.reg.adapt.VarEff<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =1000)
print(fit.reg.adapt.VarEff,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e"))
plot(precis(fit.reg.adapt.VarEff,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e")))
post<-extract(fit.reg.adapt.VarEff)
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
#Milestone XX - partially pooled multi-simmap code
setwd("~/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("~/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_multisimmap.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_multisimmap.stan")
#First run
fit.ms.regimes<- rstan::sampling(object = stan_model,data = dat,chains = 1,cores=1,iter =1000)
print(fit.ms.regimes,pars = c("hl_bar","hl_sigma","vy_bar","vy_sigma","optima_bar","optima_sigma"))

post<-extract(fit.ms.regimes)
########################################################################################################
#Milestone XX - partially pooled multi-simmap code, non-centered version
setwd("~/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("~/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_multisimmap_nc.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_multisimmap_nc.stan")
#First run
fit.ms.nc.regimes<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=4,iter =2000)
print(fit.ms.nc.regimes,pars = c("hl","vy","optima_bar","sigma"))

post<-extract(fit.ms.nc.regimes)
########################################################################################################






#Regimme/SIMMAP model with multi-simmap tree
path<-"~/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_multisimmap.stan"
file <- file.path(path)
mod <- cmdstan_model(file)











#SR V/CV
########################################################################################################
#Testing multi-level regime model with Statistical Rethinking approach to calculating V/CV
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg_mli_SR.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_mli_SR.stan")

fit.mli.regimes<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter = 2000, cores=2)

print(fit.mli.regimes,pars = c("hl","var_anc","optima","optima_bar","sigma"))
plot(precis(fit.mli.regimes,depth=2,pars = c("hl","var_anc","optima","optima_bar","sigma")))
post<-extract(fit.mli.regimes)
#pairs(fit.mli.regimes,pars = c("hl","var_anc","optima","optima_bar","sigma"))

########################################################################################################
#Testing multi-level regime model with Statistical Rethinking approach to calculating V/CV
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg_mli_nc_SR.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_reg_mli_nc_SR.stan")

fit.mli.nc.regimes<- rstan::sampling(object = stan_model,data = dat,chains = 1,iter = 2000, cores=1)

print(fit.mli.nc.regimes,pars = c("hl","var_anc","optima","optima_bar","sigma"))
plot(precis(fit.mli.nc.regimes,depth=2,pars = c("hl","var_anc","optima","optima_bar","sigma")))
post<-extract(fit.mli.nc.regimes)
#pairs(fit.mli.nc.regimes,pars = c("hl","var_anc","optima","optima_bar","sigma"))
########################################################################################################












#Milestone 8
#SIMMAP Regimes - single tree
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Milestones Reached/blouchOU_reg.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")
stan_model <- stan_model("blouchOU_reg.stan")

fit.npi.SIMMAP<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000)

print(fit.npi.SIMMAP,pars = c("hl","vy","optima","sigma2_y"))
post<-extract(fit.npi.SIMMAP)















########################################################################################################
#Milestone 4
#Regime model with direct effect variables
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_direct_adaptive_ME.stan")

stan_model <- stan_model("blouchOU_direct_adaptive_ME.stan")

fit.npi.mixed<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.npi.mixed,pars = c("hl","vy","alpha","beta","beta_e","sigma2_y"))
########################################################################################################
#Milestone 4
#Regime model with adaptive variables
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_direct_adaptive_ME.stan")

stan_model <- stan_model("blouchOU_direct_adaptive_ME.stan")

fit.npi.mixed<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000,control=list(adapt_delta=0.80))

print(fit.npi.mixed,pars = c("hl","vy","alpha","beta","beta_e","sigma2_y"))
########################################################################################################



#Using Simulated Regime data from: SBR1 - Subsample Primate Regimes

#Basic direct stan model - no regimes
intercept_test<-rnorm(100,stan_sim_data$ols_intercept,0.1)
slope_test<-rnorm(100,stan_sim_data$ols_slope,0.1)

setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/1. Standard MV Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/1. Standard MV Versions/blouchOU_v1_9.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_v1_9.stan")

fit.npi.direct<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.80))

print(fit.npi.direct,pars = c("hl","beta","beta_e","vy"))

#For downstream analysis and plots




Basic adaptive stan model - no regimes
Priors
intercept_test<-rnorm(100,stan_sim_data$ols_intercept,0.3)
slope_test<-rnorm(100,stan_sim_data$ols_slope,0.1)

```{r}
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/1. Standard MV Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/1. Standard MV Versions/blouchOU_v1_9.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_v1_9.stan")

fit.npi.adaptive<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.80))

print(fit.npi.adaptive,pars = c("hl","beta","beta_e","vy"))

#For downstream analysis and plots
```


Basic no varying effects regimes model
```{r}
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/1. Standard MV Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/1. Standard MV Versions/blouchOUReg_v1_5.stan") #Macbook Pro

stan_model <- stan_model("blouchOUReg_v1_5.stan")

#setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_v1_6.stan")

#stan_model <- stan_model("blouchOUReg_v1_6.stan")

#setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_v1_7.stan")

#stan_model <- stan_model("blouchOUReg_v1_7.stan")


fit.npi.regimes<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.80))

print(fit.npi.regimes,pars = c("hl","beta","beta_e","vy"))

#For downstream analysis and plots
```





Plotting
```{r}
library(rethinking)

stan_dens(fit.npi.regimes,pars = c("hl","beta","vy"))
trankplot(fit.npi.regimes,pars = c("hl","beta","vy"))
plot(precis(fit.npi.regimes,depth=1,pars = c("hl","beta","vy")))

#For downstream analysis and plots
ext.fit <- rstan::extract(fit.npi.regimes)


#Plot correlation between posterior distributions
library(ggplot2)
hl<-ext.fit$hl
vy<-ext.fit$vy
beta<-ext.fit$beta

new.data<-data.frame(hl,vy,beta)
ggplot(new.data,aes(x=hl,y=vy))+
  geom_point()

ggplot(new.data,aes(x=hl,y=beta))+
  geom_point()

pairs(fit.npi.regimes,pars=c("hl","beta","vy"))

```



Multilevel model - Varying Intercepts, centered priors
```{r}
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI.stan") #Macbook Pro

#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")#Mac Studio
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI.stan") #Mac Studio

stan_model <- stan_model("blouchOUReg_MLI.stan")

fit.mli.regimes<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.80))

print(fit.mli.regimes,pars = c("hl","beta_bar","sigma","beta","beta_e","vy"))

```


Plotting
```{r}
stan_dens(fit.mli.regimes,pars = c("hl","beta_bar","sigma","beta","beta_e","vy"))

#For downstream analysis and plots
ext.fit.adaptive.test <- rstan::extract(fit.mli.regimes)

library(ggplot2)
hl<-ext.fit.adaptive.test$hl
vy<-ext.fit.adaptive.test$vy
beta<-ext.fit.adaptive.test$beta

new.data<-data.frame(hl,vy,beta)
ggplot(new.data,aes(x=hl,y=vy))+
  geom_point()

ggplot(new.data,aes(x=hl,y=beta))+
  geom_point()


pairs(fit.mli.regimes,pars=c("hl","beta","vy"))

```

Multilevel model - Varying Intercepts, non-centered priors
```{r}
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI_nc.stan") #Macbook Pro

#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")#Mac Studio
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI.stan") #Mac Studio

stan_model <- stan_model("blouchOUReg_MLI_nc.stan")

fit.mli.nc.regimes<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.9))

print(fit.mli.nc.regimes,pars = c("hl","beta_bar","sigma","beta","beta_e","vy","log_lik"))

#For downstream analysis and plots
```



Multilevel model - Varying Effects, centered priors
```{r}
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLM.stan") #Macbook Pro

#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")#Mac Studio
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI.stan") #Mac Studio

stan_model <- stan_model("blouchOUReg_MLM.stan")

fit.mli.regimes<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.80))

print(fit.mli.regimes,pars = c("hl","beta_bar","sigma","beta","beta_e","vy","log_lik"))

#For downstream analysis and plots
```


Plotting
```{r}
stan_dens(fit.mli.regimes,pars = c("hl","beta_bar","sigma","beta","beta_e","vy"))

#For downstream analysis and plots
ext.fit.adaptive.test <- rstan::extract(fit.mli.regimes)

library(ggplot2)
hl<-ext.fit.adaptive.test$hl
vy<-ext.fit.adaptive.test$vy
beta<-ext.fit.adaptive.test$beta

new.data<-data.frame(hl,vy,beta)
ggplot(new.data,aes(x=hl,y=vy))+
  geom_point()

ggplot(new.data,aes(x=hl,y=beta))+
  geom_point()


pairs(fit.mli.regimes,pars=c("hl","beta","vy"))

```

Multilevel model - Varying Effects, non-centered priors
```{r}
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI_nc.stan") #Macbook Pro

#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")#Mac Studio
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI.stan") #Mac Studio

stan_model <- stan_model("blouchOUReg_MLI_nc.stan")

fit.mli.nc.regimes<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.9))

print(fit.mli.nc.regimes,pars = c("hl","beta_bar","sigma","beta","beta_e","vy","log_lik"))

#For downstream analysis and plots
```

Model Comparisons
```{r}
log_lik_M.mli<-extract_log_lik(fit.mli.regimes, merge_chains = FALSE)
log_lik_M.npi<-extract_log_lik(fit.npi.regimes, merge_chains = FALSE)

log_ratios <- -1 * log_lik_M.mli
r_eff <- relative_eff(exp(-log_ratios))
psis_result <- psis(log_ratios, r_eff = r_eff)
str(psis_result)
plot(psis_result)


r_eff <- relative_eff(exp(log_lik_M.mli), cores = 2) 
loo_1 <- loo(log_lik_M.mli, r_eff = r_eff, cores = 2)

r_eff <- relative_eff(exp(log_lik_M.npi), cores = 2) 
loo_2 <- loo(log_lik_M.npi, r_eff = r_eff, cores = 2)

loo_compare(waic(log_lik_M.mli),waic(log_lik_M.npi))
```




Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

stan.OU.data.setup<-function(OUBMparameters,tree,n){
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  name.response.trait<-c("trait_1")
  names.direct.traits<-c("trait_2")
  stan.data<-list()
  for(i in 1:n){
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
    OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
    OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
    #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X
    stan.data[[length(stan.data)+1]]<-blouchOU.setup.mv(OUBM.trdata,name.response.trait,names.direct.traits)
  }
  return(stan.data)
}

stan.BM.data.setup<-function(BMparameters,tree,n){
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  name.response.trait<-c("trait_1")
  names.direct.traits<-c("trait_2")
  stan.data<-list()
  for(i in 1:n){
    BMdata<-simulBMProcPhylTree(cervid.trdata$phy,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
    BMdata<-data.frame(BMdata,"tip.names"=row.names(BMdata))
    BM.trdata <- make.treedata(tree, BMdata,name_column="tip.names") #Combine data with 
    #BM.trdata$dat$trait_2<-BM.trdata$dat$trait_2-mean(BM.trdata$dat$trait_2) #Mean standardize X
    stan.data[[length(stan.data)+1]]<-blouchOU.setup.v1(BM.trdata,names.traits)
  }
  return(stan.data)
}

stan.OU.simulation.reps<-function(OUBMparameters,tree,n){
  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOU_v1_7.stan")
  stan_model <- stan_model("blouchOU_v1_7.stan")
  
  #######################
  stan_data<-stan.OU.data.setup(OUBMparameters,tree,n)
  #################
  saved.parameters<-data.frame(matrix(NA,n,16))
  parameter.names<-c("direct.reg.median","direct.reg.mean","direct.reg.025","direct.reg.975",
                     "t12.median","t12.mean","t12.025","t12.975","vy.median","vy.mean","vy.025","vy.975","sigma2y.median","sigma2y.mean","sigma2y.025","sigma2y.975")
  names(saved.parameters)<-parameter.names
  #######################
  saved.parameters.mode<-data.frame(matrix(NA,n,4))
  parameter.names<-c("direct.reg.mode","t12.mode","vy.mode","sigma2y.mode")
  names(saved.parameters.mode)<-parameter.names
  ###############
  
  for (i in 1:n){
    print(paste("Iteration ",i))
    stan_data_sub<-stan_data[[i]]
    stan_data_sub<-stan_data_sub[[1]]
    #fit.random<- stan(file = stan.model,data = stan_data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
    fit<- rstan::sampling(object = stan_model,data = stan_data_sub,chains = 1,iter = 2000,verbose=F,control=list(adapt_delta=0.975))
    par.summary <- summary(fit, pars = c("beta","hl","vy","sigma2_y"), probs = c(0.025,0.5, 0.975))$summary
    print(par.summary)
    par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
    saved.parameters[i,1:4]<-par.summary[2,]
    saved.parameters[i,5:8]<-par.summary[3,]
    saved.parameters[i,9:12]<-par.summary[4,]
    saved.parameters[i,13:16]<-par.summary[5,]

    fit.ext<-extract(fit)
    saved.parameters.mode[i,1]<-Mode(fit.ext$beta[2])
    saved.parameters.mode[i,2]<-Mode(fit.ext$hl)
    saved.parameters.mode[i,3]<-Mode(fit.ext$vy)
    saved.parameters.mode[i,4]<-Mode(fit.ext$sigma2_y)
    }
  return(list(saved.parameters,saved.parameters.mode))
  }

stan.BM.simulation.reps<-function(BMparameters,tree,n){
  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1/blouchBM_v1.stan")
  stan.model <- "blouchBM_v1.stan"
  #######################
  stan.data<-stan.BM.data.setup(BMparameters,tree,n)
  #######################
  saved.parameters<-data.frame(matrix(NA,n,8))
  parameter.names<-c("evol.reg.median","evol.reg.mean","evol.reg.025","evol.reg.975","sigma2y.median","sigma2y.mean","sigma2y.025","sigma2y.975")
  names(saved.parameters)<-parameter.names
  ###############
  saved.parameters.mode<-data.frame(matrix(NA,n,2))
  parameter.names<-c("evol.reg.mode","sigma2y.mode")
  names(saved.parameters)<-parameter.names
  ###############
  
  for (i in 1:n){
    print(paste("Iteration ",i))
    stan.adaptive.data<-stan.data[[i]]
    stan.constraint.data<-stan.adaptive.data[[1]]
    stan.adaptive.data<-stan.adaptive.data[[2]]
    fit.BM<- stan(file = stan.model,data = stan.constraint.data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
    
    par.summary <- summary(fit.BM, pars = c("beta","sigma2_y"), probs = c(0.025,0.5, 0.975))$summary
    par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
    saved.parameters[i,1:4]<-par.summary[1,]
    saved.parameters[i,5:8]<-par.summary[2,]

    fit.BM.ext<-extract(fit.BM)
    saved.parameters.mode[i,1]<-Mode(fit.BM.ext$beta)
    saved.parameters.mode[i,2]<-Mode(fit.BM.ext$sigma2_y)

      }
  return(list(saved.parameters,saved.parameters.mode))
}

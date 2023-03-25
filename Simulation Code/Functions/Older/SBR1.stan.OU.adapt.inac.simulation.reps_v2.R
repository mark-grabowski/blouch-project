#SBR1 OU Simulation Code - SBR1.stan.OU.adapt.inac.simulation.reps.R
#Calculates bias, precision, and innacuracy in optimal slopes and half-lives
#Works with blouchOU_v1_7.stan and blouchOU.setup.mv_v1_3.R



Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

stan.OU.data.setup<-function(OUBMparameters,tree,n){
  #source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  #names.traits<-c("trait_1","trait_2",NA,NA)
  name.response.trait<-c("trait_1")
  names.random.traits<-c("trait_2")
  
  stan.data<-list()
  for(i in 1:n){
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
    OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
    OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
    OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X SBR1
    OUBM.trdata$dat$trait_1<-OUBM.trdata$dat$trait_1-mean(OUBM.trdata$dat$trait_1) #Mean standardize Y SBR1
    stan.data[[length(stan.data)+1]]<-blouchOU.setup.mv(OUBM.trdata,name.response.trait,names.direct.traits=NULL,names.random.traits)
    
    
    }
  return(stan.data)
}

stan.BM.data.setup<-function(BMparameters,tree,n){
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1/Blouch Setup Files/blouchOU.setup.v1.R")
  #n<-50
  names.traits<-c("trait_1","trait_2",NA,NA)
  stan.data<-list()
  for(i in 1:n){
    BMdata<-simulBMProcPhylTree(cervid.trdata$phy,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
    BMdata<-data.frame(BMdata,"tip.names"=row.names(BMdata))
    BM.trdata <- make.treedata(tree, BMdata,name_column="tip.names") #Combine data with 
    BM.trdata$dat$trait_2<-BM.trdata$dat$trait_2-mean(BM.trdata$dat$trait_2) #Mean standardize X SBR1
    BM.trdata$dat$trait_1<-BM.trdata$dat$trait_1-mean(BM.trdata$dat$trait_1) #Mean standardize Y SBR1
    stan.data[[length(stan.data)+1]]<-blouchOU.setup.v1(BM.trdata,names.traits)
  }
  return(stan.data)
}

stan.OU.simulation.reps<-function(tree,n,opt.slope,hl,vy,vY0,vX0){
  #Simulations setup to allow for OUBMparameters to be estimated in simulation - SBR1
  
  #############################
  #Setup for simulations
  #1% - Added in sigma2y from Hansen et al. 2008 - vy scaled by A
  #A = log(2)/0.01 = 69.31472
  #B = -A * slope = -69.31472 * 0.25 = -17.32868
  
  Sxx <- 1
  A<-log(2)/hl
  B<- -A*opt.slope
  Syy<-sqrt(vy*(2*A))
  
  tip.num<-length(tree$tip.label)
  times <- ape::node.depth.edgelength(tree)
  T.term <- times[1:tip.num]              
  
  #Evolutionary regression based on optimal regression and A  
  rho<- (1-(1 - exp(-A * T.term))/(A * T.term))
  evol.slope<-(rho*opt.slope)[1]
    
  OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                       B=matrix(B,ncol=1,nrow=1),mPsi=matrix(0,ncol=1,nrow=1),
                       Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                       Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))

  ############################################
  #setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  #stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOU_v1_8.stan")

  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOU_v1_8.stan")
  stan_model <- stan_model("blouchOU_v1_8.stan")
  
  #######################
  stan_data<-stan.OU.data.setup(OUBMparameters,tree,n)
  #################
  saved.parameters<-data.frame(matrix(NA,n,20))
  parameter.names<-c("opt.reg","evol.reg","hl","vy","opt.reg.median","opt.reg.mean","opt.reg.025","opt.reg.975",
                     "evol.reg.median","evol.reg.mean","evol.reg.025","evol.reg.975",
                     "hl.median","hl.mean","hl.025","hl.975","vy.median","vy.mean","vy.025","vy.975")
  names(saved.parameters)<-parameter.names
  #######################
  saved.bpa<-data.frame(matrix(NA,1,10))
  parameter.names<-c("opt.slope","evol.slope","hl","vy","opt.bias","opt.imp","opt.inac","hl.bias","hl.imp","hl.inac")
  names(saved.bpa)<-parameter.names
  #######################
  
  

  for (i in 1:n){
    print(paste("Iteration ",i,"hl",hl))
    stan_data_sub<-stan_data[[i]]
    stan_data_sub<-stan_data_sub[[1]]
    #fit.random<- stan(file = stan.model,data = stan_data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
    fit<- rstan::sampling(object = stan_model,data = stan_data_sub,chains = 1,iter = 2000,verbose=FALSE,control=list(adapt_delta=0.95))
    par.summary <- summary(fit, pars = c("beta", "beta_e","hl","vy"), probs = c(0.025,0.5, 0.975))$summary
    par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
    saved.parameters[i,1]<-opt.slope
    saved.parameters[i,2]<-evol.slope
    saved.parameters[i,3]<-hl
    saved.parameters[i,4]<-vy

    saved.parameters[i,5:8]<-par.summary[2,]
    saved.parameters[i,9:12]<-par.summary[3,]
    saved.parameters[i,13:16]<-par.summary[4,]
    saved.parameters[i,17:20]<-par.summary[5,]
    }
  #Using median as measurement
  saved.bpa[1,1]<-opt.slope
  saved.bpa[1,2]<-evol.slope
  saved.bpa[1,3]<-hl
  saved.bpa[1,4]<-vy
  saved.bpa[1,5]<-mean(saved.parameters[,5]-opt.slope)
  saved.bpa[1,6]<-var(saved.parameters[,5])
  saved.bpa[1,7]<-mean((saved.parameters[,5]-opt.slope)^2)
  saved.bpa[1,8]<-mean(saved.parameters[,13]-hl)
  saved.bpa[1,9]<-var(saved.parameters[,13])
  saved.bpa[1,10]<-mean((saved.parameters[,13]-hl)^2)
  
  return(list(saved.parameters,saved.bpa))
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
    stan.adaptive.data<-
    fit.random<- stan(file = stan.model,data = stan.adaptive.data,chains = 1,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
    
    par.summary <- summary(fit.random, pars = c("beta","sigma2_y"), probs = c(0.025,0.5, 0.975))$summary
    par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
    saved.parameters[i,1:4]<-par.summary[1,]
    saved.parameters[i,5:8]<-par.summary[2,]

    fit.random.ext<-extract(fit.random)
    saved.parameters.mode[i,1]<-Mode(fit.random.ext$beta)
    saved.parameters.mode[i,2]<-Mode(fit.random.ext$sigma2_y)

      }
  return(list(saved.parameters,saved.parameters.mode))
}

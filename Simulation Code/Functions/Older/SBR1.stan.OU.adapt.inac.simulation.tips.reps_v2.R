stan.OU.tree.data.setup<-function(OUBMparameters,tips,n,A,opt.slope,sym){
  #source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  stan.data<-list()
  name.response.trait<-c("trait_1")
  names.random.traits<-c("trait_2")
  
  
  for(i in 1:n){
    tips.evol.slope<-NULL
    if(sym != 1){
      #Following Uyeda and Harmon, 2014
      #tree<-sim.bd.taxa(tips, numbsim=1, lambda=0.1, mu=0, frac = 1, complete = TRUE,
       #               stochsampling = FALSE)[[1]]
      
      #Following bayou tutorial: https://github.com/uyedaj/bayou/blob/master/tutorial.md
      tree <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = tips)
      #tree <- sim.bdtree(b = 1, d = 0.1, stop = "taxa", n = tips)
      
    }
    else{
      #Ultrametric tree - sym = 1
      tree <- sim.taxa(1, n=tips, waitsp=function()2,symmetric=TRUE, complete=TRUE)
      tree<-tree[[1]]
      }
    l.tree<-max(branching.times(tree))
    tree$edge.length<-tree$edge.length/l.tree ## rescale tree to height 1
    
    n<-length(tree$tip.label)
    mrca1 <- ape::mrca(tree)
    times <- ape::node.depth.edgelength(tree)
    T.term <- times[1:n]              
    
    rho<- (1-(1 - exp(-A * T.term))/(A * T.term))
    tips.evol.slope<-c(tips.evol.slope,(rho*opt.slope)[1])
    
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
    OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
    OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
    stan.data[[length(stan.data)+1]]<-c(blouchOU.setup.mv(OUBM.trdata,name.response.trait,names.direct.traits=NULL,names.random.traits),tips.evol.slope)
    
  }
  
  return(stan.data)
}

stan.BM.tree.data.setup<-function(BMparameters,tips,n,sym){
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  #n<-50
  name.response.trait<-c("trait_1")
  names.random.traits<-c("trait_2")
  stan.data<-list()
  for(i in 1:n){
    if(sym != 1){
      #Following Uyeda and Harmon, 2014
      #tree<-sim.bd.taxa(tips, numbsim=1, lambda=0.1, mu=0, frac = 1, complete = TRUE,
      #               stochsampling = FALSE)[[1]]
      
      #Following bayou tutorial: https://github.com/uyedaj/bayou/blob/master/tutorial.md
      tree <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = tips)
    }
    else{
      #Ultrametric tree
      tree <- sim.taxa(1, n=tips, waitsp=function()2,symmetric=TRUE, complete=TRUE)
      tree<-tree[[1]]
    }

    l.tree<-max(branching.times(tree))
    tree$edge.length<-tree$edge.length/l.tree ## rescale tree to height 1
    
    #########################
    BMdata<-simulBMProcPhylTree(tree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
    BMdata<-data.frame(BMdata,"tip.names"=row.names(BMdata))
    BM.trdata <- make.treedata(tree, BMdata,name_column="tip.names") #Combine data with 
    #BM.trdata$dat$trait_2<-BM.trdata$dat$trait_2-mean(BM.trdata$dat$trait_2) #Mean standardize X
    stan.data[[length(stan.data)+1]]<-c(blouchOU.setup.mv(BM.trdata,name.response.trait,names.direct.traits=NULL,names.random.traits),tips.evol.slope)
  }
  return(stan.data)
}

stan.OU.simulation.tips.reps<-function(tips.list,n,opt.slope,sym,hl,vy,vY0,vX0){
  #OUBM Parameters calculated in code - SBR1
  library(TreeSimGM)
  #############################
  #Setup for simulations
  #1% - Added in sigma2y from Hansen et al. 2008 - vy scaled by A
  #A = log(2)/0.01 = 69.31472
  #B = -A * slope = -69.31472 * 0.25 = -17.32868
  
  Sxx <- 1
  A<-log(2)/hl
  B<- -A*opt.slope
  Syy<-sqrt(vy*(2*A))
  
  OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                       B=matrix(B,ncol=1,nrow=1),mPsi=matrix(0,ncol=1,nrow=1),
                       Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                       Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
  
  ############################################
  
  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOU_v1_8.stan")
  #setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  #stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOU_v1_7.stan")

  stan_model <- stan_model("blouchOU_v1_8.stan")
  

  #######################
  #################
  saved.parameters<-data.frame(matrix(NA,n,22))
  all.saved.parameters<-NULL
  
  parameter.names<-c("tips","sym","opt.reg","evol.reg","hl","vy","opt.reg.median","opt.reg.mean","opt.reg.025","opt.reg.975",
                     "evol.reg.median","evol.reg.mean","evol.reg.025","evol.reg.975",
                     "hl.median","hl.mean","hl.025","hl.975","vy.median","vy.mean","vy.025","vy.975")
  names(saved.parameters)<-parameter.names

    #######################
  saved.bpa<-data.frame(matrix(NA,length(tips.list),12))
  parameter.names<-c("tips","sym","opt.slope","evol.slope","hl","vy","opt.bias","opt.imp","opt.inac","hl.bias","hl.imp","hl.inac")
  names(saved.bpa)<-parameter.names
  #######################
  
  for(j in 1:length(tips.list)){
    stan_data<-stan.OU.tree.data.setup(OUBMparameters,tips.list[j],n,as.numeric(OUBMparameters$A),opt.slope,sym)
    for (i in 1:n){
      print(paste("Tips",tips.list[j],"Iteration ",i))
      stan_data_sub<-stan_data[[i]]
      evol.slope<-stan_data_sub[[3]]
      #print(stan_data_sub)
      stan_data_sub<-stan_data_sub[[1]]
      
      fit<- rstan::sampling(object = stan_model,data = stan_data_sub,chains = 1,iter = 2000,verbose=F,control=list(adapt_delta=0.95))
      par.summary <- summary(fit, pars = c("beta", "beta_e","hl","vy"), probs = c(0.025,0.5, 0.975))$summary
      #print(par.summary)
      par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
      saved.parameters[i,1]<-tips.list[j]
      saved.parameters[i,2]<-paste(tips.list[j],if(sym==1){"-sym"})
      saved.parameters[i,3]<-opt.slope
      saved.parameters[i,4]<-evol.slope
      saved.parameters[i,5]<-hl
      saved.parameters[i,6]<-vy
      saved.parameters[i,7:10]<-par.summary[2,]
      saved.parameters[i,11:14]<-par.summary[3,]
      saved.parameters[i,15:18]<-par.summary[4,]
      saved.parameters[i,19:22]<-par.summary[5,]
      }
    all.saved.parameters<-rbind(all.saved.parameters,saved.parameters)
    #Using median as measurement
    saved.bpa[j,1]<-tips.list[j]
    saved.bpa[j,2]<-paste(tips.list[j],if(sym==1){"-sym"},sep="")
    
    saved.bpa[j,3]<-opt.slope
    saved.bpa[j,4]<-mean(saved.parameters[,4])
    saved.bpa[j,5]<-hl
    saved.bpa[j,6]<-vy
    saved.bpa[j,7]<-mean(saved.parameters[,7]-opt.slope)
    saved.bpa[j,8]<-var(saved.parameters[,7])
    saved.bpa[j,9]<-mean((saved.parameters[,7]-opt.slope)^2)
    saved.bpa[j,10]<-mean(saved.parameters[,15]-hl)
    saved.bpa[j,11]<-var(saved.parameters[,15])
    saved.bpa[j,12]<-mean((saved.parameters[,15]-hl)^2)
  }
  return(list(all.saved.parameters,saved.bpa))
}

stan.BM.simulation.tips.reps<-function(BMparameters,tips,n,sym){
  library(TreeSimGM)
  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1/blouchBM_v1.stan")
  stan.model <- "blouchBM_v1.stan"
  #######################
  stan.data<-stan.BM.tree.data.setup(BMparameters,tips,n,sym)
  
  #######################
  saved.parameters<-data.frame(matrix(NA,n,8))
  parameter.names<-c("evol.reg.median","evol.reg.mean","evol.reg.025","evol.reg.975","sigma2y.median","sigma2y.mean","sigma2y.025","sigma2y.975")
  names(saved.parameters)<-parameter.names
  
  #######################
  saved.parameters.mode<-data.frame(matrix(NA,n,2))
  parameter.names<-c("evol.reg.mode","sigma2y.mode")
  names(saved.parameters.mode)<-parameter.names
  ###############
  
  for (i in 1:n){
    print(paste("Iteration ",i))
    stan.adaptive.data<-stan.data[[i]]
    stan.adaptive.data<-stan.adaptive.data[[2]]
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

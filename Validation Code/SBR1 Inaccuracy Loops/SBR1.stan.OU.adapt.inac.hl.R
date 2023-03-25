#SBR1 OU - Calculating Inacuracy in hl code
#1. Calculate innacuracies given a list of hl values to explore - 10 reps per half life value
#Calculates bias, precision, and innacuracy in optimal slopes and half-lives
#Works with blouchOU_v1_7.stan and blouchOU.setup.mv_v1_3.R
#Adaptive model including estimation of Ya and b0

stan.adapt.OU.data.setup<-function(OUBMparameters,tree,N){
  #source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  name.response.trait<-c("trait_1")
  names.random.traits<-c("trait_2")
  
  stan.data<-list()
  for(i in 1:N){
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
    OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
    OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
    #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X SBR1 - not using because estimating Ya
    #OUBM.trdata$dat$trait_1<-OUBM.trdata$dat$trait_1-mean(OUBM.trdata$dat$trait_1) #Mean standardize Y SBR1 - not using because estimating Ya
    stan.data[[length(stan.data)+1]]<-blouchOU.setup.mv(OUBM.trdata,name.response.trait,names.direct.traits=NULL,names.random.traits)
    
    
    }
  return(stan.data)
}

stan.direct.OU.data.setup<-function(OUBMparameters,tree,N){
  #source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  name.response.trait<-c("trait_1")
  names.direct.traits<-c("trait_2")
  
  stan.data<-list()
  for(i in 1:N){
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
    OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
    OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
    #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X SBR1
    #OUBM.trdata$dat$trait_1<-OUBM.trdata$dat$trait_1-mean(OUBM.trdata$dat$trait_1) #Mean standardize Y SBR1
    stan.data[[length(stan.data)+1]]<-blouchOU.setup.mv(OUBM.trdata,name.response.trait,names.direct.traits)
    
    
  }
  return(stan.data)
}

stan.BM.data.setup<-function(BMparameters,tree,N){
  #source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  name.response.trait<-c("trait_1")
  names.random.traits<-c("trait_2")
  
  stan.data<-list()
  for(i in 1:N){
    BMdata<-simulBMProcPhylTree(tree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
    
    BMdata<-data.frame(BMdata,"tip.names"=row.names(BMdata))
    BM.trdata <- make.treedata(tree, BMdata,name_column="tip.names") #Combine data with 
    #BM.trdata$dat$trait_2<-BM.trdata$dat$trait_2-mean(BM.trdata$dat$trait_2) #Mean standardize X SBR1
    #BM.trdata$dat$trait_1<-BM.trdata$dat$trait_1-mean(BM.trdata$dat$trait_1) #Mean standardize Y SBR1
    stan.data[[length(stan.data)+1]]<-blouchOU.setup.mv(BM.trdata,name.response.trait,names.direct.traits=NULL,names.random.traits)
  }
  return(stan.data)
}
create.tree<-function(tips,sym){
    if(sym != 1){
      tree <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = tips)
    }
    else{
      #Ultrametric tree - sym = 1, symmetric
      tree <- sim.taxa(1, n=tips, waitsp=function()2,symmetric=TRUE, complete=TRUE)
      tree<-tree[[1]]
    }
  l.tree<-max(branching.times(tree))
  tree$edge.length<-tree$edge.length/l.tree ## rescale tree to height 1
  tree$edge.length[tree$edge.length==0] <- .Machine$double.eps
  
  return(tree)
  }

subsample.tree<-function(tips,tree){
  #sampled.tips<-length(tree$tip.label)-tips
  #print(length(tree$tip.label)-tips)
  tree <- keep.tip(tree,sample(tree$tip.label)[1:tips]) 
  
  l.tree<-max(branching.times(tree))
  tree$edge.length<-tree$edge.length/l.tree ## rescale tree to height 1
  tree$edge.length[tree$edge.length==0] <- .Machine$double.eps
  return(tree)
}


stan.adapt.OU.inac.hl<-function(tips.list=NULL,tree=NULL,N,opt.slope,hl.list,vy,vY0,vX0,theta,sym=0){
  #Using adaptive OU model, allows for set tree, sampled trees with a certain number of tips, and new trees with a certain number of tips
  #Simulations setup to allow for OUBMparameters to be estimated in simulation - SBR1
  
  ############################################
  #setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  #stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOU_v1_8.stan")

  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOU_v1_9.stan")
  stan_model <- stan_model("blouchOU_v1_9.stan")
  
  #######################
  #################
  saved.parameters<-data.frame(matrix(NA,N,36))
  all.saved.parameters<-NULL
  all.saved.bpa<-NULL
  
  parameter.names<-c("tips","opt.slope","evol.slope","hl","vy","Ya","theta",
                     "theta.median","theta.mean","theta.025","theta.975",
                     "opt.slope.median","opt.slope.mean","opt.slope.025","opt.slope.975",
                     "evol.slope.median","evol.slope.mean","evol.slope.025","evol.slope.975",
                     "hl.median","hl.mean","hl.025","hl.975","vy.median","vy.mean","vy.025","vy.975",
                     "Ya.median","Ya.mean","Ya.025","Ya.975",
                     "b0.median","b0.mean","b0.025","b0.975",
                     "sym")
  names(saved.parameters)<-parameter.names
  
  #######################
  saved.bpa<-data.frame(matrix(NA,length(hl.list),20))
  parameter.names<-c("tips","opt.slope","evol.slope","hl","vy",
                     "Ya","theta",
                     "theta.bias","theta.imp","theta.inac",
                     "opt.bias","opt.imp","opt.inac",
                     "hl.bias","hl.imp","hl.inac",
                     "Ya.bias","Ya.imp","Ya.inac",
                     "sym")
  names(saved.bpa)<-parameter.names
  
  #######################
  Sxx <- 1
  if(is.null(tips.list)){
      for (j in 1:length(hl.list)){
        A<-log(2)/hl.list[j]
        B<- -A*opt.slope
        Syy<-sqrt(vy*(2*A))

        tip.num<-length(tree$tip.label)
        times <- ape::node.depth.edgelength(tree)
        T.term <- times[1:tip.num]              
        
        #Evolutionary sloperession based on optimal sloperession and A  
        rho<- (1-(1 - exp(-A * T.term))/(A * T.term))
        evol.slope<-(rho*opt.slope)[1]
        
        OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                             B=matrix(B,ncol=1,nrow=1),mPsi=matrix(theta,ncol=1,nrow=1),
                             Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                             Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
        
        stan_data<-stan.adapt.OU.data.setup(OUBMparameters,tree,N)
        
        for (i in 1:N){
          print(paste("Iteration ",i,"hl",hl.list[j],"tips",tip.num,"vy",vy))
          stan_data_sub<-stan_data[[i]]
          stan_data_sub<-stan_data_sub[[1]]
          #fit.random<- stan(file = stan.model,data = stan_data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
          fit<- rstan::sampling(object = stan_model,data = stan_data_sub,chains = 1,iter = 2000,verbose=FALSE,control=list(adapt_delta=0.95))
          par.summary <- summary(fit, pars = c("beta", "beta_e","hl","vy","Ya","b0"), probs = c(0.025,0.5, 0.975))$summary
          par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
          print(par.summary)
          saved.parameters[i,1]<-tip.num
          saved.parameters[i,2]<-opt.slope
          saved.parameters[i,3]<-evol.slope
          saved.parameters[i,4]<-hl.list[j]
          saved.parameters[i,5]<-vy
          saved.parameters[i,6]<-vY0
          saved.parameters[i,7]<-theta
          saved.parameters[i,8:11]<-par.summary[1,] #Theta
          saved.parameters[i,12:15]<-par.summary[2,] #optimal slope
          saved.parameters[i,16:19]<-par.summary[3,] #evolutionary slope
          saved.parameters[i,20:23]<-par.summary[4,] #half-life
          saved.parameters[i,24:27]<-par.summary[5,] #Vy
          saved.parameters[i,28:31]<-par.summary[6,] #Ya
          saved.parameters[i,32:35]<-par.summary[6,] #b0
          saved.parameters[i,36]<-sym
          
                    }
        all.saved.parameters<-rbind(all.saved.parameters,saved.parameters)
        
        #Using median as measurement
        saved.bpa[j,1]<-tip.num
        saved.bpa[j,2]<-opt.slope
        saved.bpa[j,3]<-evol.slope
        saved.bpa[j,4]<-hl.list[j]
        saved.bpa[j,5]<-vy
        saved.bpa[j,6]<-vY0
        saved.bpa[j,7]<-theta
        saved.bpa[j,8]<-mean(saved.parameters[,8]-theta)
        saved.bpa[j,9]<-var(saved.parameters[,8])
        saved.bpa[j,10]<-mean((saved.parameters[,8]-theta)^2)
        saved.bpa[j,11]<-mean(saved.parameters[,12]-opt.slope)
        saved.bpa[j,12]<-var(saved.parameters[,12])
        saved.bpa[j,13]<-mean((saved.parameters[,12]-opt.slope)^2)
        saved.bpa[j,14]<-mean(saved.parameters[,20]-hl.list[j])
        saved.bpa[j,15]<-var(saved.parameters[,20])
        saved.bpa[j,16]<-mean((saved.parameters[,20]-hl.list[j])^2)
        saved.bpa[j,17]<-mean(saved.parameters[,28]-vY0)
        saved.bpa[j,18]<-var(saved.parameters[,28])
        saved.bpa[j,19]<-mean((saved.parameters[,28]-vY0)^2)
        saved.bpa[j,20]<-sym
        
      }
      all.saved.bpa<-rbind(all.saved.bpa,saved.bpa)
  }
  
  if(!is.null(tips.list)){
    for (k in 1:length(tips.list)){
      if(sym==1){sub.tree<-create.tree(tips.list[k],sym=1)}
      else{sub.tree<-subsample.tree(tips.list[k],tree)} #Subsanple existing tree (10K)
      for (j in 1:length(hl.list)){
        A<-log(2)/hl.list[j]
        B<- -A*opt.slope
        Syy<-sqrt(vy*(2*A))
        #Vy=(Syy^2)/(2alpha)
        
        
        tip.num<-length(sub.tree$tip.label)
        times <- ape::node.depth.edgelength(sub.tree)
        T.term <- times[1:tip.num]              
        
        #Evolutionary sloperession based on optimal sloperession and A  
        rho<- (1-(1 - exp(-A * T.term))/(A * T.term))
        evol.slope<-(rho*opt.slope)[1]
        
        OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                             B=matrix(B,ncol=1,nrow=1),mPsi=matrix(theta,ncol=1,nrow=1),
                             Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                             Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
        
        stan_data<-stan.adapt.OU.data.setup(OUBMparameters,sub.tree,N)
        
        for (i in 1:N){
          print(paste("Iteration ",i,"hl",hl.list[j],"tips",tips.list[k],"vy",vy))
          stan_data_sub<-stan_data[[i]]
          stan_data_sub<-stan_data_sub[[1]]
          #fit.random<- stan(file = stan.model,data = stan_data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
          fit<- rstan::sampling(object = stan_model,data = stan_data_sub,chains = 1,iter = 1000,verbose=FALSE,control=list(adapt_delta=0.95))
          par.summary <- summary(fit, pars = c("beta", "beta_e","hl","vy","Ya","b0"), probs = c(0.025,0.5, 0.975))$summary
          par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
          print(par.summary)
          
          saved.parameters[i,1]<-tip.num
          saved.parameters[i,2]<-opt.slope
          saved.parameters[i,3]<-evol.slope
          saved.parameters[i,4]<-hl.list[j]
          saved.parameters[i,5]<-vy
          saved.parameters[i,6]<-vY0
          saved.parameters[i,7]<-theta
          saved.parameters[i,8:11]<-par.summary[1,] #Theta
          saved.parameters[i,12:15]<-par.summary[2,] #optimal slope
          saved.parameters[i,16:19]<-par.summary[3,] #evolutionary slope
          saved.parameters[i,20:23]<-par.summary[4,] #half-life
          saved.parameters[i,24:27]<-par.summary[5,] #Vy
          saved.parameters[i,28:31]<-par.summary[6,] #Ya
          saved.parameters[i,32:35]<-par.summary[6,] #b0
          saved.parameters[i,36]<-sym
          
        }
        all.saved.parameters<-rbind(all.saved.parameters,saved.parameters)
        
        #Using median as measurement
        saved.bpa[j,1]<-tip.num
        saved.bpa[j,2]<-opt.slope
        saved.bpa[j,3]<-evol.slope
        saved.bpa[j,4]<-hl.list[j]
        saved.bpa[j,5]<-vy
        saved.bpa[j,6]<-vY0
        saved.bpa[j,7]<-theta
        saved.bpa[j,8]<-mean(saved.parameters[,8]-theta)
        saved.bpa[j,9]<-var(saved.parameters[,8])
        saved.bpa[j,10]<-mean((saved.parameters[,8]-theta)^2)
        saved.bpa[j,11]<-mean(saved.parameters[,12]-opt.slope)
        saved.bpa[j,12]<-var(saved.parameters[,12])
        saved.bpa[j,13]<-mean((saved.parameters[,12]-opt.slope)^2)
        saved.bpa[j,14]<-mean(saved.parameters[,20]-hl.list[j])
        saved.bpa[j,15]<-var(saved.parameters[,20])
        saved.bpa[j,16]<-mean((saved.parameters[,20]-hl.list[j])^2)
        saved.bpa[j,17]<-mean(saved.parameters[,28]-vY0)
        saved.bpa[j,18]<-var(saved.parameters[,28])
        saved.bpa[j,19]<-mean((saved.parameters[,28]-vY0)^2)
        saved.bpa[j,20]<-sym
        }
      all.saved.bpa<-rbind(all.saved.bpa,saved.bpa)
      }
    }
  return(list(all.saved.parameters,all.saved.bpa))
}


##

stan.direct.OU.inac.hl<-function(tips.list=NULL,tree=NULL,N,direct.slope,hl.list,vy,vY0,vX0,theta,sym=0){
  #Using direct effect OU model
  #Simulations setup to allow for OUBMparameters to be estimated in simulation - SBR1
  
  ############################################
  #setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  #stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOU_v1_8.stan")
  
  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOU_v1_8.stan")
  
  stan_model <- stan_model("blouchOU_v1_8.stan")
  
  #######################
  #################
  saved.parameters<-data.frame(matrix(NA,N,23))
  all.saved.parameters<-NULL
  all.saved.bpa<-NULL
  
  parameter.names<-c("tips","direct.slope","hl","vy","theta",
                     "theta.median","theta.mean","theta.025","theta.975",
                     "direct.slope.median","direct.slope.mean","direct.slope.025","direct.slope.975",
                     "hl.median","hl.mean","hl.025","hl.975",
                     "vy.median","vy.mean","vy.025","vy.975",
                     "sym","expected.slope")
  names(saved.parameters)<-parameter.names
  
  #######################
  saved.bpa<-data.frame(matrix(NA,length(hl.list),15))
  parameter.names<-c("tips","direct.slope","hl","vy","theta",
                     "theta.bias","theta.imp","theta.inac",
                     "direct.bias","direct.imp","direct.inac",
                     "hl.bias","hl.imp","hl.inac","sym")
  names(saved.bpa)<-parameter.names
  
  #######################
  Sxx <- 1
  if(is.null(tips.list)){
    for (j in 1:length(hl.list)){
      A<-log(2)/hl.list[j]
      B<- -A*direct.slope
      Syy<-sqrt(vy*(2*A))

      
      tip.num<-length(tree$tip.label)
      times <- ape::node.depth.edgelength(tree)
      T.term <- times[1:tip.num]              
      #print(T.term)
      #expected.slope<-direct.slope*(sqrt(Syy)+sqrt(Sxx))*(1-exp(-A*T.term[1])/(2*sqrt(Sxx)*T.term[1]))
      expected.slope<-direct.slope*(Syy^2+Sxx)*(1-exp(-A*T.term[1]))/(2*Sxx*T.term[1])
      #Slope = b (sy +sx) (1 - Exp[-at])/(2sx t),
      
      #OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),B=matrix(0,ncol=1,nrow=1),mPsi=matrix(0,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),Syx=matrix(Syx,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1))
      OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                           B=matrix(B,ncol=1,nrow=1),mPsi=matrix(theta,ncol=1,nrow=1),
                           Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                           Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
      
      stan_data<-stan.direct.OU.data.setup(OUBMparameters,tree,N)
      
      for (i in 1:N){
        print(paste("Iteration ",i,"hl",hl.list[j],"tips",tip.num,"vy",vy))
        stan_data_sub<-stan_data[[i]]
        stan_data_sub<-stan_data_sub[[1]]
        #fit.random<- stan(file = stan.model,data = stan_data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
        fit<- rstan::sampling(object = stan_model,data = stan_data_sub,chains = 1,iter = 2000,verbose=FALSE,control=list(adapt_delta=0.95))
        par.summary <- summary(fit, pars = c("beta","hl","vy"), probs = c(0.025,0.5, 0.975))$summary
        par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
        print(par.summary)
        saved.parameters[i,1]<-tip.num
        saved.parameters[i,2]<-direct.slope
        saved.parameters[i,3]<-hl.list[j]
        saved.parameters[i,4]<-vy
        saved.parameters[i,5]<-theta
        saved.parameters[i,6:9]<-par.summary[1,] #Theta
        saved.parameters[i,10:13]<-par.summary[2,] #Direct effect slope
        saved.parameters[i,14:17]<-par.summary[3,] #Half-life
        saved.parameters[i,18:21]<-par.summary[4,] #Vy
        saved.parameters[i,22]<-sym
        saved.parameters[i,23]<-expected.slope
        
      }
      all.saved.parameters<-rbind(all.saved.parameters,saved.parameters)
      
      #Using median as measurement
      saved.bpa[j,1]<-tip.num
      saved.bpa[j,2]<-direct.slope
      saved.bpa[j,3]<-hl.list[j]
      saved.bpa[j,4]<-vy
      saved.bpa[j,5]<-theta
      saved.bpa[j,6]<-mean(saved.parameters[,6]-theta)
      saved.bpa[j,7]<-var(saved.parameters[,6])
      saved.bpa[j,8]<-mean((saved.parameters[,6]-theta)^2)
      saved.bpa[j,9]<-mean(saved.parameters[,10]-direct.slope)
      saved.bpa[j,10]<-var(saved.parameters[,10])
      saved.bpa[j,11]<-mean((saved.parameters[,10]-direct.slope)^2)
      saved.bpa[j,12]<-mean(saved.parameters[,14]-hl.list[j])
      saved.bpa[j,13]<-var(saved.parameters[,14])
      saved.bpa[j,14]<-mean((saved.parameters[,14]-hl.list[j])^2)
      saved.bpa[j,15]<-sym
      
          }
    all.saved.bpa<-rbind(all.saved.bpa,saved.bpa)
  }
  if(!is.null(tips.list)){
    for (k in 1:length(tips.list)){
      if(sym==1){sub.tree<-create.tree(tips.list[k],sym=1)}
      else{sub.tree<-subsample.tree(tips.list[k],tree)} #Subsanple existing tree (10K)
      for (j in 1:length(hl.list)){
        A<-log(2)/hl.list[j]
        B<- -A*direct.slope
        Syy<-sqrt(vy*(2*A))
        tip.num<-length(sub.tree$tip.label)
        times <- ape::node.depth.edgelength(sub.tree)
        T.term <- times[1:tip.num]              
        
        #expected.slope<-direct.slope*(sqrt(Syy)+sqrt(Sxx))*(1-exp(-A*T.term[1])/(2*sqrt(Sxx)*T.term[1]))
        expected.slope<-direct.slope*(Syy^2+Sxx)*(1-exp(-A*T.term[1]))/(2*Sxx*T.term[1])
        
        
        OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                             B=matrix(B,ncol=1,nrow=1),mPsi=matrix(theta,ncol=1,nrow=1),
                             Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                             Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
        
        stan_data<-stan.direct.OU.data.setup(OUBMparameters,sub.tree,N)
        
        for (i in 1:N){
          print(paste("Iteration ",i,"hl",hl.list[j],"tips",tips.list[k],"vy",vy))
          stan_data_sub<-stan_data[[i]]
          stan_data_sub<-stan_data_sub[[1]]
          #fit.random<- stan(file = stan.model,data = stan_data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
          fit<- rstan::sampling(object = stan_model,data = stan_data_sub,chains = 1,iter = 1500,verbose=FALSE,control=list(adapt_delta=0.95))
          par.summary <- summary(fit, pars = c("beta","hl","vy"), probs = c(0.025,0.5, 0.975))$summary
          par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
          print(par.summary)
          saved.parameters[i,1]<-tip.num
          saved.parameters[i,2]<-direct.slope
          saved.parameters[i,3]<-hl.list[j]
          saved.parameters[i,4]<-vy
          saved.parameters[i,5]<-theta
          saved.parameters[i,6:9]<-par.summary[1,] #Theta
          saved.parameters[i,10:13]<-par.summary[2,] #Direct effect slope
          saved.parameters[i,14:17]<-par.summary[3,] #Half-life
          saved.parameters[i,18:21]<-par.summary[4,] #Vy
          saved.parameters[i,22]<-sym
          saved.parameters[i,23]<-expected.slope
          
        }
        all.saved.parameters<-rbind(all.saved.parameters,saved.parameters)
        
        #Using median as measurement
        saved.bpa[j,1]<-tip.num
        saved.bpa[j,2]<-direct.slope
        saved.bpa[j,3]<-hl.list[j]
        saved.bpa[j,4]<-vy
        saved.bpa[j,5]<-theta
        saved.bpa[j,6]<-mean(saved.parameters[,6]-theta)
        saved.bpa[j,7]<-var(saved.parameters[,6])
        saved.bpa[j,8]<-mean((saved.parameters[,6]-theta)^2)
        saved.bpa[j,9]<-mean(saved.parameters[,10]-direct.slope)
        saved.bpa[j,10]<-var(saved.parameters[,10])
        saved.bpa[j,11]<-mean((saved.parameters[,10]-direct.slope)^2)
        saved.bpa[j,12]<-mean(saved.parameters[,14]-hl.list[j])
        saved.bpa[j,13]<-var(saved.parameters[,14])
        saved.bpa[j,14]<-mean((saved.parameters[,14]-hl.list[j])^2)
        saved.bpa[j,15]<-sym
        
      }
      all.saved.bpa<-rbind(all.saved.bpa,saved.bpa)
    }
  }
  return(list(all.saved.parameters,all.saved.bpa))
}

stan.BM.inac.hl<-function(tips.list=NULL,tree=NULL,N,brown.slope,vX0,sigma2y.list,sym=0){
  ############################################
  #setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  #stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchBM_v1.2.stan")
  
  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchBM_v1.2.stan")
  stan_model <- stan_model("blouchBM_v1.2.stan")
  
  #######################
  #################
  saved.parameters<-data.frame(matrix(NA,N,12))
  all.saved.parameters<-NULL
  all.saved.bpa<-NULL
  
  parameter.names<-c("tips","brownian.slope","sigma2y","brown.slope.median","brown.slope.mean","brown.slope.025","brown.slope.0975",
                     "sigma2y.median","sigma2y.mean","sigma2y.025","sigma2y.095","sym")
  names(saved.parameters)<-parameter.names
  
  #######################
  saved.bpa<-data.frame(matrix(NA,length(sigma2y.list),10))
  parameter.names<-c("tips","brown.slope","sigma2y","brown.slope.bias","brown.slope.imp","brown.slope.inac","sigma2y.bias","sigma2y.imp","sigma2y.inac","sym")
  names(saved.bpa)<-parameter.names
  
  #######################
  
  if(is.null(tips.list)){
    for (j in 1:length(sigma2y.list)){
      
      tip.num<-length(tree$tip.label)

      #Brownian slope
      Sxx<-rbind(c(sqrt(sigma2y.list[j]),brown.slope),c(0,1))

      BMparameters<-list(vX0=matrix(vX0,nrow=2,ncol=1),Sxx=matrix(Sxx,nrow=2,ncol=2))
      stan_data<-stan.BM.data.setup(BMparameters,tree,N)
      
      for (i in 1:N){
        print(paste("Iteration ",i,"sigma2y",sigma2y.list[j],"tips",tip.num))
        stan_data_sub<-stan_data[[i]]
        stan_data_sub<-stan_data_sub[[1]]
        #fit.random<- stan(file = stan.model,data = stan_data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
        fit<- rstan::sampling(object = stan_model,data = stan_data_sub,chains = 1,iter = 1500,verbose=FALSE,control=list(adapt_delta=0.95))
        par.summary <- summary(fit, pars = c("beta_evol","sigma2_y"), probs = c(0.025,0.5, 0.975))$summary
        par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
        saved.parameters[i,1]<-tip.num
        saved.parameters[i,2]<-brown.slope
        saved.parameters[i,3]<-sigma2y.list[j]

        saved.parameters[i,4:7]<-par.summary[2,]
        saved.parameters[i,8:11]<-par.summary[3,]
        saved.parameters[j,12]<-sym
      }
      all.saved.parameters<-rbind(all.saved.parameters,saved.parameters)
      
      #Using median as measurement
      saved.bpa[j,1]<-tip.num
      saved.bpa[j,2]<-brown.slope
      saved.bpa[j,3]<-sigma2y.list[j]
      saved.bpa[j,4]<-mean(saved.parameters[,4]-brown.slope)
      saved.bpa[j,5]<-var(saved.parameters[,4])
      saved.bpa[j,6]<-mean((saved.parameters[,4]-brown.slope)^2)
      saved.bpa[j,7]<-mean(saved.parameters[,8]-sigma2y.list[j])
      saved.bpa[j,8]<-var(saved.parameters[,8])
      saved.bpa[j,9]<-mean((saved.parameters[,8]-sigma2y.list[j])^2)
      saved.bpa[j,10]<-sym
    }
    all.saved.bpa<-rbind(all.saved.bpa,saved.bpa)
  }
  if(!is.null(tips.list)){
    for (k in 1:length(tips.list)){
      if(sym==1){sub.tree<-create.tree(tips.list[k],sym=1)}
      else{sub.tree<-subsample.tree(tips.list[k],tree)} #Subsanple existing tree (10K)
      
      for (j in 1:length(sigma2y.list)){

        Sxx<-rbind(c(sqrt(sigma2y.list[j]),brown.slope),c(0,1))
        BMparameters<-list(vX0=matrix(vX0,nrow=2,ncol=1),Sxx=matrix(Sxx,nrow=2,ncol=2))
        stan_data<-stan.BM.data.setup(BMparameters,sub.tree,N)
        
        for (i in 1:N){
          print(paste("Iteration ",i,"sigma2y",sigma2y.list[j],"tips",tips.list[k]))
          stan_data_sub<-stan_data[[i]]
          stan_data_sub<-stan_data_sub[[1]]
          #fit.random<- stan(file = stan.model,data = stan_data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
          fit<- rstan::sampling(object = stan_model,data = stan_data_sub,chains = 1,iter = 2000,verbose=FALSE,control=list(adapt_delta=0.95))
          par.summary <- summary(fit, pars = c("beta_evol","sigma2_y"), probs = c(0.025,0.5, 0.975))$summary
          par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
          saved.parameters[i,1]<-tips.list[k]
          saved.parameters[i,2]<-brown.slope
          saved.parameters[i,3]<-sigma2y.list[j]
          
          saved.parameters[i,4:7]<-par.summary[2,]
          saved.parameters[i,8:11]<-par.summary[3,]
          saved.parameters[i,12]<-sym
          
        }
        all.saved.parameters<-rbind(all.saved.parameters,saved.parameters)
        
        #Using median as measurement
        saved.bpa[j,1]<-tips.list[k]
        saved.bpa[j,2]<-brown.slope
        saved.bpa[j,3]<-sigma2y.list[j]
        saved.bpa[j,4]<-mean(saved.parameters[,4]-brown.slope)
        saved.bpa[j,5]<-var(saved.parameters[,4])
        saved.bpa[j,6]<-mean((saved.parameters[,4]-brown.slope)^2)
        saved.bpa[j,7]<-mean(saved.parameters[,8]-sigma2y.list[j])
        saved.bpa[j,8]<-var(saved.parameters[,8])
        saved.bpa[j,9]<-mean((saved.parameters[,8]-sigma2y.list[j])^2)
        saved.bpa[j,10]<-sym
        
      }
      all.saved.bpa<-rbind(all.saved.bpa,saved.bpa)
    }
  }
  return(list(all.saved.parameters,all.saved.bpa))
}



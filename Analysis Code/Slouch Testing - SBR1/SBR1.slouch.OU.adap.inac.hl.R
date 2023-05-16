#SBR1 OU - Calculating Inacuracy in hl code for slouch
#1. Calculate innacuracies given a list of hl values to explore - 10 reps per half life value
#Calculates bias, precision, and innacuracy in optimal slopes and half-lives

slouch.OU.data.setup<-function(OUBMparameters,tree,n){
  name.response.trait<-c("trait_1")
  names.random.traits<-c("trait_2")
  
  OUBM.trdata.list<-list()
  for(i in 1:n){
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
    OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
    OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
    #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X SBR1
    #OUBM.trdata$dat$trait_1<-OUBM.trdata$dat$trait_1-mean(OUBM.trdata$dat$trait_1) #Mean standardize Y SBR1
    OUBM.trdata.list[[length(OUBM.trdata.list)+1]]<-OUBM.trdata
  }
  return(OUBM.trdata.list)
}

create.tree<-function(n,sym){
  if(sym != 1){
    tree <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = n)
  }
  else{
    #Ultrametric tree - sym = 1, symmetric
    tree <- sim.taxa(1, n=tips, waitsp=function()2,symmetric=TRUE, complete=TRUE)
    tree<-tree[[1]]
  }
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

slouch.adapt.OU.inac.hl<-function(tips.list=NULL,tree=NULL,n,opt.slope,hl.list,vy,vY0,vX0,theta,sym=0){
  #Simulations setup to allow for OUBMparameters to be estimated in simulation - SBR1
  #Can take either a list of tree tips - e.g (32,64,128) and will simulate individual trees for each size and then iterate
  #through the half-life list for them based on n simulations of data, or using one set tree (e.g. tree=something)
  #############################
  #Setup for simulations
  #1% - Added in sigma2y from Hansen et al. 2008 - vy scaled by A
  #A = log(2)/0.01 = 69.31472
  #B = -A * slope = -69.31472 * 0.25 = -17.32868
  library(matrixcalc)
  #######################
  #################
  saved.parameters<-data.frame(matrix(NA,n,12))
  all.saved.parameters<-NULL
  all.saved.bpa<-NULL
  
  parameter.names<-c("tips","opt.slope","evol.slope","hl","vy",
                     "theta","theta.mean","opt.slope.mean","evol.slope.mean","hl.mean","vy.mean","sym")
  names(saved.parameters)<-parameter.names
  
  #######################
  saved.bpa<-data.frame(matrix(NA,length(hl.list),16))
  parameter.names<-c("tips","opt.slope","evol.slope","hl","vy","theta",
                     "theta.bias","theta.imp","theta.inac",
                     "opt.bias","opt.imp","opt.inac",
                     "hl.bias","hl.imp","hl.inac","sym")
  names(saved.bpa)<-parameter.names
  #######################
  Sxx <- 1
  saved.trees<-list()
  saved.data<-list()
  ###############################
  if(is.null(tips.list)){
    for (j in 1:length(hl.list)){
      A<-log(2)/hl.list[j]
      B<- -A*opt.slope
      Syy<-sqrt(vy*(2*A))
      #Syy<-(vy*(2*A))
      
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
      
      
      for (i in 1:n){
        print(paste("Iteration ",i,"hl",hl.list[j],"tips",tip.num,"vy",vy))
        OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
        OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
        OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
        #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X SBR1
        #OUBM.trdata$dat$trait_1<-OUBM.trdata$dat$trait_1-mean(OUBM.trdata$dat$trait_1) #Mean standardize Y SBR1
        
        try(OUBM.slouch<-slouch.fit(phy = OUBM.trdata$phy,
                                    species = OUBM.trdata$phy$tip.label,
                                    response = (OUBM.trdata$dat$trait_1),
                                    random.cov = (OUBM.trdata$dat$trait_2),
                                    mv.response = NULL,
                                    mv.random.cov = NULL,
                                    hl_values = seq(0.00001, 3, length.out = 50),
                                    vy_values = seq(0.00001, 3, length.out = 50),
                                    hillclimb = TRUE,convergence = 150,
                                    lower = c(0.00001, 0.00001)),silent=TRUE)
        
        if(OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 5){
        #if(OUBM.slouch$beta_primary$coefficients_bias_corr[2]> 1.0 || OUBM.slouch$beta_primary$coefficients_bias_corr[2]< -1.0 || OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 1.0){
          saved.trees<-list(saved.trees,tree)
          saved.data<-list(saved.data,OUBM.trdata)
        }
        
        while(OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 5){
          #while(OUBM.slouch$beta_primary$coefficients_bias_corr[2]> 1.0 || OUBM.slouch$beta_primary$coefficients_bias_corr[2]< -1.0 || OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 1.0){
            
          OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
          OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
          OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
          #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X SBR1
          #OUBM.trdata$dat$trait_1<-OUBM.trdata$dat$trait_1-mean(OUBM.trdata$dat$trait_1) #Mean standardize Y SBR1
          
          try(OUBM.slouch<-slouch.fit(phy = OUBM.trdata$phy,
                                      species = OUBM.trdata$phy$tip.label,
                                      response = (OUBM.trdata$dat$trait_1),
                                      random.cov = (OUBM.trdata$dat$trait_2),
                                      mv.response = NULL,
                                      mv.random.cov = NULL,
                                      hl_values = seq(0.00001, 3, length.out = 50),
                                      vy_values = seq(0.00001, 3, length.out = 50),
                                      hillclimb = TRUE,convergence = 150,
                                      lower = c(0.00001, 0.00001)),silent=TRUE)
        }
        
        #print(summary(OUBM.slouch))
        saved.parameters[i,1]<-tip.num
        saved.parameters[i,2]<-opt.slope
        saved.parameters[i,3]<-evol.slope
        saved.parameters[i,4]<-hl.list[j]
        saved.parameters[i,5]<-vy
        saved.parameters[i,6]<-theta
        saved.parameters[i,7]<-OUBM.slouch$beta_primary$coefficients_bias_corr[1]
        saved.parameters[i,8]<-OUBM.slouch$beta_primary$coefficients_bias_corr[2]
        saved.parameters[i,9]<-OUBM.slouch$beta_evolutionary$coefficients[2]
        saved.parameters[i,10]<-OUBM.slouch$evolpar[1]
        saved.parameters[i,11]<-OUBM.slouch$evolpar[2]
        saved.parameters[i,12]<-sym
      }
      all.saved.parameters<-rbind(all.saved.parameters,saved.parameters)
      
      #Using median as measurement
      saved.bpa[j,1]<-tip.num
      saved.bpa[j,2]<-opt.slope
      saved.bpa[j,3]<-evol.slope
      saved.bpa[j,4]<-hl.list[j]
      saved.bpa[j,5]<-vy
      saved.bpa[j,6]<-theta
      
      saved.bpa[j,7]<-mean(saved.parameters[,7]-theta)
      saved.bpa[j,8]<-var(saved.parameters[,7])
      saved.bpa[j,9]<-mean((saved.parameters[,7]-theta)^2)
      
      saved.bpa[j,10]<-mean(saved.parameters[,8]-opt.slope)
      saved.bpa[j,11]<-var(saved.parameters[,8])
      saved.bpa[j,12]<-mean((saved.parameters[,8]-opt.slope)^2)
      
      saved.bpa[j,13]<-mean(saved.parameters[,10]-hl.list[j])
      saved.bpa[j,14]<-var(saved.parameters[,10])
      saved.bpa[j,15]<-mean((saved.parameters[,10]-hl.list[j])^2)
      
      saved.bpa[j,16]<-sym
      
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
        #Syy<-(vy*(2*A))
        
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
      

        for (i in 1:n){
          print(paste("Iteration ",i,"hl",hl.list[j],"tips",tips.list[k],"vy",vy))
          OUBMdata<-simulMVSLOUCHProcPhylTree(sub.tree,OUBMparameters,regimes=NULL) #Simulate data
          OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
          OUBM.trdata <- make.treedata(sub.tree, OUBMdata,name_column="tip.names") #Combine data with 
          #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X SBR1
          #OUBM.trdata$dat$trait_1<-OUBM.trdata$dat$trait_1-mean(OUBM.trdata$dat$trait_1) #Mean standardize Y SBR1
          
          try(OUBM.slouch<-slouch.fit(phy = OUBM.trdata$phy,
                                  species = OUBM.trdata$phy$tip.label,
                                  response = (OUBM.trdata$dat$trait_1),
                                  random.cov = (OUBM.trdata$dat$trait_2),
                                  mv.response = NULL,
                                  mv.random.cov = NULL,
                                  hl_values = seq(0.00001, 3, length.out = 50),
                                  vy_values = seq(0.00001, 3, length.out = 50),
                                  hillclimb = TRUE,convergence = 150,
                                  lower = c(0.00001, 0.00001)),silent=TRUE)
          
          if(OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 5){
            #if(OUBM.slouch$beta_primary$coefficients_bias_corr[2]> 1.0 || OUBM.slouch$beta_primary$coefficients_bias_corr[2]< -1.0 || OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 1.0){
              
            saved.trees<-list(saved.trees,sub.tree)
            saved.data<-list(saved.data,OUBM.trdata)
          }
            
          while(OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 5){
            #while(OUBM.slouch$beta_primary$coefficients_bias_corr[2]> 1.0 || OUBM.slouch$beta_primary$coefficients_bias_corr[2]< -1.0 || OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 1.0){
              OUBMdata<-simulMVSLOUCHProcPhylTree(sub.tree,OUBMparameters,regimes=NULL) #Simulate data
            OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
            OUBM.trdata <- make.treedata(sub.tree, OUBMdata,name_column="tip.names") #Combine data with 
            #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X SBR1
            #OUBM.trdata$dat$trait_1<-OUBM.trdata$dat$trait_1-mean(OUBM.trdata$dat$trait_1) #Mean standardize Y SBR1
            
            try(OUBM.slouch<-slouch.fit(phy = OUBM.trdata$phy,
                                    species = OUBM.trdata$phy$tip.label,
                                    response = (OUBM.trdata$dat$trait_1),
                                    random.cov = (OUBM.trdata$dat$trait_2),
                                    mv.response = NULL,
                                    mv.random.cov = NULL,
                                    hl_values = seq(0.00001, 3, length.out = 50),
                                    vy_values = seq(0.00001, 3, length.out = 50),
                                    hillclimb = TRUE,convergence = 150,
                                    lower = c(0.00001, 0.00001)),silent=TRUE)
          }
          
          #print(summary(OUBM.slouch))
          saved.parameters[i,1]<-tip.num
          saved.parameters[i,2]<-opt.slope
          saved.parameters[i,3]<-evol.slope
          saved.parameters[i,4]<-hl.list[j]
          saved.parameters[i,5]<-vy
          saved.parameters[i,6]<-theta
          saved.parameters[i,7]<-OUBM.slouch$beta_primary$coefficients_bias_corr[1]
          saved.parameters[i,8]<-OUBM.slouch$beta_primary$coefficients_bias_corr[2]
          saved.parameters[i,9]<-OUBM.slouch$beta_evolutionary$coefficients[2]
          saved.parameters[i,10]<-OUBM.slouch$evolpar[1]
          saved.parameters[i,11]<-OUBM.slouch$evolpar[2]
          saved.parameters[i,12]<-sym
        }
        all.saved.parameters<-rbind(all.saved.parameters,saved.parameters)
        
        #Using median as measurement
        saved.bpa[j,1]<-tip.num
        saved.bpa[j,2]<-opt.slope
        saved.bpa[j,3]<-evol.slope
        saved.bpa[j,4]<-hl.list[j]
        saved.bpa[j,5]<-vy
        saved.bpa[j,6]<-theta
        
        saved.bpa[j,7]<-mean(saved.parameters[,7]-theta)
        saved.bpa[j,8]<-var(saved.parameters[,7])
        saved.bpa[j,9]<-mean((saved.parameters[,7]-theta)^2)
        
        saved.bpa[j,10]<-mean(saved.parameters[,8]-opt.slope)
        saved.bpa[j,11]<-var(saved.parameters[,8])
        saved.bpa[j,12]<-mean((saved.parameters[,8]-opt.slope)^2)
        saved.bpa[j,13]<-mean(saved.parameters[,10]-hl.list[j])
        saved.bpa[j,14]<-var(saved.parameters[,10])
        saved.bpa[j,15]<-mean((saved.parameters[,10]-hl.list[j])^2)
        saved.bpa[j,16]<-sym
        
      }
      all.saved.bpa<-rbind(all.saved.bpa,saved.bpa)
  }
  }
  return(list(all.saved.parameters,all.saved.bpa,saved.trees,saved.data))
}


##
slouch.direct.OU.inac.hl<-function(tips.list=NULL,tree=NULL,n,direct.slope,hl.list,vy,vY0,vX0,theta,sym=0){
  #Simulations setup to allow for OUBMparameters to be estimated in simulation - SBR1
  #Can take either a list of tree tips - e.g (32,64,128) and will simulate individual trees for each size and then iterate
  #through the half-life list for them based on n simulations of data, or using one set tree (e.g. tree=something)
  #############################
  #Setup for simulations
  #1% - Added in sigma2y from Hansen et al. 2008 - vy scaled by A
  #A = log(2)/0.01 = 69.31472
  #B = -A * slope = -69.31472 * 0.25 = -17.32868
  library(matrixcalc)
  #######################
  #################
  saved.parameters<-data.frame(matrix(NA,n,11))
  all.saved.parameters<-NULL
  all.saved.bpa<-NULL
  
  parameter.names<-c("tips","direct.slope","hl","vy","theta","theta.mean",
                     "direct.slope.mean","hl.mean","vy.mean","sym","expected.slope")
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
  saved.trees<-list()
  saved.data<-list()
  ###############################
  if(is.null(tips.list)){
    for (j in 1:length(hl.list)){
      A<-log(2)/hl.list[j]
      B<- -A*direct.slope
      Syy<-sqrt(vy*(2*A))
      #Syy<-(vy*(2*A))
      
      tip.num<-length(tree$tip.label)
      times <- ape::node.depth.edgelength(tree)
      T.term <- times[1:tip.num]              
      
      #Evolutionary slope regrssion based on optimal slope regrssion and A  
      rho<- (1-(1 - exp(-A * T.term))/(A * T.term))

      expected.slope<-direct.slope*(Syy^2+Sxx)*(1-exp(-A*T.term[1])/(2*Sxx*T.term[1]))
      #Slope = b (sy +sx) (1 - Exp[-at])/(2sx t) #Estimated based on TH Conversation

      
      OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                           B=matrix(B,ncol=1,nrow=1),mPsi=matrix(theta,ncol=1,nrow=1),
                           Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                           Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
      
      
      for (i in 1:n){
        print(paste("Iteration ",i,"hl",hl.list[j],"tips",tip.num,"vy",vy))
        OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
        OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
        OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
        #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X SBR1
        #OUBM.trdata$dat$trait_1<-OUBM.trdata$dat$trait_1-mean(OUBM.trdata$dat$trait_1) #Mean standardize Y SBR1
        
        try(OUBM.slouch<-slouch.fit(phy = OUBM.trdata$phy,
                                    species = OUBM.trdata$phy$tip.label,
                                    response = (OUBM.trdata$dat$trait_1),
                                    random.cov = (OUBM.trdata$dat$trait_2),
                                    mv.response = NULL,
                                    mv.random.cov = NULL,
                                    hl_values = seq(0.00001, 3, length.out = 50),
                                    vy_values = seq(0.00001, 3, length.out = 50),
                                    hillclimb = TRUE,convergence = 150,
                                    lower = c(0.00001, 0.00001)),silent=TRUE)
        
        if(OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 5){
          #if(OUBM.slouch$beta_primary$coefficients_bias_corr[2]> 1.0 || OUBM.slouch$beta_primary$coefficients_bias_corr[2]< -1.0 || OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 1.0){
            
          saved.trees<-list(saved.trees,tree)
          saved.data<-list(saved.data,OUBM.trdata)
        }
        
        while(OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 5){
          #while(OUBM.slouch$beta_primary$coefficients_bias_corr[2]> 1.0 || OUBM.slouch$beta_primary$coefficients_bias_corr[2]< -1.0 || OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 1.0){
            
                    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
          OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
          OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
          #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X SBR1
          #OUBM.trdata$dat$trait_1<-OUBM.trdata$dat$trait_1-mean(OUBM.trdata$dat$trait_1) #Mean standardize Y SBR1
          
          try(OUBM.slouch<-slouch.fit(phy = OUBM.trdata$phy,
                                      species = OUBM.trdata$phy$tip.label,
                                      response = (OUBM.trdata$dat$trait_1),
                                      random.cov = (OUBM.trdata$dat$trait_2),
                                      mv.response = NULL,
                                      mv.random.cov = NULL,
                                      hl_values = seq(0.00001, 3, length.out = 50),
                                      vy_values = seq(0.00001, 3, length.out = 50),
                                      hillclimb = TRUE,convergence = 150,
                                      lower = c(0.00001, 0.00001)),silent=TRUE)
        }
        
        #print(summary(OUBM.slouch))
        saved.parameters[i,1]<-tip.num
        saved.parameters[i,2]<-direct.slope
        saved.parameters[i,3]<-hl.list[j]
        saved.parameters[i,4]<-vy
        saved.parameters[i,5]<-theta
        saved.parameters[i,6]<-OUBM.slouch$beta_primary$coefficients_bias_corr[1]
        saved.parameters[i,7]<-OUBM.slouch$beta_primary$coefficients_bias_corr[2]
        saved.parameters[i,8]<-OUBM.slouch$evolpar[1]
        saved.parameters[i,9]<-OUBM.slouch$evolpar[2]
        saved.parameters[i,10]<-sym
        saved.parameters[i,11]<-expected.slope
        
        
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
      
      saved.bpa[j,9]<-mean(saved.parameters[,7]-direct.slope)
      saved.bpa[j,10]<-var(saved.parameters[,7])
      saved.bpa[j,11]<-mean((saved.parameters[,7]-direct.slope)^2)
      
      saved.bpa[j,12]<-mean(saved.parameters[,8]-hl.list[j])
      saved.bpa[j,13]<-var(saved.parameters[,8])
      saved.bpa[j,14]<-mean((saved.parameters[,8]-hl.list[j])^2)
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
        
        #Evolutionary sloperession based on optimal sloperession and A  
        rho<- (1-(1 - exp(-A * T.term))/(A * T.term))
        expected.slope<-direct.slope*(Syy^2+Sxx)*(1-exp(-A*T.term[1])/(2*Sxx*T.term[1]))
        #Slope = b (sy +sx) (1 - Exp[-at])/(2sx t) #Estimated based on TH Conversation
        

        OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                             B=matrix(B,ncol=1,nrow=1),mPsi=matrix(theta,ncol=1,nrow=1),
                             Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                             Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
        
        
        for (i in 1:n){
          print(paste("Iteration ",i,"hl",hl.list[j],"tips",tips.list[k],"vy",vy))
          OUBMdata<-simulMVSLOUCHProcPhylTree(sub.tree,OUBMparameters,regimes=NULL) #Simulate data
          OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
          OUBM.trdata <- make.treedata(sub.tree, OUBMdata,name_column="tip.names") #Combine data with 
          #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X SBR1
          #OUBM.trdata$dat$trait_1<-OUBM.trdata$dat$trait_1-mean(OUBM.trdata$dat$trait_1) #Mean standardize Y SBR1
          
          try(OUBM.slouch<-slouch.fit(phy = OUBM.trdata$phy,
                                      species = OUBM.trdata$phy$tip.label,
                                      response = (OUBM.trdata$dat$trait_1),
                                      random.cov = (OUBM.trdata$dat$trait_2),
                                      mv.response = NULL,
                                      mv.random.cov = NULL,
                                      hl_values = seq(0.00001, 3, length.out = 50),
                                      vy_values = seq(0.00001, 3, length.out = 50),
                                      hillclimb = TRUE,convergence = 150,
                                      lower = c(0.00001, 0.00001)),silent=TRUE)
          
          if(OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 5){
#          if(OUBM.slouch$beta_primary$coefficients_bias_corr[2]> 1.0 || OUBM.slouch$beta_primary$coefficients_bias_corr[2]< -1.0 || OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 1.0){

                        saved.trees<-list(saved.trees,sub.tree)
            saved.data<-list(saved.data,OUBM.trdata)
          }
          while(OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 5){
            
          #while(OUBM.slouch$beta_primary$coefficients_bias_corr[2]> 1.0 || OUBM.slouch$beta_primary$coefficients_bias_corr[2]< -1.0 || OUBM.slouch$evolpar[1] > 5 || OUBM.slouch$evolpar[2] > 1.0){
            OUBMdata<-simulMVSLOUCHProcPhylTree(sub.tree,OUBMparameters,regimes=NULL) #Simulate data
            OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
            OUBM.trdata <- make.treedata(sub.tree, OUBMdata,name_column="tip.names") #Combine data with 
            #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X SBR1
            #OUBM.trdata$dat$trait_1<-OUBM.trdata$dat$trait_1-mean(OUBM.trdata$dat$trait_1) #Mean standardize Y SBR1
            
            try(OUBM.slouch<-slouch.fit(phy = OUBM.trdata$phy,
                                        species = OUBM.trdata$phy$tip.label,
                                        response = (OUBM.trdata$dat$trait_1),
                                        random.cov = (OUBM.trdata$dat$trait_2),
                                        mv.response = NULL,
                                        mv.random.cov = NULL,
                                        hl_values = seq(0.00001, 3, length.out = 50),
                                        vy_values = seq(0.00001, 3, length.out = 50),
                                        hillclimb = TRUE,convergence = 150,
                                        lower = c(0.00001, 0.00001)),silent=TRUE)
          }
          
          #print(summary(OUBM.slouch))
          saved.parameters[i,1]<-tip.num
          saved.parameters[i,2]<-direct.slope
          saved.parameters[i,3]<-hl.list[j]
          saved.parameters[i,4]<-vy
          saved.parameters[i,5]<-theta
          saved.parameters[i,6]<-OUBM.slouch$beta_primary$coefficients_bias_corr[1]
          saved.parameters[i,7]<-OUBM.slouch$beta_primary$coefficients_bias_corr[2]
          saved.parameters[i,8]<-OUBM.slouch$evolpar[1]
          saved.parameters[i,9]<-OUBM.slouch$evolpar[2]
          saved.parameters[i,10]<-sym
          saved.parameters[i,11]<-expected.slope
          
          
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
        
        saved.bpa[j,9]<-mean(saved.parameters[,7]-direct.slope)
        saved.bpa[j,10]<-var(saved.parameters[,7])
        saved.bpa[j,11]<-mean((saved.parameters[,7]-direct.slope)^2)
        
        saved.bpa[j,12]<-mean(saved.parameters[,8]-hl.list[j])
        saved.bpa[j,13]<-var(saved.parameters[,8])
        saved.bpa[j,14]<-mean((saved.parameters[,8]-hl.list[j])^2)
        saved.bpa[j,15]<-sym
        
      }
      all.saved.bpa<-rbind(all.saved.bpa,saved.bpa)
    }
  }
  return(list(all.saved.parameters,all.saved.bpa,saved.trees,saved.data))
}


slouch.trends.BM.inac.hl<-function(tips.list=NULL,tree=NULL,n,brown.slope,sigma2y.list,vY0,vX0,sym=0){
  #Simulations setup to allow for BMparameters to be estimated in simulation - SBR1
  #Can take either a list of tree tips - e.g (32,64,128) and will simulate individual trees for each size and then iterate
  #through the half-life list for them based on n simulations of data, or using one set tree (e.g. tree=something)
  #############################
  #Setup for simulations
  #1% - Added in sigma2y from Hansen et al. 2008 - vy scaled by A
  #A = log(2)/0.01 = 69.31472
  #B = -A * slope = -69.31472 * 0.25 = -17.32868
  library(matrixcalc)
  #######################
  #################
  saved.parameters<-data.frame(matrix(NA,n,6))
  all.saved.parameters<-NULL
  all.saved.bpa<-NULL
  
  parameter.names<-c("tips","brown.slope","sigma2y","brown.slope","sigma2y.mean","sym")
  names(saved.parameters)<-parameter.names
  
  #######################
  saved.bpa<-data.frame(matrix(NA,length(hl.list),10))
  parameter.names<-c("tips","brown.slope","sigma2y","brown.bias","brown.imp","brown.inac","sigma2y.bias","sigma2y.imp","sigma2y.inac","sym")
  names(saved.bpa)<-parameter.names
  #######################
  Sxx <- 1
  saved.trees<-list()
  saved.data<-list()
  ###############################
  if(is.null(tips.list)){
    for (j in 1:length(sigma2y.list)){
      
      #Brownian slope
      Sxx<-rbind(c(sqrt(sigma2y.list[j]),brown.slope),c(0,1))
      
      BMparameters<-list(vX0=matrix(vX0,nrow=2,ncol=1),Sxx=matrix(Sxx,nrow=2,ncol=2))
      
      
      for (i in 1:n){
        print(paste("Iteration ",i,"hl",hl.list[j],"tips",tip.num,"vy",vy))
        BMdata<-simulMVSLOUCHProcPhylTree(tree,BMparameters,regimes=NULL) #Simulate data
        BMdata<-data.frame(BMdata,"tip.names"=row.names(BMdata))
        BM.trdata <- make.treedata(tree, BMdata,name_column="tip.names") #Combine data with 
        #BM.trdata$dat$trait_2<-BM.trdata$dat$trait_2-mean(BM.trdata$dat$trait_2) #Mean standardize X SBR1
        #BM.trdata$dat$trait_1<-BM.trdata$dat$trait_1-mean(BM.trdata$dat$trait_1) #Mean standardize Y SBR1
        
        try(BM.slouch<-slouch.fit(phy = BM.trdata$phy,
                                  species = BM.trdata$phy$tip.label,
                                  response = (BM.trdata$dat$trait_1),
                                  random.cov = (BM.trdata$dat$trait_2),
                                  mv.response = NULL,
                                  mv.random.cov = NULL,
                                  hl_values = seq(0.00001, 3, length.out = 50),
                                  vy_values = seq(0.00001, 3, length.out = 50),
                                  hillclimb = TRUE,convergence = 150,
                                  lower = c(0.00001, 0.00001)),silent=TRUE)
        
        if(BM.slouch$evolpar[1] > 5 || BM.slouch$evolpar[2] > 5){
         # if(BM.slouch$beta_primary$coefficients_bias_corr[2]> 1.0 || BM.slouch$beta_primary$coefficients_bias_corr[2]< -1.0 || BM.slouch$evolpar[1] > 5 || BM.slouch$evolpar[2] > 1.0){
            
                    saved.trees<-list(saved.trees,tree)
          saved.data<-list(saved.data,BM.trdata)
        }
        
        while(BM.slouch$evolpar[1] > 5 || BM.slouch$evolpar[2] > 5){
         # while(BM.slouch$beta_primary$coefficients_bias_corr[2]> 1.0 || BM.slouch$beta_primary$coefficients_bias_corr[2]< -1.0 || BM.slouch$evolpar[1] > 5 || BM.slouch$evolpar[2] > 1.0){
            
                    BMdata<-simulMVSLOUCHProcPhylTree(tree,BMparameters,regimes=NULL) #Simulate data
          BMdata<-data.frame(BMdata,"tip.names"=row.names(BMdata))
          BM.trdata <- make.treedata(tree, BMdata,name_column="tip.names") #Combine data with 
          BM.trdata$dat$trait_2<-BM.trdata$dat$trait_2-mean(BM.trdata$dat$trait_2) #Mean standardize X SBR1
          BM.trdata$dat$trait_1<-BM.trdata$dat$trait_1-mean(BM.trdata$dat$trait_1) #Mean standardize Y SBR1
          
          try(BM.slouch<-slouch.fit(phy = BM.trdata$phy,
                                    species = BM.trdata$phy$tip.label,
                                    response = (BM.trdata$dat$trait_1),
                                    random.cov = (BM.trdata$dat$trait_2),
                                    mv.response = NULL,
                                    mv.random.cov = NULL,
                                    hl_values = seq(0.00001, 3, length.out = 50),
                                    vy_values = seq(0.00001, 3, length.out = 50),
                                    hillclimb = TRUE,convergence = 150,
                                    lower = c(0.00001, 0.00001)),silent=TRUE)
        }
        
        #print(summary(BM.slouch))
        saved.parameters[i,1]<-tip.num
        saved.parameters[i,2]<-brown.slope
        saved.parameters[i,3]<-sigma2y.list[j]
        saved.parameters[i,4]<-BM.slouch$beta_primary$coefficients_bias_corr[2]
        saved.parameters[i,5]<-BM.slouch$evolpar[1]
        saved.parameters[i,6]<-sym
        
      }
      all.saved.parameters<-rbind(all.saved.parameters,saved.parameters)
      
      #Using median as measurement
      saved.bpa[j,1]<-tip.num
      saved.bpa[j,2]<-brown.slope
      saved.bpa[j,3]<-sigma2y.list[j]
      saved.bpa[j,4]<-mean(saved.parameters[,4]-brown.slope)
      saved.bpa[j,5]<-var(saved.parameters[,4])
      saved.bpa[j,6]<-mean((saved.parameters[,4]-brown.slope)^2)
      saved.bpa[j,7]<-mean(saved.parameters[,5]-sigma2y.list[j])
      saved.bpa[j,8]<-var(saved.parameters[,5])
      saved.bpa[j,9]<-mean((saved.parameters[,5]-sigma2y.list[j])^2)
      saved.bpa[j,10]<-sym
      
    }
    all.saved.bpa<-rbind(all.saved.bpa,saved.bpa)
  }
  if(!is.null(tips.list)){
    for (k in 1:length(tips.list)){
      if(sym==1){sub.tree<-create.tree(tips.list[k],sym=1)}
      else{sub.tree<-subsample.tree(tips.list[k],tree)} #Subsanple existing tree (10K)
      
      for (j in 1:length(hl.list)){

        #Brownian slope
        Sxx<-rbind(c(sqrt(sigma2y.list[j]),brown.slope),c(0,1))
        
        BMparameters<-list(vX0=matrix(vX0,nrow=2,ncol=1),Sxx=matrix(Sxx,nrow=2,ncol=2))
        
        
        for (i in 1:n){
          print(paste("Iteration ",i,"hl",hl.list[j],"tips",tips.list[k],"vy",vy))
          BMdata<-simulMVSLOUCHProcPhylTree(sub.tree,BMparameters,regimes=NULL) #Simulate data
          BMdata<-data.frame(BMdata,"tip.names"=row.names(BMdata))
          BM.trdata <- make.treedata(sub.tree, BMdata,name_column="tip.names") #Combine data with 
          #BM.trdata$dat$trait_2<-BM.trdata$dat$trait_2-mean(BM.trdata$dat$trait_2) #Mean standardize X SBR1
          #BM.trdata$dat$trait_1<-BM.trdata$dat$trait_1-mean(BM.trdata$dat$trait_1) #Mean standardize Y SBR1
          
          try(BM.slouch<-slouch.fit(phy = BM.trdata$phy,
                                    species = BM.trdata$phy$tip.label,
                                    response = (BM.trdata$dat$trait_1),
                                    random.cov = (BM.trdata$dat$trait_2),
                                    mv.response = NULL,
                                    mv.random.cov = NULL,
                                    hl_values = seq(0.00001, 3, length.out = 50),
                                    vy_values = seq(0.00001, 3, length.out = 50),
                                    hillclimb = TRUE,convergence = 150,
                                    lower = c(0.00001, 0.00001)),silent=TRUE)
          
          if(BM.slouch$evolpar[1] > 5 || BM.slouch$evolpar[2] > 5){
            #if(BM.slouch$beta_primary$coefficients_bias_corr[2]> 1.0 || BM.slouch$beta_primary$coefficients_bias_corr[2]< -1.0 || BM.slouch$evolpar[1] > 5 || BM.slouch$evolpar[2] > 1.0){
              
                        saved.trees<-list(saved.trees,sub.tree)
            saved.data<-list(saved.data,BM.trdata)
          }
          
          #while(BM.slouch$beta_primary$coefficients_bias_corr[2]> 1.0 || BM.slouch$beta_primary$coefficients_bias_corr[2]< -1.0 || BM.slouch$evolpar[1] > 5 || BM.slouch$evolpar[2] > 1.0){
          while(BM.slouch$evolpar[1] > 5 || BM.slouch$evolpar[2] > 5){
              
                        BMdata<-simulMVSLOUCHProcPhylTree(sub.tree,BMparameters,regimes=NULL) #Simulate data
            BMdata<-data.frame(BMdata,"tip.names"=row.names(BMdata))
            BM.trdata <- make.treedata(sub.tree, BMdata,name_column="tip.names") #Combine data with 
            BM.trdata$dat$trait_2<-BM.trdata$dat$trait_2-mean(BM.trdata$dat$trait_2) #Mean standardize X SBR1
            BM.trdata$dat$trait_1<-BM.trdata$dat$trait_1-mean(BM.trdata$dat$trait_1) #Mean standardize Y SBR1
            
            try(BM.slouch<-slouch.fit(phy = BM.trdata$phy,
                                      species = BM.trdata$phy$tip.label,
                                      response = (BM.trdata$dat$trait_1),
                                      random.cov = (BM.trdata$dat$trait_2),
                                      mv.response = NULL,
                                      mv.random.cov = NULL,
                                      hl_values = seq(0.00001, 3, length.out = 50),
                                      vy_values = seq(0.00001, 3, length.out = 50),
                                      hillclimb = TRUE,convergence = 150,
                                      lower = c(0.00001, 0.00001)),silent=TRUE)
          }
          
          #print(summary(BM.slouch))
          saved.parameters[i,1]<-tip.num
          saved.parameters[i,2]<-brown.slope
          saved.parameters[i,3]<-sigma2y.list[j]
          saved.parameters[i,4]<-BM.slouch$beta_primary$coefficients_bias_corr[2]
          saved.parameters[i,5]<-BM.slouch$evolpar[1]
          saved.parameters[i,6]<-sym
          
        }
        all.saved.parameters<-rbind(all.saved.parameters,saved.parameters)
        
        #Using median as measurement
        saved.bpa[j,1]<-tip.num
        saved.bpa[j,2]<-opt.slope
        saved.bpa[j,3]<-sigma2y.list[j]
        saved.bpa[j,4]<-mean(saved.parameters[,4]-brown.slope)
        saved.bpa[j,5]<-var(saved.parameters[,4])
        saved.bpa[j,6]<-mean((saved.parameters[,4]-brown.slope)^2)
        saved.bpa[j,7]<-mean(saved.parameters[,5]-sigma2y.list[j])
        saved.bpa[j,8]<-var(saved.parameters[,5])
        saved.bpa[j,9]<-mean((saved.parameters[,5]-sigma2y.list[j])^2)
        saved.bpa[j,10]<-sym
      }
      all.saved.bpa<-rbind(all.saved.bpa,saved.bpa)
    }
  }
  return(list(all.saved.parameters,all.saved.bpa,saved.trees,saved.data))
}

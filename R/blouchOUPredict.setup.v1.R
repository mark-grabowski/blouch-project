# Save this file as `R/blouchOUPredict.setup.v1.R`

#' This is the R setup file for blouchOUPredict_v1.stan
#' It includes the function sigma.X.estimate, taken from the R Pacakge Slouch (Kopperud et al. 2020).
#'
#' @export
#' @param trdata Data formatted by make.treedata function from the R Pacakge treeplyr - only extant species
#' @param trdata.fos Data formatted by make.treedata function from the R Pacakge treeplyr - extant and fossil
#' @param names.traits Vector of trait names
#' @param classical Numeric denoting whether classical regression (=1) or inverse regression (=0) should be used in prediction
#' @return An object of class list for use in blouchOU_v1.stan
#'


sigma.X.estimate <-
  function (phy, ta, predictor, mv.predictor) {
    predictor <- matrix(predictor, nrow = length(phy$tip.label))

    N <- length(phy$tip.label)
    v <- ta # Time from root to most recent ancestor
    w <- matrix(data = 1, nrow = N, ncol = 1)
    me <- diag(mv.predictor)
    dat <- predictor
    beta1 <- solve(t(w) %*% solve(v) %*% w) %*% (t(w) %*% solve(v) %*% dat)
    e <- dat - c(beta1)
    sigma_squared <- as.numeric((t(e) %*% solve(v) %*% e) / (N-1))
    repeat{
      beta1 <- solve(t(w) %*% solve(v + me/sigma_squared) %*% w) %*% (t(w) %*% solve(v + me/sigma_squared) %*% dat)
      e <- dat - c(beta1)
      sigma_squared1 <- (t(e) %*% solve(v + me/sigma_squared) %*% e) / (N-1)
      if (abs(as.numeric(sigma_squared1) - sigma_squared) <= 0.0000001 * sigma_squared){
        break
      }
      sigma_squared <- as.numeric(sigma_squared1)
    }
    return(list(mean = as.numeric(beta1),
                sigma_squared = as.numeric(sigma_squared)))
  }


blouchOUPredict.setup.v1<-function(trdata,trdata.fos,names.traits,classical){
  #Trait names are in the order trait 1 is adapting to trait 2, data should already be logged
  #Names.traits = trait 1, trait 2, mv trait 1, mv trait 2
  fos.index<-which(trdata.fos$dat$Status=="Extinct")
  print(paste("Fossil Species #",fos.index))
  extant.index<-which(trdata.fos$dat$Status=="Extant")
  #Define tree parameters
  phy<-trdata$phy
  #max(branching.times(phy))
  phy$edge.length<-phy$edge.length/max(branching.times(phy)) ## rescale tree to height 1

  n<-length(phy$tip.label)
  mrca1 <- ape::mrca(phy)
  times <- ape::node.depth.edgelength(phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(phy$tip.label, phy$tip.label))
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  if(is.na(names.traits[3])!=TRUE){
    mv.response<-as.vector(pull(trdata$dat,names.traits[3]))}
  else{
    mv.response<-as.vector(rep(0,n))}

  if(is.na(names.traits[4])!=TRUE){
    mv.pred<-matrix(pull(trdata$dat,names.traits[4]),nrow=n,ncol=1)}
  else{
    mv.pred<-matrix(0,nrow=n,ncol=1)}

  test<-sigma.X.estimate(phy,ta, predictor = (pull(trdata$dat,names.traits[2])), mv.predictor = mv.pred)

  brownian_mean<-test[1]
  sigma_squared_x<-test[2]

  coef<-lm(Y~X,data=data.frame(Y=(pull(trdata$dat,names.traits[1])),X=(pull(trdata$dat,names.traits[2]))))$coef
  ols.intercept<-coef[1]
  ols.slope<-coef[2]
  ###############################################################################################################
  #Fossil prediction code
  phy.fos<-trdata.fos$phy
  max(branching.times(phy.fos))
  phy.fos$edge.length<-phy.fos$edge.length/max(branching.times(phy.fos)) ## rescale tree to height 1

  n.fos<-length(phy.fos$tip.label)
  mrca1 <- ape::mrca(phy.fos)
  times <- ape::node.depth.edgelength(phy.fos)
  ta.fos <- matrix(times[mrca1], nrow=n.fos, dimnames = list(phy.fos$tip.label, phy.fos$tip.label))
  T.term.fos <- times[1:n.fos]
  tia.fos <- times[1:n.fos] - ta.fos
  tja.fos <- t(tia.fos)
  tij.fos <- tja.fos + tia.fos

  if(is.na(names.traits[3])!=TRUE){
    mv.response.fos<-as.vector(pull(trdata.fos$dat,names.traits[3]))}
  else{
    mv.response.fos<-as.vector(rep(0,n.fos))}
  if(is.na(names.traits[4])!=TRUE){
    mv.pred.fos<-matrix(pull(trdata.fos$dat,names.traits[4]),nrow=n.fos,ncol=1)}
  else{
    mv.pred.fos<-matrix(0,nrow=n.fos,ncol=1)}

  #Using same as extant species for predictors - ignoring fossil data for now
  brownian_mean_fos<-test[1]
  sigma_squared_x_fos<-test[2]



  ####################

  stan_constraint_data<-list(N = n,
                             Z = 1,
                             Y = as.vector((pull(trdata$dat,names.traits[1]))),
                             mv_response = mv.response,
                             random_cov = matrix(0,nrow = n,ncol=1),
                             mv_random_cov = matrix(0,nrow=n,ncol = 1),
                             direct_cov = as.matrix((pull(trdata$dat,names.traits[2]))),
                             mv_direct_cov = mv.pred,
                             ta = as.matrix(ta),
                             T_term = as.vector(T.term),
                             tia = as.matrix(tia),
                             tja = as.matrix(tja),
                             tij = as.matrix(tij),
                             brownian_mean = brownian_mean,
                             sigma_squared_x = sigma_squared_x,
                             ols_intercept = ols.intercept,
                             ols_slope = ols.slope,

##################################################
  #Extant mean predictions
                            extant_count = 1:n,
  ###############################################
  #Fossil species means
                            N_fos = n.fos,
                            N_fos_only = n.fos-n,
                            Y_fos_means = as.vector((pull(trdata.fos$dat,names.traits[1]))),
                            mv_response_fos = mv.response.fos,
                            random_cov_fos = matrix(0,nrow = n.fos,ncol=1),
                            mv_random_cov_fos = matrix(0,nrow=n.fos,ncol = 1),
                            direct_cov_fos = as.matrix((pull(trdata.fos$dat,names.traits[2]))),
                            mv_direct_cov_fos = mv.pred.fos,
                            ta_fos = as.matrix(ta.fos),
                            T_term_fos = as.vector(T.term.fos),
                            tia_fos = as.matrix(tia.fos),
                            tja_fos = as.matrix(tja.fos),
                            tij_fos = as.matrix(tij.fos),
                            brownian_mean_fos=brownian_mean_fos,
                            sigma_squared_x_fos=sigma_squared_x_fos,
                            fos_index = as.array(fos.index),
                            extant_index = as.vector(extant.index),
                            classical = classical
  )


  #############################################################################################
  #Setting up data - random covariates
  stan_adaptive_data<-list(N = n,
                           Z = 1,
                           Y = as.vector((pull(trdata$dat,names.traits[1]))),
                           mv_response = mv.response,
                           random_cov = as.matrix((pull(trdata$dat,names.traits[2]))),
                           mv_random_cov = mv.pred,
                           direct_cov =  matrix(0,nrow = n,ncol=1),
                           mv_direct_cov = matrix(0,nrow = n,ncol=1),
                           ta = as.matrix(ta),
                           T_term = as.vector(T.term),
                           tia = as.matrix(tia),
                           tja = as.matrix(tja),
                           tij = as.matrix(tij),
                           brownian_mean = brownian_mean,
                           sigma_squared_x = sigma_squared_x,
                           ols_intercept = ols.intercept,
                           ols_slope = ols.slope,

                           ##################################################
                           #Extant mean predictions
                           extant_count = 1:n,
                           ###############################################
                           #Fossil species means
                           N_fos = n.fos,
                           N_fos_only = n.fos-n,
                           Y_fos_means = as.vector((pull(trdata.fos$dat,names.traits[1]))),
                           mv_response_fos = mv.response.fos,

                           random_cov_fos = as.matrix((pull(trdata.fos$dat,names.traits[2]))),
                           mv_random_cov_fos = mv.pred.fos,
                           direct_cov_fos = matrix(0,nrow = n.fos,ncol=1),
                           mv_direct_cov_fos = matrix(0,nrow = n.fos,ncol=1),

                           ta_fos = as.matrix(ta.fos),
                           T_term_fos = as.vector(T.term.fos),
                           tia_fos = as.matrix(tia.fos),
                           tja_fos = as.matrix(tja.fos),
                           tij_fos = as.matrix(tij.fos),
                           brownian_mean_fos=brownian_mean_fos,
                           sigma_squared_x_fos=sigma_squared_x_fos,
                           fos_index = as.array(fos.index),
                           extant_index = as.vector(extant.index),
                           classical = classical
                           )




  return(list(stan_constraint_data,stan_adaptive_data))
}



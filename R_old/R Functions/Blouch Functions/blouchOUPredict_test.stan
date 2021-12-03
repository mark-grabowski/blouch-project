functions {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Blouch v1 - 11.23.2020
//This version only recreates Slouch output, with no fossil prediction or within species data
//To be renamed Blouch_OU for ms
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code - Adaptive
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    matrix design_matrix(int N, int evol, real a, vector T_term, matrix direct_cov, matrix random_cov, int Z){
    matrix[N,Z] rho;
    matrix[N,Z] X;

    rho = to_matrix(1 - (1 - exp(-a * T_term))./(a * T_term)); //For OU model

    if(sum(random_cov)==0){
      X = direct_cov;
      }
    else if(sum(direct_cov)==0){
      X = random_cov .* rho;
    }else{
      X = append_col(direct_cov, random_cov .* rho);
      }
    //print("Design Matrix",X);
    return(X);
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code - Evolutionary
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    matrix design_matrix_evol(int N, int evol, real a, vector T_term, matrix direct_cov, matrix random_cov, int Z){
    matrix[N,Z] rho;
    matrix[N,Z+1] X;
   
    rho=to_matrix(rep_vector(1,N)); //For Evolutionary regression
    
    if(sum(random_cov)==0){
      X = append_col(rep_vector(1,N), direct_cov);
      }
    else if(sum(direct_cov)==0){
      X = append_col(rep_vector(1,N), random_cov .* rho);

    }else{
      X = append_col(rep_vector(1,N),direct_cov);
      X = append_col(X,random_cov .* rho);
      }
    //print("Design Matrix",X);
    return(X);
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Vt function
///////////////////////////////////////////////////////////////////////////////////////////////////////////
  matrix varcov_model(int N, matrix tij, matrix tja, matrix ta, matrix random_cov,int Z, real sigma2_y, real a, matrix x0, matrix sigma_squared_x, vector beta1, vector T_term){
  
  vector[N] sigma2s;
  matrix[N,N] ti;
  matrix[N,N] term0;
  matrix[N,N] term1;
  matrix[N,N] term2;
  real s1;
  vector[Z] beta1sq;
  matrix[N,N] Vt;

  
  for (i in 1:Z)
    beta1sq[i] = beta1[i]^2;
    if (a < 1e-14){
    if(sum(random_cov) == 0){ //Direct covariates and small a
        Vt = sigma2_y * ta;
    }else{
      s1 = sum(sigma_squared_x * beta1sq);
      Vt = sigma2_y * ta + s1 * ta .* ((ta .* ta) ./ 12 + tja .* (tja') ./ 4);
    }
  }
  else{ //Random + Direct covariates, or just random
    if(sum(random_cov) != 0){
      s1 = sum(sigma_squared_x * beta1sq);
      ti = rep_matrix(T_term,N);

      term0 = ((s1 + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) .* exp(-a * tij);
      term1 = (1 - exp(-a * ti)) ./ (a * ti); 
      term2 = exp(-a * tja) .* (1 - exp(-a * ti)) ./ (a * ti);

      Vt = term0 + s1 * (ta .* term1 .* (term1') - ((1 - exp(-a * ta)) ./ a) .* (term2 + (term2')));
      //print("Random Only");
    }
    else{ //Direct covariates only
      //print("Direct Only")
      Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij));
    }
  }    
  
  
  return(Vt);
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//V_me function
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrix varcov_measurement(int N, int Z, matrix ta, matrix direct_cov, matrix mv_direct_cov, matrix mv_random_cov, matrix sigma_squared_x, vector beta){
  //matrix[N,N] Vxt[Z];
  matrix[N,N] Vxtr;
  matrix[N,N] Vxtd;
  matrix[N,N] Vxt;
  matrix[N,N] Vx;
  matrix[N,N] Vu_given_x;
  matrix[N,N] beta2_Vu_given_x;
  matrix[N,N] Vu;
  matrix[N,N] Vur;
  matrix[N,N] Vud;
  matrix[N,Z] sigma_squared_x_rep;
  matrix[N,Z] P;
  Vur = rep_matrix(0,N,N);
  Vud = rep_matrix(0,N,N);
  Vxtd = rep_matrix(0,N,N);
  Vxtr = rep_matrix(0,N,N);

  if(sum(mv_random_cov) != 0){
    Vur = diag_matrix(to_vector(mv_random_cov));
    Vxtr = ta * sigma_squared_x[1,1];
    }
    
  if(sum(mv_direct_cov) != 0){
    Vud = diag_matrix(to_vector(mv_direct_cov[,1])); //Puts 
    Vxtd = diag_matrix(rep_vector(variance(direct_cov[,1]),N))-Vud;
    }

  Vxt = Vxtd+Vxtr;
  Vu = Vud+Vur;
  Vx = Vxt + Vu;

//Vu_given_x
  Vu_given_x = Vu - Vu * inverse(Vx) * Vu;

//beta2_Vu_given_x
  beta2_Vu_given_x = Vu_given_x * square(beta[1]);
  
  //if(is_nan(beta2_Vu_given_x[1,1])==1){
//    beta2_Vu_given_x=rep_matrix(0,N,N);
//  }
    
  return(beta2_Vu_given_x);
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
data {
  //Extant data
  int N; //species number
  int Z; //number of traits
  vector[N] Y; //y variable
  vector[N] mv_response;
  matrix[N,Z] direct_cov;
  matrix[N,Z] mv_direct_cov;
  matrix[N,Z] random_cov;
  matrix[N,Z] mv_random_cov;
  matrix[N,N] ta; //The following calculated in R based on the phylogeny
  vector[N] T_term; 
  matrix[N,N] tia;
  matrix[N,N] tja;
  matrix[N,N] tij;
  matrix[1,Z] brownian_mean;
  matrix[1,Z] sigma_squared_x;
  real ols_intercept;
  real ols_slope;

 //Extant prediction
  int extant_count[N];

  //Fossil species data
  int N_fos; //Total species number + fossils
  int N_fos_only; //Fossil species number
  vector[N_fos] Y_fos_means; //y variable
  vector[N_fos] mv_response_fos; //y variable
  matrix[N_fos,Z] direct_cov_fos;
  matrix[N_fos,Z] mv_direct_cov_fos;
  matrix[N_fos,Z] random_cov_fos;
  matrix[N_fos,Z] mv_random_cov_fos;

  matrix[N_fos,N_fos] ta_fos; //The following calculated in R based on the phylogeny
  vector[N_fos] T_term_fos; 
  matrix[N_fos,N_fos] tia_fos;
  matrix[N_fos,N_fos] tja_fos;
  matrix[N_fos,N_fos] tij_fos;
  matrix[1,Z] brownian_mean_fos;
  matrix[1,Z] sigma_squared_x_fos;
  int fos_index[N_fos_only];
  int extant_index[N];
  
  int classical;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
parameters {
  real <lower = 0> a;
  real <lower = 0, upper = variance(Y)*4> sigma2_y; //Added to limit the variance based on Kjetil's suggestion
  //real <lower = 0> sigma2_y; //Added to limit the variance based on Kjetil's suggestion
  real alpha; //OU alpha
  vector[Z] beta; //OU beta

  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
transformed parameters {
  matrix[N,Z+1] X_evol;
  vector[Z+1] beta_evol;
  matrix[N,N] V_ev;
  matrix[N,N] V_me_ev;
  matrix[N,N] Vt_ev;

  Vt_ev = varcov_model(N,  tij,  tja,  ta, random_cov, Z,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term);
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    //print("sum!=0");
    V_me_ev = varcov_measurement(N, Z, ta, direct_cov, mv_direct_cov, mv_random_cov, sigma_squared_x, beta);
    V_ev = Vt_ev + V_me_ev + diag_matrix(mv_response);

    }
  else{
    //print("sum=0");
    V_me_ev = rep_matrix(0,N,N);
    V_ev = Vt_ev;
    }
  //Calculate evolutionary regression slope
  X_evol = design_matrix_evol(N,  1,  a,  T_term, direct_cov,random_cov, Z);
  //X_evol = random_cov;
  //ev.beta1.var <- pseudoinverse(t(X0)%*%V.inverse%*%X0)

  //ev.beta1 <- ev.beta1.var%*%(t(X0)%*%V.inverse%*%Y)
  beta_evol = inverse(X_evol'*inverse(V_ev)*X_evol)*(X_evol'*inverse(V_ev)*Y); //Hansen et al. 2008
  //beta_evol = inverse(X_evol'*inverse(V_final)*X_evol)*X_evol'*inverse(V_final)*Y; //Hansen et al. 2008
  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
model {
//Declare variables  
  //real a;
  vector[N] mu;
  matrix[N,N] V;
  matrix[N,N] Vt;
  matrix[N,N] V_me;
  matrix[N,Z] X;
  matrix[N, N] L_V;
  
//Priors
//Priors
  a ~ lognormal(1.0,1.0); //Possibly better for half-lives, but slow for smaller half-lives
  //sigma2_y ~ exponential(0.1);
  alpha ~ normal(ols_intercept,0.5); //Changed from 0.5 for MBM ms
  beta ~ normal(ols_slope,0.4); //Changed to 0.5 for MBM ms, Cervidae - latter because OLS slope is like 6, so needs more
//////////////////////////////////////////////////////////////////////////////////////////////////////
  //a = log(2)/hl;
////////////////////////////////////////////////////////////////////////////////////////////////////////
//Regression - either constraint for direct cov or adaptive for random cov

//Set up X matrix
  X = design_matrix( N,  0,  a,  T_term, direct_cov,random_cov, Z);
  
//Set up V matix
  Vt = varcov_model(N,  tij,  tja,  ta,  random_cov, Z,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term);
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    V_me = varcov_measurement(N, Z, ta, direct_cov, mv_direct_cov,mv_random_cov, sigma_squared_x, beta);
  }
  else{
    V_me = rep_matrix(0,N,N);
  }
  //print(Vt,V_me,diag_matrix(mv_response));
  V = Vt + V_me + diag_matrix(mv_response);
  L_V = cholesky_decompose(V);
  
//OU with random covariates
  mu = X*beta+alpha;
  Y ~ multi_normal_cholesky(mu , L_V);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////


generated quantities {
  real <lower = 0> vy;
  real <lower = 0> hl;
  //real <lower = 0> a;
  vector[N] pred_mean;
  real grand_mean;
  real sst;
  real sse;
  real r_squared;
  
  /////////////////////////////////////////////////////////////////////
  //R2 and evolutionary regression
  matrix[N,Z] X_opt;
  matrix[N,N] V_final;
  matrix[N,N] V_me_final;
  matrix[N,N] Vt_final;

  ///////////////////////////////////////////////////////////
  //Fossil Predictions
  matrix[N_fos,N_fos] V_fos;
  matrix[N_fos,N_fos] Vt_fos;
  matrix[N_fos_only,N_fos_only] L_V_fos;
  matrix[N_fos_only, N_fos_only] fos_V;
  matrix[N_fos-N_fos_only, N_fos-N_fos_only] extant_V;
  matrix[N_fos_only, N_fos-N_fos_only] fos_extant_V;
  matrix[N_fos_only,N_fos_only] fos_V_new;

  vector[N_fos] X_mu_evol;
  vector[N_fos_only] X_mu_fos;
  vector[N] X_mu_extant;
  matrix[N_fos_only,Z] X_fos_new;
  vector[N_fos_only] X_pred_fos_means;

  vector[N_fos] Y_mu_evol;
  vector[N_fos_only] Y_mu_fos;
  vector[N] Y_mu_extant;
  vector[N_fos_only] Y_fos_new;
  vector[N_fos_only] Y_pred_fos_means;

  //Extant predictions
  int entant_rest_index[N-1];
  matrix[N-1,N-1] extant_drop_V;
  real drop_V;
  matrix[1,N-1] drop_extant_V;
  
  real X_mu_extant_drop;
  vector[N-1] X_mu_extant_rest;
  real X_extant_new;
  real extant_V_new;
  vector[N] X_pred_extant_means;
  
  real Y_mu_extant_drop;
  vector[N-1] Y_mu_extant_rest;
  real Y_extant_new;
  vector[N] Y_pred_extant_means;
  
  real e_MSE[N];
  real MSE;
  real RMSE;
  real e_mu_MSE[N];
  real MSE_mu;
  real RMSE_mu;
//////////

//////////
  hl = log(2)/a;
  //a = log(2)/hl;
  vy = sigma2_y/(2*a);
  
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Calculate r2 based on constraint or adaptive regression
  X_opt = design_matrix(N,  0,  a,  T_term, direct_cov,random_cov, Z);
  Vt_final = varcov_model(N,  tij,  tja,  ta,  random_cov, Z,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term);
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    V_me_final = varcov_measurement(N, Z, ta, direct_cov, mv_direct_cov, mv_random_cov, sigma_squared_x, beta);
    }
  else{
    V_me_final = rep_matrix(0,N,N);
    }

  V_final = Vt_final + V_me_final + diag_matrix(mv_response);
  
  pred_mean = (X_opt*beta+alpha);
  grand_mean = ((rep_vector(1,N))' * inverse(V_final) * Y) / sum(inverse(V_final));
  sst = ((Y - grand_mean)' * inverse(V_final) * (Y - grand_mean));
  sse = ((Y - pred_mean)' * inverse(V_final) * (Y - pred_mean));
  r_squared = (sst - sse) / sst;
  
  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Mean Extant prediction - always use evolutionary regression for prediction
  if(classical==1){ //For classical vs. inverse regression - predicting X vs. predicting Y
   if(sum(direct_cov) != 0){//X_mu = mean body mass predictions
    X_mu_extant = inverse(beta[1] * inverse(V_final) * beta[1]) * beta[1] * inverse(V_final) * (Y - alpha); //Non-evolutionary for regression
    //X_mu_extant = inverse(beta[1] * inverse(V_final) * beta[1]') * beta[1] * inverse(V_final) * (Y - alpha); //Non-evolutionary for regression

  }else{
    //X_mu_extant = inverse(beta_evol[1] * inverse(V_final) * beta_evol[1]') * beta_evol[1] * inverse(V_final) * (Y - alpha); //Evolutionary regression for optimal analysis
    X_mu_extant = inverse(beta_evol[2] * inverse(V_final) * beta_evol[2]) * beta_evol[2] * inverse(V_final) * (Y - alpha); //Evolutionary regression for optimal analysis

  }
  }
  else{
    Y_mu_extant = direct_cov*beta+alpha;
  }

  for(i in 1:N){
      if(i<N && i>1){
        entant_rest_index = append_array(head(extant_count,i-1),tail(extant_count,N-i));
      }
      if(i==N){
        entant_rest_index = head(extant_count,i-1); //Beginning to i
      }
      if(i==1){
        entant_rest_index = tail(extant_count,N-1); //Beginning to i
      }
      extant_drop_V = V_final[entant_rest_index,entant_rest_index];
      drop_V = V_final[i,i];
      drop_extant_V = to_matrix(V_final[i,entant_rest_index]); //Variance of missing species based on covariance matrix
   
   if(classical==1){ //For classical vs. inverse regression - predicting X vs. predicting Y
      X_mu_extant_drop = X_mu_extant[i];
      X_mu_extant_rest = X_mu_extant[entant_rest_index];
      X_extant_new = (X_mu_extant_drop + drop_extant_V * inverse(extant_drop_V) * (direct_cov[entant_rest_index] - to_matrix(X_mu_extant_rest)))[1,1];
      extant_V_new = (drop_V - drop_extant_V * inverse(extant_drop_V) * (drop_extant_V'))[1,1]; //From Martins and Hansen, 1997
      X_pred_extant_means[i] = normal_rng(X_extant_new,sqrt(extant_V_new));
      
      e_MSE[i] = ((X_pred_extant_means[i])-(direct_cov[i,1]))^2;
      e_mu_MSE[i] = ((X_mu_extant[i])-(direct_cov[i,1]))^2;
  
      
      Y_mu_extant_drop = 0;
      Y_mu_extant_rest = rep_vector(0,N-1);
      Y_extant_new = 0;
      Y_pred_extant_means = rep_vector(0,N);
   }
   else{
      Y_mu_extant_drop = Y_mu_extant[i];
      Y_mu_extant_rest = Y_mu_extant[entant_rest_index];
      Y_extant_new = (Y_mu_extant_drop + drop_extant_V * inverse(extant_drop_V) * (Y[entant_rest_index] - Y_mu_extant_rest))[1];
      extant_V_new = (drop_V - drop_extant_V * inverse(extant_drop_V) * (drop_extant_V'))[1,1]; //From Martins and Hansen, 1997
      Y_pred_extant_means[i] = normal_rng(Y_extant_new,sqrt(extant_V_new));
      
      e_MSE[i] = ((Y_pred_extant_means[i])-(Y[i]))^2;
      e_mu_MSE[i] = ((Y_mu_extant[i])-(Y[i]))^2;
      
      X_mu_extant_drop = 0;
      X_mu_extant_rest = rep_vector(0,N-1);
      X_extant_new = 0;
      X_pred_extant_means = rep_vector(0,N);
   }
  }

  //MSE = sum(e_MSE)/(N-Z-1);
  MSE = sum(e_MSE)/(N-2);
  RMSE=sqrt(MSE);
  //MSE_mu=sum(e_mu_MSE)/(N-Z-1);
  MSE_mu=sum(e_mu_MSE)/(N-2);
  RMSE_mu =sqrt(MSE_mu);

//
//Fossil species mean prediction code
  //X_fos = design_matrix(N_fos,  evol,  a,  T_term_fos, direct_cov_fos,random_cov_fos, Z);
  Vt_fos = varcov_model(N_fos,  tij_fos,  tja_fos,  ta_fos,  random_cov, Z,  sigma2_y,  a, brownian_mean_fos,  sigma_squared_x_fos,  beta,  T_term_fos);
  //V_me_fos = varcov_measurement(N_fos, Z, ta_fos,direct_cov_fos, mv_direct_cov_fos, mv_random_cov_fos, sigma_squared_x, beta);
  if(classical == 1){ //For classical vs. inverse regression - predicting X vs. predicting Y
    V_fos = Vt_fos + diag_matrix(mv_response_fos);
  }else{
    V_fos = Vt_fos;}

  if(classical == 1){ //For classical vs. inverse regression - predicting X vs. predicting Y
  if(sum(direct_cov) != 0){//X_mu = mean body mass predictions
    //X_mu_evol = inverse(beta[1] * inverse(V_fos) * beta[1]') * beta[1] * inverse(V_fos) * (Y_fos_means - alpha); //Non-evolutionary for regression
    X_mu_evol = inverse(beta[1] * inverse(V_fos) * beta[1]) * beta[1] * inverse(V_fos) * (Y_fos_means - alpha); //Non-evolutionary for regression

    X_mu_fos = X_mu_evol[fos_index];
    X_mu_extant = X_mu_evol[extant_index]; //Mean of Y/body suze for non missing species (including missing species?)
    Y_mu_evol = rep_vector(0,N_fos);
    Y_mu_fos = rep_vector(0,N_fos_only);
    Y_mu_extant = rep_vector(0,N); //Mean of Y/body suze for non missing species (including missing species?)

  }else{
    //X_mu_evol = inverse(beta_evol[1] * inverse(V_fos) * beta_evol[1]') * beta_evol[1] * inverse(V_fos) * (Y_fos_means - alpha); //Evolutionary regression for optimal analysis
    X_mu_evol = inverse(beta_evol[2] * inverse(V_fos) * beta_evol[2]) * beta_evol[2] * inverse(V_fos) * (Y_fos_means - alpha); //Evolutionary regression for optimal analysis

    X_mu_fos = X_mu_evol[fos_index];
    X_mu_extant = X_mu_evol[extant_index]; //Mean of Y/body suze for non missing species (including missing species?)
    Y_mu_evol = rep_vector(0,N_fos);
    Y_mu_fos = rep_vector(0,N_fos_only);
    Y_mu_extant = rep_vector(0,N); //Mean of Y/body suze for non missing species (including missing species?)

    }
  }
  else{
    X_mu_evol = rep_vector(0,N_fos);
    X_mu_fos  = rep_vector(0,N_fos_only);
    X_mu_extant = rep_vector(0,N);
    X_pred_fos_means = rep_vector(0,N_fos_only);
    X_fos_new = rep_matrix(0,N_fos_only,Z);
    
    Y_mu_evol = direct_cov_fos*beta+alpha;
    Y_mu_fos = Y_mu_evol[fos_index];
    Y_mu_extant = Y_mu_evol[extant_index]; //Mean of Y/body suze for non missing species (including missing species?)
    //print(Y_mu_extant);
  }
  fos_V = V_fos[fos_index,fos_index]; //Variance of missing species based on covariance matrix
  extant_V = V_fos[extant_index,extant_index]; //Variance of non-missing species
  fos_extant_V = V_fos[fos_index,extant_index]; //Covariance of missing with non-missing species
  
    if(classical == 1){ //For classical vs. inverse regression - predicting X vs. predicting Y
      X_fos_new = to_matrix(X_mu_fos) + fos_extant_V * inverse(extant_V) * (direct_cov - to_matrix(X_mu_extant));
      Y_fos_new = rep_vector(0,N_fos_only);
    }
    else{//Inverse regression
      Y_fos_new = (Y_mu_fos) + fos_extant_V * inverse(extant_V) * (Y - (Y_mu_extant));
      X_fos_new = rep_matrix(0,N_fos_only,Z);
    }
  
  if(N_fos_only == 1){
    fos_V_new = fos_V - fos_extant_V * inverse(extant_V) * (fos_extant_V'); //From Martins and Hansen, 1997
}
  if(N_fos_only > 1){
    fos_V_new = fos_V - fos_extant_V * inverse(extant_V) * (fos_extant_V'); //From Martins and Hansen, 1997
  }
  
  L_V_fos = cholesky_decompose(fos_V_new);
  if(classical == 1){ //For classical vs. inverse regression - predicting X vs. predicting Y
    X_pred_fos_means = multi_normal_cholesky_rng(to_vector(X_fos_new),L_V_fos);
    Y_pred_fos_means = rep_vector(0,N_fos_only);}
  else{
    Y_pred_fos_means = multi_normal_cholesky_rng(to_vector(Y_fos_new),L_V_fos);
    X_pred_fos_means = rep_vector(0,N_fos_only);}
}

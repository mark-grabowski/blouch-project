functions {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Blouch v1.8 SBC Version - 12.07.2022
//This code is based on Blouch v1.8, but includes code to do:
//Simulation Based Calibration - half-life, vy, Ya, B0, B1
//This code estimates evolutionary parameters without prediction
//Adding in code to do multivariate analyses (multiple xs) for SB R1
//Allow priors to be set on optimal slope as well as evolutionary slope
//Allow priors on Ya and b0 along with k
//Works with blouchOU.setup.mv_v1_3.R function
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code - Adaptive
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    matrix design_matrix(int N,real a, vector T_term, matrix direct_cov, matrix random_cov, int Z, int Z_direct, int Z_random){
    matrix[N,Z+1] X;
    matrix[N,Z_direct+1] X_part;
    vector[N] rho;
    matrix[N,Z_random] rhos;
    
    rho = (1 - (1 - exp(-a * T_term))./(a * T_term)); //For OU model
    rhos = rep_matrix(rho,Z_random);

    if(num_elements(random_cov)==0){
      X = append_col(rep_vector(1,N), direct_cov);
      }
    else if(num_elements(direct_cov)==0){
      X = append_col(rep_vector(1,N), random_cov .* rhos);
    }else if((num_elements(direct_cov)!=0)&&(num_elements(random_cov)!=0)){
      X_part = append_col(rep_vector(1,N),direct_cov);
      X = append_col(X_part,random_cov .* rhos);}
    return(X);
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code - Evolutionary
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    matrix design_matrix_evol(int N,real a, vector T_term, matrix direct_cov, matrix random_cov, int Z, int Z_direct, int Z_random){
    matrix[N,Z_direct+1] X_part;
    matrix[N,Z+1] X;
    
    if(num_elements(random_cov)==0){
      X = append_col(rep_vector(1,N), direct_cov);
      }
    else if(num_elements(direct_cov)==0){
      X = append_col(rep_vector(1,N), random_cov);
    }else if((num_elements(direct_cov)!=0)&&(num_elements(random_cov)!=0)){
      X_part = append_col(rep_vector(1,N),direct_cov);
      X = append_col(X_part,random_cov);}
    return(X);
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Vt function - Biological variance based on phylogeny and rate of adaptation - Based on Hansen et al. (2008)
///////////////////////////////////////////////////////////////////////////////////////////////////////////
  matrix varcov_model(int N, matrix tij, matrix tja, matrix ta, matrix random_cov,real sigma2_y, real a, matrix x0, matrix sigma_squared_x, vector beta, vector T_term, int Z, int Z_direct, int Z_random){
  
  matrix[N,N] ti;
  matrix[N,N] term0;
  matrix[N,N] term1;
  matrix[N,N] term2;
  real s1;
  matrix[N,N] Vt;
  vector[Z_random] ones;
  
  ones = rep_vector(1,Z_random);
  
  //Random + Direct covariates, or just random
  if(num_elements(random_cov) != 0){
      if(Z_random>1){
        s1 = beta[2+Z_direct:1+Z_direct+Z_random]' * sigma_squared_x * ones;
        }else{
            s1 = sigma_squared_x[1,1] * (beta[1+Z_direct+Z_random]^2);
            }
      ti = rep_matrix(T_term,N);

      term0 = ((s1 + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) .* exp(-a * tij);
      term1 = (1 - exp(-a * ti)) ./ (a * ti); 
      term2 = exp(-a * tja) .* (1 - exp(-a * ti)) ./ (a * ti);

      Vt = term0 + s1 * (ta .* term1 .* (term1') - ((1 - exp(-a * ta)) ./ a) .* (term2 + (term2'))); //From Hansen et al. (2008)
    }
    else{ //Direct covariates only
      Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)
    }
  
  return(Vt);
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//V_me function Measurement variance in the predictors - Based on Hansen and Bartoszek (2012), Equation 10 and associated text
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrix varcov_measurement(int N, matrix ta, matrix direct_cov, matrix mv_direct_cov, matrix mv_random_cov, matrix sigma_squared_x, 
vector beta, int Z, int Z_direct, int Z_random, real a, vector T_term){
  
  matrix[N,Z_direct] direct_rho2;
  matrix[N,Z_random] random_rho2;
  matrix[N,Z_direct+Z_random] rho2s;
  
  matrix[N,N] Vur[Z_random];
  matrix[N,N] Vxtr[Z_random];
  matrix[N,N] Vud[Z_direct];
  matrix[N,N] Vxtd[Z_direct];
  
  matrix[N,N] Vxt[Z];
  matrix[N,N] Vu[Z];
  matrix[N,N] Vx[Z];
  
  vector[N] Vu_given_x[Z];
  vector[N] beta2_Vu_given_x[Z];
  vector[N] beta2_Vu_given_x_sum;

  Vur = rep_array(rep_matrix(0,N,N),Z_random);
  Vxtr = rep_array(rep_matrix(0,N,N),Z_random);

  Vud = rep_array(rep_matrix(0,N,N),Z_direct);
  Vxtd = rep_array(rep_matrix(0,N,N),Z_direct);

  direct_rho2 = rep_matrix(1,N,Z_direct);
  random_rho2 = rep_matrix((1 - (1 - exp(-a * T_term))./(a * T_term))^2,Z_random); //For OU model
  
  rho2s = append_col(direct_rho2,random_rho2);

  if(num_elements(mv_random_cov) != 0){
    for (i in 1:Z_random){
      Vur[i] = diag_matrix(mv_random_cov[,i]);
      Vxtr[i] = ta * sigma_squared_x[1,i];
      }
    }
  if(num_elements(mv_direct_cov) != 0){
    for (i in 1:Z_direct){
      Vud[i] = diag_matrix(mv_direct_cov[,i]);
      Vxtd[i] = diag_matrix(rep_vector(variance(direct_cov[,i]),N))-Vud[i];
      }
  }
  
  Vxt = append_array(Vxtd,Vxtr);
  Vu = append_array(Vud,Vur);
  
  for (i in 1:Z){ //Adding variance matrices together for individual traits 
    Vx[i] = Vxt[i] + Vu[i];
    Vu_given_x[i] = diagonal(Vu[i] - Vu[i] * inverse(Vx[i]) * Vu[i]);
    beta2_Vu_given_x[i] = to_vector(Vu_given_x[i] * square(beta[i+1]) .* rho2s[,i]); //Changed to element-wise, appears to be same as R
    }  

 
  
  for (i in 1:N){
        beta2_Vu_given_x_sum[i] = sum(beta2_Vu_given_x[,i]);
      }

  return(diag_matrix(beta2_Vu_given_x_sum));
  }

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Data Block
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
data {
  //Extant data
  int N; //species number
  int Z; //total number of traits
  int Z_direct; //number of direct traits
  int Z_random; //number of adaptive traits
  int Z_me_direct;
  int Z_me_random;
  
  vector[N] Y; //y variable
  vector[N] mv_response;
  matrix[N,Z_direct] direct_cov;
  matrix[N,Z_me_direct] mv_direct_cov;
  matrix[N,Z_random] random_cov;
  matrix[N,Z_me_random] mv_random_cov;
  
  matrix[N,N] ta; //The following calculated in R based on the phylogeny
  vector[N] T_term; 
  matrix[N,N] tia;
  matrix[N,N] tja;
  matrix[N,N] tij;
  matrix[1,Z_random] brownian_mean;
  matrix[Z_random,Z_random] sigma_squared_x;
  real ols_intercept;
  vector[Z] ols_slope;
  real Ya_prior;
}
/////////////////////////////////////////////////////
transformed data{
  real <lower = 0> hl_sim = lognormal_rng(-1.4,1.1);
  real <lower = 0> vy_sim = fabs(cauchy_rng(0,0.01));
  //real vy_sim = uniform_rng(0,variance(Y)*4);
  //real <lower = 0, upper = variance(Y)*4> sigma2_y; //Added to limit the variance based on Kjetil's suggestion

  vector[Z+1] beta_sim; //OU beta
  real beta1_sim; //OU beta
  real beta2_sim; //OU beta

  vector[Z+1] beta_e_sim;
  real Ya_sim = normal_rng(Ya_prior,0.2);
  real b0_sim = normal_rng(ols_intercept,0.025);
  real a_sim = log(2)/hl_sim;
  real sigma2_y_sim = vy_sim*(2*a_sim);
  matrix[N,N] V_sim;
  matrix[N,N] Vt_sim;
  matrix[N,N] V_me_sim;

  vector[N] rho_sim = (1 - (1 - exp(-a_sim * T_term))./(a_sim * T_term)); //For OU model

  for (i in 1:Z){
    beta_sim[i+1] = normal_rng(ols_slope[i], 0.075); //Slope originally 0.4
  }
  
  real k_sim  = exp(-a_sim * T_term[1])*Ya_sim+(1-exp(-a_sim * T_term[1]))*b0_sim+(1-exp(-a_sim*T_term[1])-rho_sim[1])*beta_sim[2]*brownian_mean[1,1];
  beta_sim[1] = normal_rng(k_sim, 0.025); //Slope originally 0.4
  beta1_sim = beta_sim[1];
  beta2_sim = beta_sim[2];
  
//Set up V matix
  Vt_sim = varcov_model(N,  tij,  tja,  ta, random_cov, sigma2_y_sim,  a_sim, brownian_mean,  sigma_squared_x,  beta_sim,  T_term, Z, Z_direct, Z_random);
  if((num_elements(mv_direct_cov)!= 0) || (num_elements(mv_random_cov)!= 0)){
    V_me_sim = varcov_measurement(N, ta, direct_cov, mv_direct_cov,mv_random_cov, sigma_squared_x, beta_sim, Z, Z_direct, Z_random, a_sim, T_term);
    if(num_elements(mv_response)!= 0){
      V_sim = Vt_sim + V_me_sim + diag_matrix(mv_response); //Hansen & Bartoszek (2012) - Eq 10
    }
    else{
        V_sim = Vt_sim + V_me_sim; //Hansen & Bartoszek (2012) - Eq 10
      }
  }else{
    if(num_elements(mv_response)!= 0){
       V_sim = Vt_sim + diag_matrix(mv_response);
       V_me_sim = rep_matrix(0,N,N);
    }else{
      V_sim = Vt_sim;
      V_me_sim = rep_matrix(0,N,N);
    }
  }
  
  beta_e_sim[1] = beta_sim[1];
  for(i in 1:Z_random){
    beta_e_sim[Z_random+i] = beta_sim[1+Z_direct+i]* rho_sim[1]; //Only random covariates get evolutionary regressions
    }  


  matrix[N,N] L_V_sim = cholesky_decompose(V_sim);
  matrix[N,Z+1] X_sim = design_matrix_evol(N, a_sim,  T_term, direct_cov,random_cov,  Z,  Z_direct,  Z_random);
  vector[N] Y_sim = multi_normal_cholesky_rng(X_sim*beta_e_sim, L_V_sim);
  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Parameter Block
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
parameters {
  real <lower=0> hl;
  //real <lower=0,upper=variance(Y)*4> vy;
  real <lower = 0> vy; //To use with non-uniform priors
  vector[Z+1] beta; //OU beta
  vector[Z_random] beta_e; //OU beta
  real Ya;
  real b0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Transformed Parameter Block
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
transformed parameters {

  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Model Block
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
model {
  //Declare variables
  matrix[N,N] V;
  matrix[N,N] Vt;
  matrix[N,N] V_me;
  vector[N] mu;
  matrix[N,Z+1] X;
  matrix[N, N] L_V;
  real a;
  real sigma2_y;
  real k;
  vector[N] rho;
  //Priors
  hl ~ lognormal(-1.4,1.1); //Tree length = 1, for original Cervidae tree = 17 Ma - simulations
  vy ~ cauchy(0,0.01);
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  a = log(2)/hl;
  sigma2_y = vy*(2*a);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  //rho = (1 - (1 - exp(-a * T_term[1]))./(a * T_term[1])); //For OU model
  if(num_elements(random_cov)!=0){
    rho = (1 - (1 - exp(-a * T_term))./(a * T_term)); //For OU model
    k  = exp(-a * T_term[1])*Ya+(1-exp(-a * T_term[1]))*b0+(1-exp(-a*T_term[1])-rho[1])*beta[2]*brownian_mean[1,1];
    //z = exp(-a * T_term)*Ya+(1-exp(-a * T_term))*b0+(1-exp(-a*T_term)-rho)*beta[2+Z_direct:Z_direct+Z_random+1]'*brownian_mean[1,1:Z_random]'; //MV Case
    //print(z);
  }

  else{
    k = b0;
    rho = rep_vector(1,N);
    }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  Ya ~ normal(Ya_prior,0.2);
  //Ya ~ normal(Ya_prior,sqrt(vy));

  b0 ~ normal(ols_intercept,0.025);

  if(num_elements(random_cov)!=0){
    beta[1] ~ normal(k, 0.025); //Slope originally 0.4
    beta[2:Z+1] ~ normal(ols_slope, 0.075); //Slope originally 0.4
    //beta[1] ~ normal(ols_intercept, 0.4); //Slope originally 0.4
    //beta[2:Z+1] ~ normal(ols_slope, 0.4); //Slope originally 0.4

  }
  else{
    beta[1] ~ normal(ols_intercept, 0.025); //Slope originally 0.4
    beta[2:Z+1] ~ normal(ols_slope, 0.075); //Slope originally 0.4
    }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
//Regression - either constraint for direct cov or adaptive for random cov
//Set up X matrix
  X = design_matrix(N, a,  T_term, direct_cov,random_cov, Z, Z_direct, Z_random);
//Set up V matix
  Vt = varcov_model(N,  tij,  tja,  ta, random_cov, sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term, Z, Z_direct, Z_random);
  if((num_elements(mv_direct_cov)!= 0) || (num_elements(mv_random_cov)!= 0)){
    V_me = varcov_measurement(N, ta, direct_cov, mv_direct_cov,mv_random_cov, sigma_squared_x, beta, Z, Z_direct, Z_random, a, T_term);
    if(num_elements(mv_response)!= 0){
      V = Vt + V_me + diag_matrix(mv_response); //Hansen & Bartoszek (2012) - Eq 10
    }
    else{
        V = Vt + V_me; //Hansen & Bartoszek (2012) - Eq 10
      }
  }else{
    if(num_elements(mv_response)!= 0){
       V = Vt + diag_matrix(mv_response);
       V_me = rep_matrix(0,N,N);
    }else{
      V = Vt;
      V_me = rep_matrix(0,N,N);
    }
  }
  
  L_V = cholesky_decompose(V);
  mu = X*beta;
  Y_sim ~ multi_normal_cholesky(mu , L_V);
  
  for(i in 1:Z_random){
    beta_e[i] ~ normal(beta[1+Z_direct+i]* rho,0.1); //Only random covariates get evolutionary regressions
    }  
  }


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//Generated Quantities Block
////////////////////////////////////////////////////////////////////////////////////////////////////////////
generated quantities {
  real sigma2_y;
  real a;

  vector[6] pars_;
  //////////
  a = log(2)/hl;
  sigma2_y = vy*(2*a);

  pars_[1]=hl_sim;
  pars_[2]=vy_sim;
  pars_[3]=Ya_sim;
  pars_[4]=b0_sim;
  pars_[5]=beta1_sim;
  pars_[6]=beta2_sim;
  
  array[6] int <lower=0,upper=1> ranks_ = {hl<hl_sim,vy<vy_sim,Ya<Ya_sim,b0<b0_sim,beta[1]<beta1_sim,beta[2]<beta2_sim};
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////


}

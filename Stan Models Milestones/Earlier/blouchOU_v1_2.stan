functions {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Blouch v1.1 - 07.20.2022
//This code estimates evolutionary parameters without prediction
//Adding in code to do multivariate analyses (multiple xs) for SB R1
//Works with blouchOU.setup.mv.R function
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code - Adaptive
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    matrix design_matrix(int N,real a, vector T_term, matrix direct_cov, matrix random_cov, int Z, int Z_direct, int Z_random){
    matrix[N,Z] X;
    vector[N] rho;
    matrix[N,Z_random] rhos;
    
    rho = (1 - (1 - exp(-a * T_term))./(a * T_term)); //For OU model
    rhos = rep_matrix(rho,Z_random);

    if(sum(random_cov)==0){
      X = direct_cov;
      }
    else if(sum(direct_cov)==0){
      X = (random_cov .* rhos);
    }else{
      X = append_col(direct_cov, (random_cov .* rhos));
      }
    return(X);
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code - Evolutionary
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    matrix design_matrix_evol(int N,real a, vector T_term, matrix direct_cov, matrix random_cov, int Z, int Z_direct, int Z_random){
    matrix[N,Z_direct+1] X_part;
    matrix[N,Z+1] X;
    
    if(sum(random_cov)==0){
      X = append_col(rep_vector(1,N), direct_cov);
      }
    else if(sum(direct_cov)==0){
      X = append_col(rep_vector(1,N), random_cov);
    }else{
      X_part = append_col(rep_vector(1,N),direct_cov);
      X = append_col(X_part,random_cov);
      }
    return(X);
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Vt function - Biological variance based on phylogeny and rate of adaptation - Based on Hansen et al. (2008)
///////////////////////////////////////////////////////////////////////////////////////////////////////////
  matrix varcov_model(int N, matrix tij, matrix tja, matrix ta, matrix random_cov,real sigma2_y, real a, matrix x0, matrix sigma_squared_x, vector beta, vector T_term, int Z, int Z_direct, int Z_random){
  
  vector[N] sigma2s;
  matrix[N,N] ti;
  matrix[N,N] term0;
  matrix[N,N] term1;
  matrix[N,N] term2;
  real s1;
  vector[Z_random] betasq;
  matrix[N,N] Vt;

  
  for (i in 1:Z_random){
    betasq[i] = beta[i+Z_direct]^2;
    }

  //Random + Direct covariates, or just random
  if(sum(random_cov) != 0){
      s1 = sum(sigma_squared_x * betasq);
      ti = rep_matrix(T_term,N);

      term0 = ((s1 + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) .* exp(-a * tij);
      term1 = (1 - exp(-a * ti)) ./ (a * ti); 
      term2 = exp(-a * tja) .* (1 - exp(-a * ti)) ./ (a * ti);
      //print((1 - exp(-a * ta))[1,3],(a * ta)[1,3]);
    
      Vt = term0 + s1 * (ta .* term1 .* (term1') - ((1 - exp(-a * ta)) ./ a) .* (term2 + (term2')));
      //Vt = term0 + s1 * ta .* (term1 .* (term1') - ((1 - exp(-a * ta)) ./ (a * ta)) .* (term2 + (term2')));    
        
      //print("Random Only");
    }
    else{ //Direct covariates only
      //print("Direct Only")
      Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij));
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

  if(sum(mv_random_cov) != 0){
    for (i in 1:Z_random){
      Vur[i] = diag_matrix(mv_random_cov[,i]);
      Vxtr[i] = ta * sigma_squared_x[1,i];
      }
    }
  if(sum(mv_direct_cov) != 0){
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
    //print("Rep");
    //print(i);
    //print(length(square(beta[i+n_regimes])));
    //print(dims(Vu_given_x[i] * square(beta[i+n_regimes]) .* rho2s[,i])); //Changed to element-wise, appears to be same as R
    beta2_Vu_given_x[i] = to_vector(Vu_given_x[i] * square(beta[i]) .* rho2s[,i]); //Changed to element-wise, appears to be same as R
    }  

 
  
  for (i in 1:N){
        //print(beta2_Vu_given_x[,i]);
        beta2_Vu_given_x_sum[i] = sum(beta2_Vu_given_x[,i]);
      }
  //print(beta2_Vu_given_x_sum);

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
  vector[N] Y; //y variable
  vector[N] mv_response;
  matrix[N,Z_direct] direct_cov;
  matrix[N,Z_direct] mv_direct_cov;
  matrix[N,Z_random] random_cov;
  matrix[N,Z_random] mv_random_cov;
  matrix[N,N] ta; //The following calculated in R based on the phylogeny
  vector[N] T_term; 
  matrix[N,N] tia;
  matrix[N,N] tja;
  matrix[N,N] tij;
  matrix[1,Z_random] brownian_mean;
  matrix[1,Z_random] sigma_squared_x;
  real ols_intercept;
  vector[Z] ols_slope;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Parameter Block
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
parameters {
  real <lower = 0> a;
  //real <lower = 0.1386294,upper=693.1472> a; //Lower = hl = 5, upper = hl = 0.001
  real <lower = 0, upper = variance(Y)*4> sigma2_y; //Added to limit the variance based on Kjetil's suggestion
  //real <lower = 0> sigma2_y;
  real alpha; //OU alpha
  vector[Z] beta; //OU beta

  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Transformed Parameter Block
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
transformed parameters {
  matrix[N,Z+1] X_evol;
  vector[Z+1] beta_evol_complete;
  vector[Z_random] beta_evol;
  vector[Z_random] beta_optimal;
  vector[Z_direct] beta_direct;

  matrix[N,N] V_ev;
  matrix[N,N] V_me_ev;
  matrix[N,N] Vt_ev;
  
  beta_evol_complete = rep_vector(0,Z+1);
  beta_evol = rep_vector(0,Z_random);
  beta_optimal = rep_vector(0,Z_random);
  beta_direct = rep_vector(0,Z_direct);

  Vt_ev = varcov_model(N,  tij,  tja,  ta, random_cov, sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term, Z, Z_direct, Z_random);
  
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    V_me_ev = varcov_measurement(N, ta, direct_cov, mv_direct_cov, mv_random_cov, sigma_squared_x, beta, Z, Z_direct, Z_random, a, T_term);
    
    
    if(sum(mv_random_cov)!= 0){
      V_ev = Vt_ev + V_me_ev + diag_matrix(mv_response); //Hansen & Bartoszek (2012) - Eq 10
    }else{
        V_ev = Vt_ev + V_me_ev; //Hansen & Bartoszek (2012) - Eq 10
    }
  }else if(sum(mv_response)!= 0){
      V_me_ev = rep_matrix(0,N,N);
      V_ev = Vt_ev + diag_matrix(mv_response); //Hansen & Bartoszek (2012) - Eq 10
    }else{
      V_me_ev = rep_matrix(0,N,N);
      V_ev = Vt_ev;  
    }
  
    
  //Calculate evolutionary regression slope
  X_evol = design_matrix_evol(N, a,  T_term, direct_cov, random_cov, Z, Z_direct, Z_random);
  beta_evol_complete = inverse(X_evol'*inverse(V_ev)*X_evol)*(X_evol'*inverse(V_ev)*Y); //Hansen et al. 2008
  
    
  if(sum(direct_cov) != 0){
    beta_direct = beta[1:Z_direct];}

  if(sum(random_cov) != 0){
    beta_evol = beta_evol_complete[Z_direct+2:Z_direct+Z_random+1];
    beta_optimal = beta[Z_direct+1:Z_direct+Z_random];}
  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Model Block
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
model {
  //Declare variables
  //real a;
  //real sigma2_y;
  matrix[N,N] V;
  matrix[N,N] Vt;
  matrix[N,N] V_me;

  vector[N] mu;
  matrix[N,Z] X;
  matrix[N, N] L_V;
  //array[1] matrix[N,N];
  //matrix[N,Z] X_evol;
  //vector[Z] beta_evol;

//Priors
  a ~ lognormal(1.0,1.0); //a = log(2)/half-life
  //sigma2_y ~ exponential(0.1); //
  alpha ~ normal(ols_intercept,0.5); //Intercept - orginally 0.5
  beta ~ normal(ols_slope, 0.4); //Slope originally 0.4

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //a = log(2)/hl;
  //sigma2_y = vy*(2*a);

  //////////////////////////////////////////////////////////////////////////////////////////////////////
//Regression - either constraint for direct cov or adaptive for random cov
//
//Set up X matrix
  //print("OK");
  X = design_matrix(N, a,  T_term, direct_cov,random_cov, Z, Z_direct, Z_random);
  
//Set up V matix
  Vt = varcov_model(N,  tij,  tja,  ta, random_cov, sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term, Z, Z_direct, Z_random);
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    V_me = varcov_measurement(N, ta, direct_cov, mv_direct_cov,mv_random_cov, sigma_squared_x, beta, Z, Z_direct, Z_random, a, T_term);
    if(sum(mv_random_cov)!= 0){
      V = Vt + V_me + diag_matrix(mv_response); //Hansen & Bartoszek (2012) - Eq 10
    }else{
        V = Vt + V_me; //Hansen & Bartoszek (2012) - Eq 10
      }
  }else if(sum(mv_response)!= 0){
      V_me = rep_matrix(0,N,N);
      V = Vt + diag_matrix(mv_response); //Hansen & Bartoszek (2012) - Eq 10
    }else{
      V_me = rep_matrix(0,N,N);
      V = Vt;  
    }
  L_V = cholesky_decompose(V);
  
//OU with random covariates
  mu = X*beta+alpha;
  Y ~ multi_normal_cholesky(mu , L_V);

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//Generated Quantities Block
////////////////////////////////////////////////////////////////////////////////////////////////////////////
generated quantities {
  real <lower = 0> vy;
  //real <lower = 0> sigma2_y;
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
  
//////////
  hl = log(2)/a;
  //a = log(2)/hl;
  vy = sigma2_y/(2*a);
  //sigma2_y = vy*(2*a);
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Calculate r2 based on constraint or adaptive regression

  X_opt = design_matrix(N, a,  T_term, direct_cov,random_cov,  Z,  Z_direct,  Z_random);
  Vt_final = varcov_model(N,  tij,  tja,  ta, random_cov, sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term, Z, Z_direct, Z_random);
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    //print("sum!=0");
    V_me_final = varcov_measurement(N, ta, direct_cov, mv_direct_cov, mv_random_cov, sigma_squared_x, beta, Z, Z_direct, Z_random, a, T_term);
    V_final = Vt_final + V_me_final + diag_matrix(mv_response);

    }
  else{
    //print("sum=0");
    V_me_final = rep_matrix(0,N,N);
    V_final = Vt_final;
    }

  
  pred_mean = (X_opt*beta+alpha);
  grand_mean = ((rep_vector(1,N))' * inverse(V_final) * Y) / sum(inverse(V_final));
  sst = ((Y - grand_mean)' * inverse(V_final) * (Y - grand_mean));
  sse = ((Y - pred_mean)' * inverse(V_final) * (Y - pred_mean));
  r_squared = (sst - sse) / sst;

}

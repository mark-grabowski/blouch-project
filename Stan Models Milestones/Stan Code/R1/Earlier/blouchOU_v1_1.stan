functions {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Blouch v1.1 - 07.20.2022
//This code estimates evolutionary parameters without prediction
//Adding in code to do multivariate analyses (multiple xs) for SB R1
//Works with blouchOU.setup.mv.R function
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code - Adaptive
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    matrix design_matrix(int N, int Z, real a, vector T_term, matrix direct_cov, matrix random_cov){
    vector[N] rho;
    matrix[N,Z] X;
    matrix[N,Z] rhos;
    rho = (1 - (1 - exp(-a * T_term))./(a * T_term)); //For OU model
    rhos = rep_matrix(rho,Z);

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
    matrix design_matrix_evol(int N, int Z, real a, vector T_term, matrix direct_cov, matrix random_cov){
    matrix[N,Z+1] X;
   
    
    if(sum(random_cov)==0){
      X = append_col(rep_vector(1,N), direct_cov);
      }
    else if(sum(direct_cov)==0){
      X = append_col(rep_vector(1,N), random_cov);
    }else{
      X = append_col(rep_vector(1,N),direct_cov);
      X = append_col(X,random_cov);
      }
    return(X);
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Vt function - Biological variance based on phylogeny and rate of adaptation - Based on Hansen et al. (2008)
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

    if(sum(random_cov) != 0){
      s1 = sum(sigma_squared_x * beta1sq);
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
matrix varcov_measurement(int N, int Z, matrix ta,  matrix direct_cov, matrix mv_direct_cov, matrix mv_random_cov, matrix sigma_squared_x, vector beta){
  matrix[N,N] Vur[Z];
  matrix[N,N] Vxtr[Z];
  matrix[N,N] Vud[Z];
  matrix[N,N] Vxtd[Z];
  matrix[N,N] Vxt[Z];
  matrix[N,N] Vu[Z];
  matrix[N,N] Vx[Z];
  matrix[N,N] Vu_given_x[Z];
  matrix[N,N] beta2_Vu_given_x[Z];
  matrix[N,N] beta2_Vu_given_x_sum;
  matrix[N,N] dummy;
  
  dummy = rep_matrix(0,N,N);
  Vur = rep_array(dummy,Z);
  Vud = rep_array(dummy,Z);
  Vxtd = rep_array(dummy,Z);
  Vxtr = rep_array(dummy,Z);

  if(sum(mv_random_cov) != 0){
    for (i in 1:Z){
      Vur[i] = diag_matrix(mv_random_cov[,i]); //
      Vxtr[i] = ta * sigma_squared_x[1,i];
      }
    }
  if(sum(mv_direct_cov) != 0){
    for (i in 1:Z){
      Vud[i] = diag_matrix(mv_direct_cov[,i]); //
      Vxtd[i] = diag_matrix(rep_vector(variance(direct_cov[,i]),N))-Vud[i];
      }
  }
  for (i in 1:Z){
    Vxt[i] = Vxtd[i]+Vxtr[i];
    Vu[i] = Vud[i]+Vur[i];
    Vx[i] = Vxt[i] + Vu[i];
    Vu_given_x[i] = Vu[i] - Vu[i] * inverse(Vx[i]) * Vu[i];
    beta2_Vu_given_x[i] = Vu_given_x[i] * square(beta[i]);
    }
  
  for (j in 1:N){
    for (k in 1:N){
        beta2_Vu_given_x_sum[j,k] = sum(beta2_Vu_given_x[,j,k]);
    }
      
    }
  //print(beta2_Vu_given_x[1,1,1]);
  //print(beta2_Vu_given_x[2,1,1]);
  //print(beta2_Vu_given_x_sum[1,1]);
  
  
  return(beta2_Vu_given_x_sum);
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Data Block
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
data {
  //Extant data
  int N; //species number
  int Z; //number of predictors
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
  vector[Z+1] beta_evol;
  matrix[N,N] V_ev;
  matrix[N,N] V_me_ev;
  matrix[N,N] Vt_ev;

  Vt_ev = varcov_model(N,  tij,  tja,  ta, random_cov, Z,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term);
  
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    V_me_ev = varcov_measurement(N, Z, ta, direct_cov, mv_direct_cov, mv_random_cov, sigma_squared_x, beta);
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
  X_evol = design_matrix_evol(N, Z,  a,  T_term, direct_cov,random_cov);
  beta_evol = inverse(X_evol'*inverse(V_ev)*X_evol)*(X_evol'*inverse(V_ev)*Y); //Hansen et al. 2008
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
  alpha ~ normal(ols_intercept,0.5); //Intercept
  beta ~ normal(ols_slope, 0.4); //Slope

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //a = log(2)/hl;
  //sigma2_y = vy*(2*a);

  //////////////////////////////////////////////////////////////////////////////////////////////////////
//Regression - either constraint for direct cov or adaptive for random cov
//
//Set up X matrix
  //print("OK");
  X = design_matrix(N,  Z, a,  T_term, direct_cov,random_cov);
  
//Set up V matix
  Vt = varcov_model(N,  tij,  tja,  ta, random_cov, Z,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term);
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    V_me = varcov_measurement(N, Z, ta, direct_cov, mv_direct_cov,mv_random_cov, sigma_squared_x, beta);
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

  X_opt = design_matrix(N,  Z,  a,  T_term, direct_cov,random_cov);
  Vt_final = varcov_model(N,  tij,  tja,  ta, random_cov, Z,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term);
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    //print("sum!=0");
    V_me_final = varcov_measurement(N, Z, ta, direct_cov, mv_direct_cov, mv_random_cov, sigma_squared_x, beta);
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

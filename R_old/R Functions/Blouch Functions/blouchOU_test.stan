functions {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Blouch v1 - 11.23.2020
//This version only recreates Slouch output, with no fossil prediction or within species data
//To be renamed Blouch_OU for ms
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  matrix design_matrix(int N, int evol, real a, vector T_term, matrix direct_cov, matrix random_cov, int Z){
    matrix[N,Z] rho;
    matrix[N,Z] X;
    if(evol==0) {
      rho = to_matrix(1 - (1 - exp(-a * T_term))./(a * T_term)); //For OU model
    }
    else{
      rho=to_matrix(rep_vector(1.0,N)); //For Evolutionary regression
    }
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
//    if (a < 1e-14){ //Original Bjorn code - is this causing weird Slouch numbers?
//      if(sum(random_cov) == 0){ //Direct covariates and small a
//        Vt = sigma2_y * ta;
//        print("Direct - small a");
//      }else{
//        s1 = sum(sigma_squared_x * beta1sq);
//        print("Random - small a");
//        Vt = sigma2_y * ta + s1 * ta .* ((ta .* ta) ./ 12 + tja .* (tja') ./ 4);
//      }
//    }
//  else{ //Random + Direct covariates, or just random
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
//  }    
  
  
  return(Vt);
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//V_me function
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrix varcov_measurement(int N, int Z, matrix ta,  matrix direct_cov, matrix mv_direct_cov, matrix mv_random_cov, matrix sigma_squared_x, vector beta){
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


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
parameters {
  real <lower = 0> a; //Lower = hl = 4, upper = hl = 0.0010
  //real <lower = 0.1386294,upper=693.1472> a; //Lower = hl = 5, upper = hl = 0.001
  //real <lower = 0, upper = variance(Y)*4> sigma2_y; //Added to limit the variance based on Kjetil's suggestion
  real <lower = 0> sigma2_y;
  //real <lower = 0> sigma2_y;
  real alpha; //OU alpha
  vector[Z] beta; //OU beta

  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
transformed parameters {
  
  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
model {
//Declare variables
  //real a;
  //real sigma2_y;
  vector[N] mu;
  matrix[N,N] V;
  matrix[N,N] Vt;
  matrix[N,N] V_me;
  matrix[N,Z] X;
  matrix[N, N] L_V;
  
//Priors
  //a ~ lognormal(log(1.25),log(1.25)); 
  a ~ lognormal(1.0,1.15); 
  //sigma2_y ~ cauchy(0,1.0);
  sigma2_y ~ exponential(0.1);
  alpha ~ normal(ols_intercept,0.025); //Simulations PREVIOUS CODE 0.1 for both
  beta ~ normal(ols_slope,0.075); //Simulations

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //a = log(2)/hl;
  //sigma2_y = vy*(2*a);

  //////////////////////////////////////////////////////////////////////////////////////////////////////
//Regression - either constraint for direct cov or adaptive for random cov
//
//Set up X matrix
  X = design_matrix( N,  0,  a,  T_term, direct_cov,random_cov, Z);
  
//Set up V matix
  Vt = varcov_model(N,  tij,  tja,  ta, random_cov, Z,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term);
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
  //mu = X*beta;
  Y ~ multi_normal_cholesky(mu , L_V);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

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
  matrix[N,Z] X_evol;
  vector[Z] beta_evol;
  
//////////
  hl = log(2)/a;
  //a = log(2)/hl;
  vy = sigma2_y/(2*a);
  //sigma2_y = vy*(2*a);
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Calculate r2 based on constraint or adaptive regression
  X_opt = design_matrix(N,  0,  a,  T_term, direct_cov,random_cov, Z);
  Vt_final = varcov_model(N,  tij,  tja,  ta, random_cov, Z,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term);
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
  
  
//#############################################################################################################################
//Calculate evolutionary regression slope
  X_evol = design_matrix(N,  1,  a,  T_term, direct_cov,random_cov, Z);
  beta_evol = inverse(X_evol'*inverse(V_final)*X_evol)*X_evol'*inverse(V_final)*Y; //Hansen et al. 2008
}

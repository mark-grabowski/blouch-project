functions {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Blouch BM model
//12.07.2020
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  matrix design_matrix(int N, int evol, real a, vector T_term, matrix direct_cov, matrix random_cov, int Z){
    matrix[N,Z] rho;
    matrix[N,Z] X;

    rho = to_matrix(T_term ./ 2); //BM Model

    if(sum(random_cov)==0){
      X = direct_cov;
      }
    else if(sum(direct_cov)==0){
      X = random_cov .* rho;
    }else{
      X = append_col(direct_cov, random_cov .* rho);
      }
    return(X);
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

  if(sum(random_cov) == 0){ //Direct covariates only and a = 0
        Vt = sigma2_y * ta;
    }else{
      s1 = sum(sigma_squared_x * beta1sq); //Random covariates and a = 0
      Vt = sigma2_y * ta + s1 * ta .* ((ta .* ta) ./ 12 + tja .* (tja') ./ 4);
    }
  return(Vt);
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
}


parameters {
  real <lower = 0> sigma2_y;
  real alpha;
  vector[Z] beta;
}

transformed parameters {
  }

model {
//Declare variables
  vector[N] mu;
  matrix[N,N] V;
  matrix[N,N] Vt;
  matrix[N,N] V_me;
  matrix[N, N] L_V;
  real a = 0;
  matrix[N,Z] X;


//Priors
  //sigma2_y ~ exponential(0,1.0);
  alpha ~ normal(ols_intercept,0.5); //Intercept
  beta ~ normal(ols_slope,0.4); //Slope


  //Set up X matrix
  X = design_matrix( N,  0,  a,  T_term, direct_cov,random_cov, Z);

//Set up V matix
  Vt = varcov_model(N,  tij,  tja,  ta,  random_cov, Z,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term);
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    V_me = varcov_measurement(N, Z, ta, direct_cov, mv_direct_cov,mv_random_cov, sigma_squared_x, beta);
  }
  else{
    V_me = rep_matrix(0,N,N);
  }  //print(V_me[1,1],V_me[1,2]);
  V = Vt + V_me + diag_matrix(mv_response);
  L_V = cholesky_decompose(V);

//OU with random covariates
  mu = X*beta+alpha;
  Y ~ multi_normal_cholesky(mu , L_V);
}

generated quantities {
  vector[N] pred_mean;
  real grand_mean;
  real sst;
  real sse;
  real r_squared;
  real a = 0;


  matrix[N,Z] X_opt;
  matrix[N,N] V_final;
  matrix[N,N] V_me_final;
  matrix[N,N] Vt_final;
  matrix[N, N] L_V_final;

  //Calculate r2 based on regression
  X_opt = design_matrix(N,  0,  a,  T_term, direct_cov,random_cov, Z);
  Vt_final = varcov_model(N,  tij,  tja,  ta,  random_cov, Z,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term);
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    V_me_final = varcov_measurement(N, Z, ta, direct_cov, mv_direct_cov, mv_random_cov, sigma_squared_x, beta);
    }
  else{
    V_me_final = rep_matrix(0,N,N);
    }
  V_final = Vt_final + V_me_final + diag_matrix(mv_response);

  L_V_final = cholesky_decompose(V_final);
  pred_mean = (X_opt*beta+alpha);
  grand_mean = ((rep_vector(1,N))' * (V_final') * Y) / sum(inverse(V_final));
  sst = ((Y - grand_mean)' * inverse(V_final) * (Y - grand_mean));
  sse = ((Y - pred_mean)' * inverse(V_final) * (Y - pred_mean));
  r_squared = (sst - sse) / sst;
}

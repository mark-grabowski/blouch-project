//Blouch OU model reprogrammed - 2023
//Adaptive model accounts for measurement error on X and/or Y variables
//See below for which lines to comment/uncomment based on whether measurement error is present
//Prior Predictive Check Code

functions {
  matrix calc_dmX(real a, vector T_term, matrix X){
    int N = dims(X)[1];
    int Z_adapt = dims(X)[2];
    vector[N] rho = (1 - (1 - exp(-a * T_term))./(a * T_term));
    matrix[N,Z_adapt] rhos = rep_matrix(rho,Z_adapt);
    matrix[N,Z_adapt] dmX = X .* rhos;
    return(dmX);
  }
  matrix calc_V(real a,real sigma2_y,matrix ta, matrix tij, matrix tja, vector T_term, vector beta, matrix sigma2_x) {
    int N = dims(ta)[1];
    int Z_adapt = dims(beta)[1];
    vector[Z_adapt] ones = rep_vector(1,Z_adapt);
    matrix[N,N] ti = rep_matrix(T_term,N);
    matrix[N,N] term0;
    matrix[N,N] term1;
    matrix[N,N] term2;
    matrix[N,N] Vt;
    real var_opt;
    if(Z_adapt==1){var_opt = beta[1] * beta[1] * sigma2_x[1,1];
//    }else{var_opt = (beta[1:Z_adapt] .* beta[1:Z_adapt])' * sigma2_x * ones;}
    }else{var_opt = (beta[1:Z_adapt])' * sigma2_x * ones;}
    term0 = ((var_opt + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) .* exp(-a * tij);
    term1 = (1 - exp(-a * ti)) ./ (a * ti);
    term2 = exp(-a * tja) .* (1 - exp(-a * ti)) ./ (a * ti);
    Vt = term0 + var_opt * (ta .* term1 .* (term1') - ((1 - exp(-a * ta)) ./ a) .* (term2 + (term2'))); //From Hansen et al. (2008)
    return Vt;
  }
}
data {
  int N;
  int Z_adapt;
  vector[N] Y_obs;
  matrix[N,Z_adapt] X_obs;
  vector[N] Y_error;
  matrix[N,Z_adapt] X_error;
  matrix[N,N] ta;
  matrix[N,N] tij;
  matrix[N,N] tja;
  vector[N] T_term;
  matrix[Z_adapt,Z_adapt] sigma2_x;
  vector[2] hl_prior;
  real vy_prior;
  vector[2] optima_prior;
  vector[2] beta_prior;

}
parameters {

}
transformed parameters{
}
model {

}
generated quantities {
  matrix[N,N] V;
  matrix[N,N] L_v;
  matrix[N,Z_adapt] dmX;
  vector[N] mu;
  real<lower=0> hl = lognormal_rng(hl_prior[1],hl_prior[2]);
  real<lower=0> vy = exponential_rng(vy_prior);
  real sigma2_y = vy*(2*(log(2)/hl));
  real a = log(2)/hl;
  vector[N] Y_sim;
  vector[N] Y_sim_obs;
  matrix[N,Z_adapt] X_sim;
  vector[Z_adapt] beta_sim;
  real optima_sim = normal_rng(optima_prior[1],optima_prior[2]); //intercept from OLS
  for(i in 1:Z_adapt){
     beta_sim[i]= normal_rng(beta_prior[1],beta_prior[2]);
  }
  for(i in 1:(Z_adapt)){//Given measurement error in X variable, uncomment this nested statement
    for(j in 1:N){
      X_sim[j,i] = normal_rng(X_obs[j,i], X_error[j,i]);
    }
  }
  V = calc_V(a,sigma2_y,ta,tij,tja,T_term,beta_sim,sigma2_x);
  L_v = cholesky_decompose(V);
  dmX = calc_dmX(a,T_term,X_sim);//Given measurement error in X variable, uncomment this statement
  mu = optima_sim+dmX*beta_sim;
  Y_sim = multi_normal_cholesky_rng(mu , L_v);//Given measurement error in Y variable, uncomment this statement
  for(i in 1:N){
    Y_sim_obs[i] = normal_rng(Y_sim[i],Y_error[i]); //Given measurement error in Y variable, uncomment this statement
  }
}


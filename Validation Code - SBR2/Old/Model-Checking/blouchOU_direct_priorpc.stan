//Blouch OU model reprogrammed
//Direct effect model
//Using Hansen (1997)
//With Statistical Rethinking ME

functions {
  matrix calc_V( real a,real sigma2_y,matrix ta, matrix tij) {
    int N = dims(ta)[1];
    matrix[N, N] Vt;
    Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)
    return Vt;
  }
}
data {
  int N;
  int Z_direct;
  vector[N] Y_obs;
  matrix[N,Z_direct] X_obs;
  vector[N] Y_error;
  matrix[N,Z_direct] X_error;
  matrix[N,N] ta;
  matrix[N,N] tij;
  vector[2] hl_prior;
  real vy_prior;
  vector[2] optima_prior;
  vector[2] beta_prior;
}
parameters {

}
model {

}
generated quantities {
  vector[N] mu;
  matrix[N,N] V;
  matrix[N,N] L_v;
  matrix[N,Z_direct] X_sim;
  vector[N] Y_sim;
  vector[N] Y_sim_obs;
  real<lower=0> hl = lognormal_rng(hl_prior[1],hl_prior[2]);
  real<lower=0> vy = exponential_rng(vy_prior);
  real sigma2_y = vy*(2*(log(2)/hl));
  real a = log(2)/hl;
  real optima = normal_rng(optima_prior[1],optima_prior[2]);
  vector[Z_direct] beta;
  for (i in 1:Z_direct){
    beta[i]= normal_rng(beta_prior[1],beta_prior[2]);
  }
  for(i in 1:(Z_direct)){
    for(j in 1:N){
      X_sim[j,i] = normal_rng(X_obs[j,i], X_error[j,i]);
    }
  }
  mu = optima+X_sim*beta;//Given measurement error in X variable, uncomment this statement
  V = calc_V(a, sigma2_y,ta, tij);
  L_v = cholesky_decompose(V);
  Y_sim = multi_normal_cholesky_rng(mu , L_v);//Given measurement error in Y variable, uncomment this statement
  for(i in 1:N){
    Y_sim_obs[i] = normal_rng(Y_sim[i],Y_error[i]); //Given measurement error in Y variable, uncomment this statement
  }
}

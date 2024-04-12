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
  real <lower = 0> hl;
  vector[Z_direct] beta;
  real optima;
  real <lower=0> vy;
  vector[N] Y;
  matrix[N,Z_direct] X;
}
model {
  vector[N] mu;
  real a;
  matrix[N,N] V;
  matrix[N,N] L_v;
  real sigma2_y = vy*(2*(log(2)/hl));
  target += lognormal_lpdf(hl|hl_prior[1],hl_prior[2]);
  target += exponential_lpdf(vy|vy_prior);
  target += normal_lpdf(optima|optima_prior[1],optima_prior[2]);
  target += normal_lpdf(beta|beta_prior[1],beta_prior[2]);
  a = log(2)/hl;
  V = calc_V(a, sigma2_y,ta, tij);
  L_v = cholesky_decompose(V);
  for(i in 1:Z_direct){//Given measurement error in X variable, uncomment this nested statement
    target += normal_lpdf(X[,i]|0,1);
    target += normal_lpdf(X_obs[,i]|X[,i],X_error[,i]);
  }
  mu = optima+X*beta;//Given measurement error in X variable, uncomment this statement
  target += multi_normal_cholesky_lpdf(Y | mu , L_v);
  target += normal_lpdf(Y_obs | Y, Y_error);

}
generated quantities {
  matrix[N,N] V;
  matrix[N,N] L_v;
  vector[N] mu;
  vector[N] Y_sim;
  vector[N] Y_sim_obs;
  real sigma2_y = vy*(2*(log(2)/hl));
  real a = log(2)/hl;
  V = calc_V(a, sigma2_y,ta, tij);
  L_v = cholesky_decompose(V);
  mu = optima+X*beta;//Given measurement error in X variable, uncomment this statement
  Y_sim = multi_normal_cholesky_rng(mu , L_v);//Given measurement error in Y variable, uncomment this statement
  for(i in 1:N){
    Y_sim_obs[i] = normal_rng(Y_sim[i],Y_error[i]); //Given measurement error in Y variable, uncomment this statement
  }
}

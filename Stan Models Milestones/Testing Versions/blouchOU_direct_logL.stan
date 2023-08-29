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
  int Z; 
  vector[N] Y_obs;
  matrix[N,Z] X_obs;
  vector[N] Y_error;
  matrix[N,Z] X_error;
  matrix[N,N] ta;
  matrix[N,N] tij;
}
parameters {
  real <lower = 0> hl;
  vector[Z] beta; 
  real alpha;
  //real <lower=0> sigma2_y;
  real <lower=0> vy;
  vector[N] Y;
  matrix[N,Z] X;
}
model {
  vector[N] mu;
  real a;
  matrix[N,N] V;
  matrix[N,N] L_v;
  //sigma2_y ~ exponential(1);
  real sigma2_y = vy*(2*(log(2)/hl));
  hl ~ lognormal(log(0.25),0.75);
  vy ~ exponential(1);
  alpha ~ normal(-1.179507,0.75);
  beta ~ normal(6.304451,1.5);
  a = log(2)/hl;
  V = calc_V(a, sigma2_y,ta, tij);
  L_v = cholesky_decompose(V);
  for(i in 1:Z){//Given measurement error in X variable, uncomment this nested statement
    X[,i] ~ normal(0,1);  
    X_obs[,i] ~ normal(X[,i], X_error[,i]);
  }
  mu = alpha+X*beta;//Given measurement error in X variable, uncomment this statement
  //mu = alpha+X_obs*beta;//Given nomeasurement error in X variable, uncomment this statement
  Y ~ multi_normal_cholesky(mu , L_v);//Given measurement error in Y variable, uncomment this statement
  Y_obs ~ normal(Y,Y_error); //Given measurement error in Y variable, uncomment this statement
  //Y_obs ~ multi_normal_cholesky(mu , L_v); //Given no measurement error in Y variable, uncomment this statement
}
generated quantities {
  //real vy = sigma2_y/(2*(log(2)/hl));
  matrix[N,N] V;
  matrix[N,N] inv_V;
  vector[N] mu;
  real g_i;
  real sigma_ii;
  real sigma_i;
  real u_i;
  vector[N] log_lik;
  real sigma2_y = vy*(2*(log(2)/hl));
  real a = log(2)/hl;
  //Based on https://cran.r-project.org/web/packages/loo/vignettes/loo2-non-factorized.html#loo-cv-for-multivariate-normal-models
  //LOO-CV for multivariate normal models
  V = calc_V(a, sigma2_y,ta, tij);
  inv_V = inverse(V);
  mu = alpha+X*beta;//Given measurement error in X variable, uncomment this statement

  for(i in 1:N){
      g_i = (inv_V*(Y_obs-mu))[i];
      sigma_ii = inv_V[i,i];
      u_i = Y_obs[i]-g_i/sigma_ii;
      sigma_i = 1/sigma_ii;
      
      log_lik[i] = -0.5*log(2*pi()*sigma_i)-0.5*(square(Y_obs[i]-u_i)/sigma_i);
      }
}

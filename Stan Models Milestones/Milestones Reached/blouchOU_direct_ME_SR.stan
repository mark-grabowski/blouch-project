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
  vector<lower=0>[Z] beta; 
  real alpha;
  real <lower=0> sigma2_y;
  vector[N] Y;
  matrix[N,Z] X;
}
model {
  vector[N] mu;
  real a;
  matrix[N,N] V;
  matrix[N,N] L_v;
  hl ~ lognormal(log(0.4),1); 
  sigma2_y ~ exponential(1);
  alpha ~ normal(4,0.2);
  beta ~ lognormal(log(0.4),1); 
  a = log(2)/hl;
  V = calc_V(a, sigma2_y,ta, tij);
  L_v = cholesky_decompose(V);
  for(i in 1:Z){
    X[,i] ~ normal(0,1);  
    X_obs[,i] ~ normal(X[,i], X_error[,i]);
  }
  mu = alpha+X*beta;  
  Y ~ multi_normal_cholesky(mu , L_v);
  Y_obs ~ normal(Y,Y_error);
}
generated quantities {
  real vy = sigma2_y/(2*(log(2)/hl));
}
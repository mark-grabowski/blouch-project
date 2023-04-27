//Blouch OU model reprogrammed
//Direct effect model
//Using Hansen (1997) 

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
  vector[N] Y;
  vector[N] X;
  matrix[N,N] ta;
  matrix[N,N] tij;
  }


parameters {
  real <lower = 0> hl;
  real beta; 
  real alpha;
  real <lower=0> sigma2_y;
}


model {
  matrix[N,N] V;
  vector[N] mu;
  real a;
  matrix[N,N] K_v;
  matrix[N,N] L_v;
  
  hl ~ lognormal(log(0.4),1); 
  sigma2_y ~ exponential(1);
  alpha ~ normal(4,0.1);
  beta ~ normal(0,0.2);
  
  a = log(2)/hl;
  K_v = calc_V(a, sigma2_y,ta, tij);
  
  L_v = cholesky_decompose(K_v);
  mu = alpha+beta*X;  
  Y ~ multi_normal_cholesky(mu , L_v);
}
generated quantities {
  real vy = sigma2_y/(2*(log(2)/hl));
}

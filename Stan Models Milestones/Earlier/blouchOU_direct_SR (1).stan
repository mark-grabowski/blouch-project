functions {
      matrix calc_V(matrix Dmat, real var_anc, real alpha) {
        int N = dims(Dmat)[1];
        matrix[N, N] K;
        for (i in 1:N) {
          for (j in 1:N) {
            K[i, j] = var_anc * exp(-alpha * Dmat[i,j] );
            K[j, i] = K[i, j];
          }
        }
        return K;
    }
}
  
data {
  int N; 
  int Z; 
  vector[N] Y;
  matrix[N,Z] X;
  matrix[N,N] Dmat;
  }


parameters {
  real <lower = 0> hl;
  vector[Z] beta; 
  real alpha;
  real <lower=0> var_anc;
}


model {
  matrix[N,N] V;
  vector[N] mu;
  real a;
  matrix[N,N] K_v;
  matrix[N,N] L_v;
  
  hl ~ lognormal(log(0.25),0.75); 
  //vy ~ exponential(1);
  var_anc ~ exponential(0.1);
  alpha ~ normal(4,0.2);
  beta ~ normal(0,0.1);
  
  a = log(2)/hl;
  K_v = calc_V(Dmat, var_anc, a);
  L_v = cholesky_decompose(K_v);
  mu = alpha+X*beta;  
  Y ~ multi_normal_cholesky(mu , L_v);
}
generated quantities {

}

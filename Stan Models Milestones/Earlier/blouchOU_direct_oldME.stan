//Blouch OU model reprogrammed
//Direct effect model with Measurement Error
//Using Hansen (1997) 

functions {
  matrix calc_Vt(int N, real a,real sigma2_y,matrix ta, matrix tij) {
    //int N = dims(ta)[1];
    matrix[N, N] Vt;
    Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)
    return Vt;
  }
  matrix calc_Ve(int N, vector Y_error){
    matrix[N,N] V_e= diag_matrix(Y_error);
    return(V_e);
  }
  matrix calc_beta2V_ux(int N, vector beta, matrix X_error){
      int Z = dims(X_error)[2];
      matrix[N,N] Vud[Z] = rep_array(rep_matrix(0,N,N),Z);
      matrix[N,N] Vxtd[Z] = rep_array(rep_matrix(0,N,N),Z);
      matrix[N,N] Vxd[Z] = rep_array(rep_matrix(0,N,N),Z);
      vector[N] Vu_given_x[Z];
      vector[N] beta2_Vu_given_x[Z];
      vector[N] beta2_Vu_given_x_sum;
      for (i in 1:Z){
        Vud[i] = diag_matrix(X_error[,i]);
        Vxtd[i] = diag_matrix(rep_vector(variance(X_error[,i]),N))-Vud[i];
        Vxd[i] = Vxtd[i] + Vud[i];
        Vu_given_x[i] = diagonal(Vud[i] - Vud[i] * inverse(Vxd[i]) * Vud[i]);
        beta2_Vu_given_x[i] = to_vector(Vu_given_x[i] * square(beta[i]));
        }
      for (i in 1:N){
        beta2_Vu_given_x_sum[i] = sum(beta2_Vu_given_x[,i]);
        }
      return(diag_matrix(beta2_Vu_given_x_sum));
    }
}

data {
  int N; 
  int Z; 
  vector[N] Y;
  matrix[N,Z] X;
  vector[N] Y_error;
  matrix[N,Z] X_error;
  matrix[N,N] ta;
  matrix[N,N] tij;
  }

parameters {
  real <lower = 0> hl;
  vector[Z] beta; 
  real alpha;
  real <lower=0> sigma2_y;
}

model {
  matrix[N,N] V;
  vector[N] mu;
  real a;
  matrix[N,N] V_t;
  matrix[N,N] V_e = rep_matrix(0,N,N);
  matrix[N,N] beta2V_ux = rep_matrix(0,N,N);
  matrix[N,N] L_v;
  hl ~ lognormal(log(0.4),1); 
  sigma2_y ~ exponential(1);
  alpha ~ normal(4,0.2);
  beta ~ normal(0,0.2);
  a = log(2)/hl;
  V_t = calc_Vt(N, a, sigma2_y,ta, tij);
  if(num_elements(Y_error)!=0){
    V_e = calc_Ve(N, Y_error);  
  }
  if(num_elements(X_error)!=0){
    beta2V_ux = calc_beta2V_ux(N, beta, X_error);
  }
  V=V_t+V_e+beta2V_ux;
  L_v = cholesky_decompose(V);
  mu = alpha+X*beta;  
  Y ~ multi_normal_cholesky(mu , L_v);
}
generated quantities {
  real vy = sigma2_y/(2*(log(2)/hl));
}

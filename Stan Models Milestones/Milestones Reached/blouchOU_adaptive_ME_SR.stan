//Blouch OU model reprogrammed
//Adaptive model
//Following Hansen et al. (2008)
//With Statistical Rethinking ME

functions {
  matrix calc_dmX(real a, vector T_term, matrix X){
    int N = dims(X)[1];
    int Z = dims(X)[2];
    vector[N] rho = (1 - (1 - exp(-a * T_term))./(a * T_term));
    matrix[N,Z] rhos = rep_matrix(rho,Z);
    matrix[N,Z] dmX = X .* rhos;
    return(dmX);
  }
  matrix calc_V(real a,real sigma2_y,matrix ta, matrix tij, matrix tja, vector T_term, vector beta, matrix sigma2_x) {
    int N = dims(ta)[1];
    int Z = dims(beta)[1];
    vector[Z] ones = rep_vector(1,Z);
    matrix[N,N] ti = rep_matrix(T_term,N);
    matrix[N,N] term0;
    matrix[N,N] term1;
    matrix[N,N] term2;
    matrix[N,N] Vt;
    real var_opt;
    if(Z==1){var_opt = beta[1] * beta[1] * sigma2_x[1,1];
//    }else{var_opt = (beta[1:Z] .* beta[1:Z])' * sigma2_x * ones;}
    }else{var_opt = (beta[1:Z])' * sigma2_x * ones;}
    term0 = ((var_opt + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) .* exp(-a * tij);
    term1 = (1 - exp(-a * ti)) ./ (a * ti); 
    term2 = exp(-a * tja) .* (1 - exp(-a * ti)) ./ (a * ti);
    Vt = term0 + var_opt * (ta .* term1 .* (term1') - ((1 - exp(-a * ta)) ./ a) .* (term2 + (term2'))); //From Hansen et al. (2008)
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
  matrix[N,N] tja;
  vector[N] T_term;
  matrix[Z,Z] sigma2_x;
  
}
parameters {
  real <lower = 0> hl;
  vector<lower=0>[Z] beta; 
  real alpha;
  real <lower=0> sigma2_y;
  vector[N] Y;
  matrix[N,Z] X;
}
transformed parameters{
}
model {
  vector[N] mu;
  real a;
  matrix[N,N] V;
  matrix[N,N] L_v;
  matrix[N,Z] dmX;
  hl ~ lognormal(log(0.4),1); 
  sigma2_y ~ exponential(1);
  alpha ~ normal(4,0.2);
  beta ~ lognormal(log(0.4),1);
  a = log(2)/hl;
  V = calc_V(a,sigma2_y,ta,tij,tja,T_term,beta,sigma2_x);
  L_v = cholesky_decompose(V);
  for(i in 1:Z){
    X[,i] ~ normal(0,1);  
    X_obs[,i] ~ normal(X[,i], X_error[,i]);
  }
  dmX = calc_dmX(a,T_term,X);
  mu = alpha+dmX*beta;  
  Y ~ multi_normal_cholesky(mu , L_v);
  Y_obs ~ normal(Y,Y_error);
}
generated quantities {
  real a = log(2)/hl;
  real vy = sigma2_y/(2*(log(2)/hl));
  vector[N] rho = (1 - (1 - exp(-a * T_term))./(a * T_term)); 
  vector[Z] beta_e;
  for(i in 1:Z){
    beta_e[i] = beta[i]* rho[i]; 
    }  
}

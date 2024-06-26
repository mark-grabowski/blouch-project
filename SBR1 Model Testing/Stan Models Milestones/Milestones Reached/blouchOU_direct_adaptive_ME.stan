//Blouch OU model reprogrammed
//Combination of Direct Adaptive Predictors
//Following Hansen et al. (2008)
//With Statistical Rethinking ME

functions {
  matrix calc_mixed_dmX(real a, vector T_term, matrix X, int Z_direct, int Z_adaptive){
    int N = dims(X)[1];
    int Z = dims(X)[2];
    vector[N] rho = (1 - (1 - exp(-a * T_term))./(a * T_term));
    matrix[N,Z_adaptive] rhos = rep_matrix(rho,Z_adaptive);
    matrix[N,Z] dmX = append_col(X[,1:Z_direct],X[,Z_direct+1:Z_adaptive+Z_direct] .* rhos);
    return(dmX);
  }
  matrix calc_direct_V( real a,real sigma2_y,matrix ta, matrix tij) {
        int N = dims(ta)[1];
        matrix[N, N] Vt;
        Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)
        return Vt;
  }
  matrix calc_adaptive_V(real a,real sigma2_y,matrix ta, matrix tij, matrix tja, vector T_term, vector beta, matrix sigma2_x) {
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
    }else{var_opt = beta[1:Z]' * sigma2_x * ones;}
    term0 = ((var_opt + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) .* exp(-a * tij);
    term1 = (1 - exp(-a * ti)) ./ (a * ti); 
    term2 = exp(-a * tja) .* (1 - exp(-a * ti)) ./ (a * ti);
    Vt = term0 + var_opt * (ta .* term1 .* (term1') - ((1 - exp(-a * ta)) ./ a) .* (term2 + (term2'))); //From Hansen et al. (2008)
    return Vt;
  }
}
data {
  int N; 
  int Z_direct;
  int Z_adaptive;
  int Z_X_error;
  vector[N] Y_obs;
  matrix[N,Z_direct+Z_adaptive] X_obs;
  vector[N] Y_error;
  matrix[N,Z_X_error] X_error;
  matrix[N,N] ta;
  matrix[N,N] tij;
  matrix[N,N] tja;
  vector[N] T_term;
  matrix[Z_adaptive,Z_adaptive] sigma2_x;
}
parameters {
  real <lower = 0> hl;
  vector<lower=0>[Z_direct+Z_adaptive] beta; //Assuming a positive relationship among traits
  real alpha;
  //real <lower=0> sigma2_y;
  real <lower=0> vy;
  vector[N] Y;
  matrix[N,Z_direct+Z_adaptive] X;
}
transformed parameters{
}
model {
  vector[N] mu;
  real a;
  matrix[N,N] V;
  matrix[N,N] L_v;
  matrix[N,Z_direct+Z_adaptive] dmX;
  real sigma2_y = vy*(2*(log(2)/hl));
  hl ~ lognormal(log(0.25),0.75);
  vy ~ exponential(20);
  alpha ~ normal(2,0.2); //intercept from OLS
  beta ~ normal(0,0.25); 
  a = log(2)/hl;
  for(i in 1:(Z_direct+Z_adaptive)){
    X[,i] ~ normal(0,1);  
    X_obs[,i] ~ normal(X[,i], X_error[,i]);
  }
  dmX = calc_mixed_dmX(a,T_term,X,Z_direct,Z_adaptive);
  V = calc_adaptive_V(a,sigma2_y,ta,tij,tja,T_term,beta[(Z_direct+1):(Z_adaptive+Z_direct)],sigma2_x);
  L_v = cholesky_decompose(V);
  mu = alpha+dmX*beta;  
  Y ~ multi_normal_cholesky(mu , L_v);
  Y_obs ~ normal(Y,Y_error);
}
generated quantities {
  real a = log(2)/hl;
  real sigma2_y = vy*(2*(log(2)/hl));
  vector[N] rho = (1 - (1 - exp(-a * T_term))./(a * T_term)); 
  vector[Z_adaptive] beta_e;
  for(i in 1:Z_adaptive){
    beta_e[i] = beta[Z_direct+i]* rho[i]; 
    }  
}

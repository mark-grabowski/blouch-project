//Blouch OU model reprogrammed - 2023
//Adaptive model accounts for measurement error on X and/or Y variables
//See below for which lines to comment/uncomment based on whether measurement error is present

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
  real <lower = 0> hl;
  vector<lower=0>[Z_adapt] beta; 
  real optima;
  //real <lower=0> sigma2_y;
  real <lower=0> vy;
  vector[N] Y;
  matrix[N,Z_adapt] X;
}
transformed parameters{
}
model {
  vector[N] mu;
  real a;
  matrix[N,N] V;
  matrix[N,N] L_v;
  matrix[N,Z_adapt] dmX;
  real sigma2_y = vy*(2*(log(2)/hl));
  //hl ~ lognormal(log(0.25),0.25);
  target += lognormal_lpdf(hl|hl_prior[1],hl_prior[2]);
  //vy ~ exponential(5);
  target += exponential_lpdf(vy|vy_prior);
  //optima ~ normal(2,0.2); //intercept from OLS
  //beta ~ normal(0,0.25); 
  target += normal_lpdf(optima|optima_prior[1],optima_prior[2]);
  target += normal_lpdf(beta|beta_prior[1],beta_prior[2]);

  a = log(2)/hl;
  V = calc_V(a,sigma2_y,ta,tij,tja,T_term,beta,sigma2_x);
  L_v = cholesky_decompose(V);
  for(i in 1:(Z_adapt)){//Given measurement error in X variable, uncomment this nested statement
    //X[,i] ~ normal(0,1);
    target += normal_lpdf(X[,i]|0,1);
    //X_obs[,i] ~ normal(X[,i], X_error[,i]);
    target += normal_lpdf(X_obs[,i]|X[,i],X_error[,i]);
  }
  dmX = calc_dmX(a,T_term,X);//Given measurement error in X variable, uncomment this statement  
  //dmX = calc_dmX(a,T_term,X_obs);//Given no measurement error in X variable, uncomment this statement
  mu = optima+dmX*beta;  
  //Y ~ multi_normal_cholesky(mu , L_v);//Given measurement error in Y variable, uncomment this statement
  //Y_obs ~ normal(Y,Y_error); //Given measurement error in Y variable, uncomment this statement
  //Y_obs ~ multi_normal_cholesky(mu , L_v); //Given no measurement error in Y variable, uncomment this statement
  target += multi_normal_cholesky_lpdf(Y | mu , L_v);
  target += normal_lpdf(Y_obs | Y, Y_error);

}
generated quantities {
  matrix[N,N] V;
  matrix[N,N] inv_V;
  matrix[N,Z_adapt] dmX;
  vector[N] mu;
  real g_i;
  real sigma_ii;
  real sigma_i;
  real u_i;
  vector[N] log_lik;

  real a = log(2)/hl;
  real sigma2_y = vy*(2*(log(2)/hl));
  vector[N] rho = (1 - (1 - exp(-a * T_term))./(a * T_term)); 
  vector[Z_adapt] beta_e;
  
  for(i in 1:Z_adapt){
    beta_e[i] = beta[i]* rho[i]; 
    }  
  //LOO-CV for multivariate normal models
  dmX = calc_dmX(a,T_term,X);//Given measurement error in X variable, uncomment this statement  
  V = calc_V(a,sigma2_y,ta,tij,tja,T_term,beta,sigma2_x);
  inv_V = inverse(V);
  mu = optima+dmX*beta;  
  for(i in 1:N){
      g_i = (inv_V*(Y_obs-mu))[i];
      sigma_ii = inv_V[i,i];
      u_i = Y_obs[i]-g_i/sigma_ii;
      sigma_i = 1/sigma_ii;
      
      log_lik[i] = -0.5*log(2*pi()*sigma_i)-0.5*(square(Y_obs[i]-u_i)/sigma_i);
      }
}


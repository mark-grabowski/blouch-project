//Blouch OU model reprogrammed - 2023
//Combination of direct effect and adaptive predictors
//Model accounts for measurement error on X and/or Y variables
//See below for which lines to comment/uncomment based on whether measurement error is present

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

}
transformed parameters{
}
model {
 
}
generated quantities {
  matrix[N,N] V;
  matrix[N,N] L_v;
  matrix[N,Z_direct+Z_adaptive] dmX;
  matrix[N,Z_direct+Z_adaptive] X_sim;
  vector[N] Y_sim;
  vector[N] Y_sim_obs;
  vector[N] mu;
  real hl = lognormal_rng(log(0.25),0.75);
  real vy = exponential_rng(1);
  real sigma2_y = vy*(2*(log(2)/hl));
  real a = log(2)/hl;
  vector[Z_direct+Z_adaptive] beta; //Assuming a positive relationship among traits
  real alpha = normal_rng(2,0.2); //intercept from OLS
  for(i in 1:(Z_direct+Z_adaptive)){
     beta[i]= normal_rng(0,0.25);
  }
  for(i in 1:(Z_direct+Z_adaptive)){//Given measurement error in X variable, uncomment this nested statement
    for(j in 1:N){
      X_sim[j,i] = normal_rng(X_obs[j,i], X_error[j,i]);  
    }
  }
  dmX = calc_mixed_dmX(a,T_term,X_sim,Z_direct,Z_adaptive); //Given measurement error in X variable, uncomment this statement
  V = calc_adaptive_V(a,sigma2_y,ta,tij,tja,T_term,beta[(Z_direct+1):(Z_adaptive+Z_direct)],sigma2_x);
  L_v = cholesky_decompose(V);
  mu = alpha+dmX*beta;  
  Y_sim = multi_normal_cholesky_rng(mu , L_v);//Given measurement error in Y variable, uncomment this statement
  for(i in 1:N){
    Y_sim_obs[i] = normal_rng(Y_sim[i],Y_error[i]); //Given measurement error in Y variable, uncomment this statement
  }

}

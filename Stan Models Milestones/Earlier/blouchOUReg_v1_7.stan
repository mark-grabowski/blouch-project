functions {
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Blouch Fixed Niches Code
//08.05.2022 - Revised for SBR1 to be use multivariate predictors, no need to use different random and direct datasets, and can use combo direct and response traits
//09.04.2022 - Added code to allow for correlated predictors - Hansen et al. (2008)
//04.12.2023 - Include new LL code
//04.12.2023 - Trying to get ME to work - may not
//04.12.2023 - Shifted ME to following SR
//////////////////////////////////////////////////////////////////////////////////////////////////////
matrix calc_optima(real a, int n_regimes, int n_lineages, int max_node_length, matrix nodes, matrix nodes_time, 
                     matrix t_end, matrix t_beginning, matrix regime_time, int[,] regimes_matrix){
    matrix[n_lineages,n_regimes] optima_matrix =rep_matrix(0,n_lineages,n_regimes);
      
      for(i in 1:n_lineages){
        int max_j = max_node_length;
        while (max_j > 0 && t_end[i,max_j] == 0) max_j -= 1;
        for(j in 1:max_j){
          optima_matrix[i, regimes_matrix[i,j]] += exp(-a * t_beginning[i,j]) - exp(-a * t_end[i,j]);
          optima_matrix[i, regimes_matrix[i,j]] += exp(-a * nodes_time[i,1]);
        }
      }
      return(optima_matrix);
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code - Adaptive
//////////////////////////////////////////////////////////////////////////////////////////////////////
matrix design_matrix(int N, real a, vector T_term, matrix direct_cov, matrix random_cov,int n_regimes, int n_lineages, 
int max_node_length, matrix nodes, matrix nodes_time, matrix t_end, matrix t_beginning, matrix regime_time, 
int[,] regimes_matrix, int Z, int Z_direct, int Z_random){
  
  vector[N] rho;
  matrix[N,Z_random] rhos;

  matrix[N,Z+n_regimes] X;
  matrix[N,n_regimes] X_reg;
  matrix[N,Z_direct+n_regimes] X_part;

  X_reg = calc_optima( a,  n_regimes,  n_lineages,  max_node_length,  nodes,  nodes_time,  t_end, t_beginning,  
  regime_time,  regimes_matrix);
  
  rho = (1 - (1 - exp(-a * T_term))./(a * T_term)); //For OU model
  rhos = rep_matrix(rho,Z_random);

  if(sum(random_cov)==0){
    X = append_col(X_reg,direct_cov);}
  else if(sum(direct_cov)==0){
    X = append_col(X_reg,random_cov .* rhos);}
  else{
    X_part = append_col(X_reg,direct_cov);
    X = append_col(X_part,random_cov .* rhos);}
  return(X);}
  
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code - Evolutionary
//////////////////////////////////////////////////////////////////////////////////////////////////////

matrix design_matrix_evol(int N, real a, vector T_term, matrix direct_cov, matrix random_cov, int n_regimes, int n_lineages,
int max_node_length, matrix nodes, matrix nodes_time, matrix t_end, matrix t_beginning,matrix regime_time, int[,] regimes_matrix, 
int Z, int Z_direct, int Z_random){
  matrix[N,Z+n_regimes] X;
  matrix[N,n_regimes] X_reg;
  matrix[N,Z_direct+n_regimes] X_part;
  
  X_reg = calc_optima( a,  n_regimes,  n_lineages,  max_node_length,  nodes,  nodes_time,  t_end, t_beginning,  regime_time,  regimes_matrix);

  
  if(sum(random_cov)==0){
    X = append_col(X_reg,direct_cov);}
  else if(sum(direct_cov)==0){
    X = append_col(X_reg,random_cov);}
  else{
    X_part = append_col(X_reg,direct_cov);
    X = append_col(X_part,random_cov);}
    
  return(X);}
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Vt function
//////////////////////////////////////////////////////////////////////////////////////////////////////
  matrix varcov_model(int N, matrix tij, matrix tja, matrix ta, matrix random_cov, real sigma2_y, real a, matrix x0,
  matrix sigma_squared_x, vector beta, vector T_term, int n_regimes, int Z, int Z_direct, int Z_random){
  
  vector[N] sigma2s;
  matrix[N,N] ti;
  matrix[N,N] term0;
  matrix[N,N] term1;
  matrix[N,N] term2;
  real s1;
  vector[Z_random] betasq;
  matrix[N,N] Vt;
  vector[Z_random] ones;
  
  ones = rep_vector(1,Z_random);

  if(sum(random_cov) != 0){
      if(Z_random>1){
        s1 = beta[2+Z_direct:1+Z_direct+Z_random]' * sigma_squared_x * ones;
        }else{
            s1 = sigma_squared_x[1,1] * (beta[1+Z_direct+Z_random]^2);
            }
      ti = rep_matrix(T_term,N);
      term0 = ((s1 + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) .* exp(-a * tij);
      term1 = (1 - exp(-a * ti)) ./ (a * ti); 
      term2 = exp(-a * tja) .* (1 - exp(-a * ti)) ./ (a * ti);

      Vt = term0 + s1 * (ta .* term1 .* (term1') - ((1 - exp(-a * ta)) ./ a) .* (term2 + (term2'))); //From Hansen et al. (2008)

    }
    else{ //Direct covariates only
      Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)
    }
  return(Vt);
  }


//////////////////////////////////////////////////////////////////////////////////////////////////////
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
data {
  //Extant data
  int N; //species number
  int Z; //total number of traits
  int Z_direct; //number of direct traits
  int Z_random; //number of adaptive traits
  int Z_me_direct;
  int Z_me_random;

  vector[N] Y; //y variable
  vector[N] mv_response;
  matrix[N,Z_direct] direct_cov;
  matrix[N,Z_me_direct] mv_direct_cov;
  matrix[N,Z_random] random_cov;
  matrix[N,Z_me_random] mv_random_cov;

  matrix[N,N] ta; //The following calculated in R based on the phylogeny
  vector[N] T_term; 
  matrix[N,N] tia;
  matrix[N,N] tja;
  matrix[N,N] tij;
  matrix[1,Z_random] brownian_mean;
  matrix[Z_random,Z_random] sigma_squared_x;
  real ols_intercept;
  vector[Z] ols_slope;
  
  //Fixed regimes
  int n_regimes;
  int n_lineages;
  int max_node_length;
  
  matrix[n_lineages,max_node_length] nodes;
  matrix[n_lineages,max_node_length] nodes_time;
  matrix[n_lineages,max_node_length] t_end;
  matrix[n_lineages,max_node_length] t_beginning;
  matrix[n_lineages,max_node_length] regime_time;
  int regimes_matrix[n_lineages,max_node_length];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////
parameters {
  real <lower = 0> hl;
  real <lower = 0> vy;
  real <lower = 0> sigma;
  vector[Z + n_regimes] beta; //OU beta
  vector[Z_random] beta_e; //OU beta
  matrix[N,Z_direct] direct_true;
  matrix[N,Z_random] random_true;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Transformed Parameter Block
//////////////////////////////////////////////////////////////////////////////////////////////////////
transformed parameters {



}
//////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////
model {
//Declare variables
  vector[N] mu;
  matrix[N,N] V;
  matrix[N,N] Vt;
  matrix[N,N] V_me;
  matrix[N,Z+n_regimes] X;
  matrix[N, N] L_V;
  real a;
  real sigma2_y;
  real rho;
  vector[N] Y_true;

//Priors
  hl ~ lognormal(log(0.4),1); //Tree length = 1 Ma
  vy ~ exponential(1);
  beta[1:n_regimes] ~ normal(ols_intercept,0.1); //Added for SBR1; 3 regimes prior
  beta[n_regimes+1:n_regimes+Z] ~ normal(ols_slope,0.1); //Added for SBR1; 3 regimes prior
  
//////////////////////////////////////////////////////////////////////////////////////////////////////
  a = log(2)/hl;
  sigma2_y = vy*(2*a);
  
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Regression - either constraint for direct cov or adaptive for random cov
//Set up X matrix
  if(Z_me_direct!= 0){// || (Z_me_random!= 0)){
    //for(i in 1:N){
      for(j in 1:Z_me_direct){
        direct_true[,j] ~ normal(0,1);
        direct_cov[,j] ~ normal(direct_true[,j],(mv_direct_cov[,j]));}
      //}
  }
  if(Z_me_random!= 0){// || (Z_me_random!= 0)){
    //for(i in 1:N){
      for(j in 1:Z_me_random){
        random_true[,j] ~ normal(0,1);
        random_cov[,j] ~ normal(random_true[,j],sqrt(mv_random_cov[,j]));}
      //}
  }
  X = design_matrix(N,  a,  T_term,  direct_true , random_true, n_regimes, n_lineages,max_node_length, nodes, nodes_time, t_end,
  t_beginning,regime_time,regimes_matrix, Z, Z_direct, Z_random);

  Vt = varcov_model(N,  tij,  tja,  ta,  random_cov, sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term, n_regimes, 
  Z, Z_direct, Z_random);

  L_V = cholesky_decompose(Vt);
  
  mu = X*beta;
  //Y_true ~ normal(mu,sqrt(mv_response));
  Y ~ multi_normal_cholesky(mu, L_V);
  
  //rho = (1 - (1 - exp(-a * T_term[1]))./(a * T_term[1])); //For OU model
  
  //for(i in 1:Z_random){
  //  beta_e[i] ~ normal(beta[n_regimes+Z_direct+i]* rho,0.2);
  //  }  

}
//////////////////////////////////////////////////////////////////////////////////////////////////////
generated quantities {
 
}

functions {
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Blouch Fixed Niches - Multilevel Model - Varying Intercepts Code, non-centered priors
//03.30.2023 - Added varying intercepts code, but centered priors
//09.04.2022 - Added code to allow for correlated predictors - Hansen et al. (2008)
//08.05.2022 - Revised for SBR1 to be use multivariate predictors, no need to use different random and direct datasets, and can use combo direct and response traits
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
  //print(dims(X));
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
//V_me function
//////////////////////////////////////////////////////////////////////////////////////////////////////
matrix varcov_measurement(int N, matrix ta, matrix direct_cov, matrix mv_direct_cov, matrix mv_random_cov, matrix sigma_squared_x, 
vector beta, int n_regimes, int Z, int Z_direct, int Z_random, real a, vector T_term){
  
  matrix[N,Z_direct] direct_rho2;
  matrix[N,Z_random] random_rho2;
  matrix[N,Z_direct+Z_random] rho2s;
  
  matrix[N,N] Vur[Z_random];
  matrix[N,N] Vxtr[Z_random];
  matrix[N,N] Vud[Z_direct];
  matrix[N,N] Vxtd[Z_direct];
  
  matrix[N,N] Vxt[Z];
  matrix[N,N] Vu[Z];
  matrix[N,N] Vx[Z];
  
  vector[N] Vu_given_x[Z];
  vector[N] beta2_Vu_given_x[Z];
  vector[N] beta2_Vu_given_x_sum;

  Vur = rep_array(rep_matrix(0,N,N),Z_random);
  Vxtr = rep_array(rep_matrix(0,N,N),Z_random);

  Vud = rep_array(rep_matrix(0,N,N),Z_direct);
  Vxtd = rep_array(rep_matrix(0,N,N),Z_direct);

  direct_rho2 = rep_matrix(1,N,Z_direct);
  random_rho2 = rep_matrix((1 - (1 - exp(-a * T_term))./(a * T_term))^2,Z_random); //For OU model
  
  rho2s = append_col(direct_rho2,random_rho2);

  if(sum(mv_random_cov) != 0){
    for (i in 1:Z_random){
      Vur[i] = diag_matrix(mv_random_cov[,i]);
      Vxtr[i] = ta * sigma_squared_x[1,i];
      }
    }
  if(sum(mv_direct_cov) != 0){
    for (i in 1:Z_direct){
      Vud[i] = diag_matrix(mv_direct_cov[,i]);
      Vxtd[i] = diag_matrix(rep_vector(variance(direct_cov[,i]),N))-Vud[i];
      }
  }
  
  Vxt = append_array(Vxtd,Vxtr);
  Vu = append_array(Vud,Vur);
  
  for (i in 1:Z){ //Adding variance matrices together for individual traits 
    Vx[i] = Vxt[i] + Vu[i];
    Vu_given_x[i] = diagonal(Vu[i] - Vu[i] * inverse(Vx[i]) * Vu[i]);
    beta2_Vu_given_x[i] = to_vector(Vu_given_x[i] * square(beta[i+n_regimes]) .* rho2s[,i]); //Changed to element-wise, appears to be same as R
    }  

 
  
  for (i in 1:N){
        beta2_Vu_given_x_sum[i] = sum(beta2_Vu_given_x[,i]);
      }

  return(diag_matrix(beta2_Vu_given_x_sum));
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
  
  int n_regimes; //Fixed regimes
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
  real <lower = 0, upper = 3> hl;
  real <lower = 0> vy;
  vector[Z + n_regimes] beta; //OU beta
  vector[Z_random] beta_e; //OU beta
  real beta_bar; //average prior for regimes
  real sigma; //standard deviation for regimes
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

//Priors
  hl ~ lognormal(-0.5,1.5); //Tree length = 1 Ma
  vy ~ exponential(1);
  beta_bar ~ normal(0,1);
  sigma ~ exponential(1);
  
  beta[1:n_regimes] ~ normal(beta_bar,sigma); //Varying intercepts model
  beta[n_regimes+1:n_regimes+Z] ~ normal(ols_slope,0.5); //Static slope

  a = log(2)/hl;
  sigma2_y = vy*(2*a);
  
  //Regression - either constraint for direct cov or adaptive for random cov
  X = design_matrix(N,  a,  T_term, direct_cov,random_cov, n_regimes, n_lineages,max_node_length, nodes, nodes_time, t_end,
  t_beginning,regime_time,regimes_matrix, Z, Z_direct, Z_random);//Set up X matrix

  Vt = varcov_model(N,  tij,  tja,  ta,  random_cov, sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term, n_regimes, 
  Z, Z_direct, Z_random);//Set up V matix
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    V_me = varcov_measurement(N, ta, direct_cov, mv_direct_cov,mv_random_cov, sigma_squared_x, beta, n_regimes, Z, Z_direct, Z_random, a, T_term);
      if(sum(mv_response)!= 0){
          V = Vt + V_me + diag_matrix(mv_response);
      }else{
          V = Vt + V_me; //Hansen & Bartoszek (2012) - Eq 10
          }
  }else{
    V_me = rep_matrix(0,N,N);
    if(sum(mv_response)!= 0){
       V = Vt + diag_matrix(mv_response);
    }else{V = Vt;}
  }

  V = Vt;
  L_V = cholesky_decompose(V);
  
  //OU with random covariates
  mu = X*beta;
  Y ~ multi_normal_cholesky(mu , L_V);

  rho = (1 - (1 - exp(-a * T_term[1]))./(a * T_term[1])); //For OU model
  for(i in 1:Z_random){
    beta_e[i] ~ normal(beta[n_regimes+Z_direct+i]* rho,0.2);
    }  
  //target += log_lik;  // automatically sums

}
//////////////////////////////////////////////////////////////////////////////////////////////////////
generated quantities {
  vector[N] mu;
  matrix[N,N] V;
  matrix[N,N] Vt;
  matrix[N,N] V_me;
  matrix[N, N] L_V;
  matrix[N,Z+n_regimes] X;
  real a;
  real sigma2_y;
  //vector[N] log_lik;
  real log_lik;
  
  a = log(2)/hl;
  sigma2_y = vy*(2*a);

  X = design_matrix(N,  a,  T_term, direct_cov,random_cov, n_regimes, n_lineages,max_node_length, nodes, nodes_time, t_end,
  t_beginning,regime_time,regimes_matrix, Z, Z_direct, Z_random);//Set up X matrix

  Vt = varcov_model(N,  tij,  tja,  ta,  random_cov, sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term, n_regimes, 
  Z, Z_direct, Z_random);//Set up V matix
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    V_me = varcov_measurement(N, ta, direct_cov, mv_direct_cov,mv_random_cov, sigma_squared_x, beta, n_regimes, Z, Z_direct, Z_random, a, T_term);
      if(sum(mv_response)!= 0){
          V = Vt + V_me + diag_matrix(mv_response);
      }else{
          V = Vt + V_me; //Hansen & Bartoszek (2012) - Eq 10
          }
  }else{
    V_me = rep_matrix(0,N,N);
    if(sum(mv_response)!= 0){
       V = Vt + diag_matrix(mv_response);
    }else{V = Vt;}
  }

  V = Vt;
  L_V = cholesky_decompose(V);
  mu = X*beta;
  log_lik = multi_normal_lpdf(Y | mu, V);

}

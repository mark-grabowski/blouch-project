functions {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Blouch Fixed Niches Code
//08.05.2022 - Revised for SBR1 to be use multivariate predictors, no need to use different random and direct datasets, and can use combo direct and response traits
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrix calc_optima(real a, int n_regimes, int n_lineages, int max_node_length, matrix nodes, matrix nodes_time, 
matrix t_end,matrix t_beginning, matrix regime_time, int[,] regimes_matrix){
  matrix[n_lineages,max_node_length+1] weights_seg_matrix = rep_matrix(0.0,n_lineages,max_node_length+1);
  matrix[n_lineages,n_regimes] optima_matrix = rep_matrix(0.0,n_lineages,n_regimes);
  real regimes_matrix_sep[n_lineages, max_node_length+1,n_regimes] = rep_array(0.0,n_lineages, max_node_length+1,n_regimes);
  
  //Each niche will get its own weighting based on length of time spend in it across each lineage
  for(i in 1:n_lineages){
    for(j in 1:max_node_length){
      if(t_end[i,j] != 0){
        weights_seg_matrix[i,j] = exp(-a * t_beginning[i,j]) - exp(-a * t_end[i,j]);
        weights_seg_matrix[i,j+1] = exp(-a * nodes_time[i,1]);
      }
    }
  }

  for(i in 1:n_lineages){
    for(j in 1:max_node_length){
      int z = regimes_matrix[i,j];
      if (z!=0){
        regimes_matrix_sep[i,j,z] = weights_seg_matrix[i,j];
      }
    }
  }

  for(z in 1:n_regimes){
    for(i in 1:n_lineages){
      for(j in 1:max_node_length){
        optima_matrix[i,z] = sum(regimes_matrix_sep[i,,z]);
      }
    }
  }

return(optima_matrix);
}
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code - Adaptive
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code - Evolutionary
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrix design_matrix_evol(int N, real a, vector T_term, matrix direct_cov, matrix random_cov, int n_regimes, int n_lineages,
int max_node_length, matrix nodes, matrix nodes_time, matrix t_end, matrix t_beginning,matrix regime_time, int[,] regimes_matrix, 
int Z, int Z_direct, int Z_random){
  vector[N] rho;
  matrix[N,Z_random] rhos;
  matrix[N,Z+n_regimes] X;
  matrix[N,n_regimes] X_reg;
  matrix[N,Z_direct+n_regimes] X_part;
  
  X_reg = calc_optima( a,  n_regimes,  n_lineages,  max_node_length,  nodes,  nodes_time,  t_end, t_beginning,  regime_time,  regimes_matrix);
  rho=rep_vector(1,N); //For Evolutionary regression
  rhos = rep_matrix(rho,Z_random);

  
  if(sum(random_cov)==0){
    X = append_col(X_reg,direct_cov);}
  else if(sum(direct_cov)==0){
    X = append_col(X_reg,random_cov .* rhos);}
  else{
    X_part = append_col(X_reg,direct_cov);
    X = append_col(X_part,random_cov .* rhos);}
    
  return(X);}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Vt function
///////////////////////////////////////////////////////////////////////////////////////////////////////////
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

  for (i in 1:Z_random){
    betasq[i] = beta[i+Z_direct+n_regimes]^2; //Only betas that are for random covariates
    }

    //Random + Direct covariates, or just random
  if(sum(random_cov) != 0){
      s1 = sum(sigma_squared_x * betasq); 
      ti = rep_matrix(T_term,N);

      term0 = ((s1 + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) .* exp(-a * tij);
      term1 = (1 - exp(-a * ti)) ./ (a * ti); 
      term2 = exp(-a * tja) .* (1 - exp(-a * ti)) ./ (a * ti);

      Vt = term0 + s1 * (ta .* term1 .* (term1') - ((1 - exp(-a * ta)) ./ a) .* (term2 + (term2')));
      //Vt = term0 + s1 * ta .* (term1 .* (term1') - ((1 - exp(-a * ta)) ./ (a * ta)) .* (term2 + (term2')));    

    }
    else{ //Direct covariates only
      Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij));
      //Vt = (exp(-a * tij)) .* ((sigma2_y /( 2 * a)) * (1 - exp(-2 * a * ta)));
    }
  return(Vt);
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//V_me function
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

  //for (i in 1:Z){ //Adding variance matrices together for individual traits Original SB Code
  //Vxt[i] = Vxtd[i]+Vxtr[i];
    //Vu[i] = Vud[i]+Vur[i];rs
    //Vx[i] = Vxt[i] + Vu[i];
    //Vu_given_x[i] = Vu[i] - Vu[i] * inverse(Vx[i]) * Vu[i];
    //beta2_Vu_given_x[i] = Vu_given_x[i] * square(beta[i+n_regimes]);
    //}
  
  Vxt = append_array(Vxtd,Vxtr);
  Vu = append_array(Vud,Vur);
  
  for (i in 1:Z){ //Adding variance matrices together for individual traits 
    Vx[i] = Vxt[i] + Vu[i];
    Vu_given_x[i] = diagonal(Vu[i] - Vu[i] * inverse(Vx[i]) * Vu[i]);
    //print("Rep");
    //print(i);
    //print(length(square(beta[i+n_regimes])));
    //print(dims(Vu_given_x[i] * square(beta[i+n_regimes]) .* rho2s[,i])); //Changed to element-wise, appears to be same as R
    beta2_Vu_given_x[i] = to_vector(Vu_given_x[i] * square(beta[i+n_regimes]) .* rho2s[,i]); //Changed to element-wise, appears to be same as R
    }  

    
  //print(Vu[1,1,1]);
  //print(Vu[2,1,1]);

  //print(Vx[1,1,1]);
  //print(Vx[2,1,1]); 

  //print(Vu_given_x[1,1,1]);
  //print(Vu_given_x[2,1,1]); //Problem - "nan"
  //print("Rep");

  //print(beta2_Vu_given_x);

  //for (j in 1:N){
  //  for (k in 1:N){
  //      beta2_Vu_given_x_sum[j,k] = sum(beta2_Vu_given_x[,j,k]);
  //    }
  //  }
  
  for (i in 1:N){
        //print(beta2_Vu_given_x[,i]);
        beta2_Vu_given_x_sum[i] = sum(beta2_Vu_given_x[,i]);
      }
  //print(beta2_Vu_given_x_sum);

  return(diag_matrix(beta2_Vu_given_x_sum));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
data {
  //Extant data
  int N; //species number
  int Z; //total number of traits
  int Z_direct; //number of direct traits
  int Z_random; //number of adaptive traits
  vector[N] Y; //y variable
  vector[N] mv_response;
  matrix[N,Z_direct] direct_cov;
  matrix[N,Z_direct] mv_direct_cov;
  matrix[N,Z_random] random_cov;
  matrix[N,Z_random] mv_random_cov;
  
  matrix[N,N] ta; //The following calculated in R based on the phylogeny
  vector[N] T_term; 
  matrix[N,N] tia;
  matrix[N,N] tja;
  matrix[N,N] tij;
  matrix[1,Z_random] brownian_mean;
  matrix[1,Z_random] sigma_squared_x;
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
parameters {
  real <lower = 0> a;
  real <lower = 0, upper = variance(Y)*4> sigma2_y; //Added to limit the variance based on Kjetil's suggestion
  vector[Z + n_regimes] beta; //OU beta
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Transformed Parameter Block
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
transformed parameters {
  matrix[N,Z+n_regimes] X_evol;
  vector[Z+n_regimes] beta_evol;
  matrix[N,N] V_ev;
  matrix[N,N] V_me_ev;
  matrix[N,N] Vt_ev;

  Vt_ev = varcov_model(N,  tij,  tja,  ta,  random_cov,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,T_term, n_regimes, Z, Z_direct, Z_random);

  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    V_me_ev = varcov_measurement(N, ta, direct_cov, mv_direct_cov, mv_random_cov, sigma_squared_x, beta, n_regimes, Z, Z_direct, Z_random, a, T_term);
    V_ev = Vt_ev + V_me_ev + diag_matrix(mv_response);
    }
  else{
    V_me_ev = rep_matrix(0,N,N);
    V_ev = Vt_ev;
    }
  //Calculate evolutionary regression slope
  X_evol = design_matrix_evol(N,  a,  T_term, direct_cov,random_cov, n_regimes, n_lineages,max_node_length, nodes, nodes_time, 
  t_end, t_beginning,regime_time,regimes_matrix, Z, Z_direct, Z_random);
  beta_evol = inverse(X_evol'*inverse(V_ev)*X_evol)*(X_evol'*inverse(V_ev)*Y); //Hansen et al. 2008
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
model {
//Declare variables
  vector[N] mu;
  matrix[N,N] V;
  matrix[N,N] Vt;
  matrix[N,N] V_me;
  matrix[N,Z+n_regimes] X;
  matrix[N, N] L_V;

//Priors
  a ~ lognormal(1.0,1.0); //a = log(2)/half-life
  beta[1:n_regimes] ~ normal(ols_intercept,0.5); //Changed from 0.5 for Original Sub
  beta[n_regimes+1:n_regimes+Z] ~ normal(ols_slope,0.4); //Changed from 0.4 for Original Sub
  
//////////////////////////////////////////////////////////////////////////////////////////////////////
  //hl = log(2)/a;
  //a = log(2)/hl;
  //sigma2_y = vy*(2*a);

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Regression - either constraint for direct cov or adaptive for random cov

//Set up X matrix
  X = design_matrix(N,  a,  T_term, direct_cov,random_cov, n_regimes, n_lineages,max_node_length, nodes, nodes_time, t_end,
  t_beginning,regime_time,regimes_matrix, Z, Z_direct, Z_random);

  //Set up V matix
  Vt = varcov_model(N,  tij,  tja,  ta,  random_cov, sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term, n_regimes, 
  Z, Z_direct, Z_random);
  V_me = varcov_measurement(N, ta, direct_cov, mv_direct_cov,mv_random_cov, sigma_squared_x, beta, n_regimes, Z, Z_direct, Z_random, a, T_term);
  V = Vt + V_me + diag_matrix(mv_response);
  V = Vt;
  //print(V[1,1]);
  //print(dims(X));
  //print(dims(beta));
  //print(dims(Vt));
  //print(dims(V_me));
  L_V = cholesky_decompose(V);
  //OU with random covariates
  mu = X*beta;
  Y ~ multi_normal_cholesky(mu , L_V);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
generated quantities {
  real <lower = 0> vy;
  //real <lower = 0> sigma2_y;
  //real <lower = 0> a;
  real <lower = 0> hl;
  vector[N] pred_mean;
  real grand_mean;
  real sst;
  real sse;
  real r_squared;
  matrix[N,N] V_final;
  matrix[N,N] V_me_final;
  matrix[N,N] Vt_final;
  matrix[N,Z+n_regimes] X;


  //a = log(2)/hl;
  hl = log(2)/a;
  vy = sigma2_y/(2*a);
  //sigma2_y = vy*(2*a);
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Set up X matrix
  X = design_matrix(N, a,  T_term, direct_cov, random_cov, n_regimes, n_lineages,max_node_length, nodes, nodes_time,
  t_end, t_beginning, regime_time, regimes_matrix, Z, Z_direct, Z_random);

//Calculate V matrix
  Vt_final = varcov_model(N,  tij,  tja,  ta,  random_cov, sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  
  T_term, n_regimes, Z, Z_direct, Z_random);
  if((sum(mv_direct_cov)!= 0) || (sum(mv_random_cov)!= 0)){
    V_me_final = varcov_measurement(N, ta, direct_cov, mv_direct_cov, mv_random_cov, sigma_squared_x, beta,n_regimes, Z, Z_direct, Z_random, a, T_term);
    V_final = Vt_final + V_me_final + diag_matrix(mv_response);
    }
  else{
   V_me_final = rep_matrix(0,N,N);
   V_final = Vt_final;
  }
  //V_final = Vt_final;

  //Calculate r2 based on constraint or adaptive regression
  pred_mean = X*beta;
  grand_mean = ((rep_vector(1,N))' * inverse(V_final) * Y) / sum(inverse(V_final));
  sst = ((Y - grand_mean)' * inverse(V_final) * (Y - grand_mean));
  sse = ((Y - pred_mean)' * inverse(V_final) * (Y - pred_mean));
  r_squared = (sst - sse) / sst;

}

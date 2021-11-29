functions {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Blouch Intercept only (OU1) Code
//02.11.2021
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
matrix calc_optima(real a, int n_regimes, int n_lineages, int max_node_length, matrix nodes, matrix nodes_time, matrix t_end,matrix t_beginning, matrix regime_time, int[,] regimes_matrix){

  //matrix[n_lineages,max_node_length+1] weights_seg_matrix;
  matrix[n_lineages,max_node_length+1] weights_seg_matrix = rep_matrix(0.0,n_lineages,max_node_length+1);
  matrix[n_lineages,n_regimes] optima_matrix = rep_matrix(0.0,n_lineages,n_regimes);
  real regimes_matrix_sep[n_lineages, max_node_length+1,n_regimes] = rep_array(0.0,n_lineages, max_node_length+1,n_regimes);
  
  //int z;
  //Each niche will get its own weighting based on length of time spend in it across each lineage
  //print(regimes_matrix);
  for(i in 1:n_lineages){
    for(j in 1:max_node_length){
      if(t_end[i,j] != 0){
        weights_seg_matrix[i,j] = exp(-a * t_beginning[i,j]) - exp(-a * t_end[i,j]);
        weights_seg_matrix[i,j+1] = exp(-a * nodes_time[i,1]);
      }
    }
  }

  //print("weights_seg_matrix",weights_seg_matrix);
  //print("regimes_matrix",regimes_matrix);

  for(i in 1:n_lineages){
    for(j in 1:max_node_length){
      int z = regimes_matrix[i,j];
      if (z!=0){
        regimes_matrix_sep[i,j,z] = weights_seg_matrix[i,j];
      }
    }
  }
  ///print(regimes_matrix_sep[,3,1]);

  for(z in 1:n_regimes){
    for(i in 1:n_lineages){
      for(j in 1:max_node_length){
        //print(regimes_matrix_sep[i,,z]);
        optima_matrix[i,z] = sum(regimes_matrix_sep[i,,z]);
      }
    }
  }

return(optima_matrix);
}
  
  

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Design Matrix Code
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  matrix design_matrix(int N, int evol, real a, vector T_term, int Z,int n_regimes, int n_lineages, int max_node_length, matrix nodes, matrix nodes_time, matrix t_end, matrix t_beginning,
  matrix regime_time, int[,] regimes_matrix){
  matrix[N,n_regimes] X_reg;

  X_reg = calc_optima( a,  n_regimes,  n_lineages,  max_node_length,  nodes,  nodes_time,  t_end, t_beginning,  regime_time,  regimes_matrix);
    
  return(X_reg);
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Vt function
///////////////////////////////////////////////////////////////////////////////////////////////////////////
  matrix varcov_model(int N, matrix tij, matrix tja, matrix ta,int Z, real sigma2_y, real a, vector beta1, vector T_term,int n_regimes){
  matrix[N,N] Vt;
  
  Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij));
  
  return(Vt);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
data {
  //Extant data
  int N; //species number
  int Z; //number of traits
  vector[N] Y; //y variable
  vector[N] mv_response;
  matrix[N,N] ta; //The following calculated in R based on the phylogeny
  vector[N] T_term; 
  matrix[N,N] tia;
  matrix[N,N] tja;
  matrix[N,N] tij;

  
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
  //real <lower = 0.1386294, upper = 100> a; //Lower based on a little more than a for a hl of 3*tree height - 6.93 * tree height
  //real <lower = 0> hl;
  real <lower = 0, upper = 1000> hl; //
  //real <lower = 0,upper = variance(Y)*2> sigma2_y; //Added to limit the variance based on Kjetil's suggestion
  //real <lower = 0,upper = variance(Y)*2> vy;
  real <lower = 0> vy;  
  vector[Z + n_regimes] beta; //OU beta
  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
transformed parameters {

  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
model {
//Declare variables
  real a;
  //real hl;
  real sigma2_y;
  vector[N] mu;
  matrix[N,N] V;
  //matrix[N,N] Vt;
  //matrix[N,N] V_me;
  matrix[N,Z+n_regimes] X;
  matrix[N, N] L_V;

//Priors
  //a ~ lognormal(-0.4,0.5);
  //hl ~ lognormal(1.75,1.25); //Tree length = 1 Ma
  //hl ~ normal(-1.5,1.0); //Tree length = 1 Ma

  //hl ~ cauchy(0,0.1); //Tree length = 1 Ma
  //vy ~ exponential(1.0);
  //vy ~ lognormal(log(0.1),0.25);

  beta ~ normal(0,2);
  
//////////////////////////////////////////////////////////////////////////////////////////////////////
  //hl = log(2)/a;
  a = log(2)/hl;
  sigma2_y = vy*(2*a);

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Regression - either constraint for direct cov or adaptive for random cov

//Set up X matrix
  X = design_matrix( N,  1,  a,  T_term, Z, n_regimes, n_lineages,
  max_node_length, nodes, nodes_time, t_end, t_beginning,regime_time,regimes_matrix);

//Set up V matix
  V = varcov_model(N,  tij,  tja,  ta, Z,  sigma2_y,  a,  beta,  T_term, n_regimes);
  //V_me = varcov_measurement(N, Z, ta, direct_cov, mv_direct_cov,mv_random_cov, sigma_squared_x, beta);
  //V = Vt + V_me + diag_matrix(mv_response);
  //V=Vt;
  L_V = cholesky_decompose(V);
  
//OU with random covariates
  mu = X*beta;
  Y ~ multi_normal_cholesky(mu , L_V);
  //print(mu);

////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

generated quantities {
  //real <lower = 0> vy;
  real <lower = 0> sigma2_y;
  real <lower = 0> a;
  //real <lower = 0> hl;
  vector[N] pred_mean;
  real grand_mean;
  real sst;
  real sse;
  real r_squared;
  matrix[N,N] V_final;
  //matrix[N,N] V_me_final;
  matrix[N,N] Vt_final;
  matrix[N,Z+n_regimes] X_evol;
  vector[Z+n_regimes] beta_evol;

  a = log(2)/hl;
  //hl = log(2)/a;
  //vy = sigma2_y/(2*a);
  sigma2_y = vy*(2*a);
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Calculate V matrix
  Vt_final = varcov_model(N,  tij,  tja,  ta, Z,  sigma2_y,  a,  beta,  T_term, n_regimes);
  //V_me_final = varcov_measurement(N, Z, ta, direct_cov, mv_direct_cov, mv_random_cov, sigma_squared_x, beta);
  //V_final = Vt_final + V_me_final + diag_matrix(mv_response);
  V_final=Vt_final;
  
//Calculate evolutionary regression slope
  X_evol = design_matrix( N,  1,  a,  T_term, Z, n_regimes, n_lineages,
  max_node_length, nodes, nodes_time, t_end, t_beginning,regime_time,regimes_matrix);

  beta_evol = inverse(X_evol'*inverse(V_final)*X_evol)*X_evol'*inverse(V_final)*Y; //Hansen et al. 2008
//Calculate r2 based on constraint or adaptive regression

  pred_mean = (X_evol*beta_evol);
  grand_mean = ((rep_vector(1,N))' * (V_final') * Y) / sum(inverse(V_final));
  sst = ((Y - grand_mean)' * inverse(V_final) * (Y - grand_mean));
  sse = ((Y - pred_mean)' * inverse(V_final) * (Y - pred_mean));
  r_squared = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////



}

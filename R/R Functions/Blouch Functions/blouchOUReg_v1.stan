functions {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Blouch Fixed Niches Code
//02.09.2021
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
  matrix design_matrix(int N, int evol, real a, vector T_term, matrix direct_cov, matrix random_cov, int Z,int n_regimes, int n_lineages, int max_node_length, matrix nodes, matrix nodes_time, matrix t_end, matrix t_beginning,
  matrix regime_time, int[,] regimes_matrix){
  matrix[N,Z] rho;
  matrix[N,Z+n_regimes] X;
  matrix[N,n_regimes] X_reg;
  matrix[n_lineages,n_regimes] optima_matrix = rep_matrix(0.0,n_lineages,n_regimes);

    if(n_regimes != 0){
      X_reg = calc_optima( a,  n_regimes,  n_lineages,  max_node_length,  nodes,  nodes_time,  t_end, t_beginning,  regime_time,  regimes_matrix);
    }
    if(evol==0) {
      rho = to_matrix(1 - (1 - exp(-a * T_term))./(a * T_term)); //For OU model
    }
    else{
      rho=to_matrix(rep_vector(1.0,N)); //For Evolutionary regression
    }

    if(sum(random_cov)==0){
     if(n_regimes != 0){
      X = append_col(X_reg,direct_cov);
     }else{
      X = direct_cov;
      }
    }
    if(sum(direct_cov)==0){
      if(n_regimes != 0){
        X = append_col(X_reg,random_cov .* rho);
      }else{
        X = random_cov .* rho;
      }
    }
    return(X);
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Vt function
///////////////////////////////////////////////////////////////////////////////////////////////////////////
  matrix varcov_model(int N, matrix tij, matrix tja, matrix ta, matrix random_cov,int Z, real sigma2_y, real a, matrix x0, matrix sigma_squared_x, vector beta1, vector T_term,int n_regimes){
  
  vector[N] sigma2s;
  matrix[N,N] ti;
  matrix[N,N] term0;
  matrix[N,N] term1;
  matrix[N,N] term2;
  real s1;
  vector[Z] beta1sq;
  matrix[N,N] Vt;

  for (i in 1:Z){
    beta1sq[i] = beta1[i+n_regimes]; //Only betas that are for random covariates
    beta1sq[i] = beta1sq[i]^2; //Square betas
    }

    //Random + Direct covariates, or just random
    if(sum(random_cov) != 0){
      s1 = sum(sigma_squared_x * beta1sq); 
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
matrix varcov_measurement(int N, int Z, matrix ta, matrix direct_cov, matrix mv_direct_cov, matrix mv_random_cov, matrix sigma_squared_x, vector beta){
  //matrix[N,N] Vxt[Z];
  matrix[N,N] Vxtr;
  matrix[N,N] Vxtd;
  matrix[N,N] Vxt;
  matrix[N,N] Vx;
  matrix[N,N] Vu_given_x;
  matrix[N,N] beta2_Vu_given_x;
  matrix[N,N] Vu;
  matrix[N,N] Vur;
  matrix[N,N] Vud;
  matrix[N,Z] sigma_squared_x_rep;
  matrix[N,Z] P;
  Vur = rep_matrix(0,N,N);
  Vud = rep_matrix(0,N,N);
  Vxtd = rep_matrix(0,N,N);
  Vxtr = rep_matrix(0,N,N);

  if(sum(mv_random_cov) != 0){
    Vur = diag_matrix(to_vector(mv_random_cov));
    Vxtr = ta * sigma_squared_x[1,1];
    }
    
  if(sum(mv_direct_cov) != 0){
    Vud = diag_matrix(to_vector(mv_direct_cov[,1])); //Puts 
    Vxtd = diag_matrix(rep_vector(variance(direct_cov[,1]),N))-Vud;
    }

  Vxt = Vxtd+Vxtr;
  Vu = Vud+Vur;
  Vx = Vxt + Vu;

//Vu_given_x
  Vu_given_x = Vu - Vu * inverse(Vx) * Vu;

//beta2_Vu_given_x
  beta2_Vu_given_x = Vu_given_x * square(beta[1]);
    
  return(beta2_Vu_given_x);
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
data {
  //Extant data
  int N; //species number
  int Z; //number of traits
  vector[N] Y; //y variable
  vector[N] mv_response;
  matrix[N,Z] direct_cov;
  matrix[N,Z] mv_direct_cov;
  matrix[N,Z] random_cov;
  matrix[N,Z] mv_random_cov;
  matrix[N,N] ta; //The following calculated in R based on the phylogeny
  vector[N] T_term; 
  matrix[N,N] tia;
  matrix[N,N] tja;
  matrix[N,N] tij;
  matrix[1,Z] brownian_mean;
  matrix[1,Z] sigma_squared_x;
  real ols_intercept;
  real ols_slope;
  
  
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
  real <lower = 0.1386294, upper = 100> a; //Lower based on a little more than a for a hl of 3*tree height - 6.93 * tree height
  //real <lower = 0.001,upper=5.0> hl;
  //real <lower = 0> hl; //Upper is 4*tree height
  //real <lower = 0,upper = variance(Y)*2> sigma2_y; //Added to limit the variance based on Kjetil's suggestion
  //real <lower = 0,upper = variance(Y)*2> vy;
  real <lower = 0, upper = variance(Y)*2> sigma2_y;  
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
  //real a;
  //real hl;
  //real sigma2_y;
  vector[N] mu;
  matrix[N,N] V;
  //matrix[N,N] Vt;
  //matrix[N,N] V_me;
  matrix[N,Z+n_regimes] X;
  matrix[N, N] L_V;

//Priors
  a ~ lognormal(-0.5,0.75);
  //hl ~ lognormal(-1.5,1.35); //Tree length = 1 Ma
  //sigma2_y~exponential(1.0);
  //sigma2_y~cauchy(0,1.0);
  //vy ~ exponential(1.0);
  beta[1:n_regimes] ~ normal(ols_intercept,1.0);
  beta[n_regimes+1] ~ normal(ols_slope,1.0);
  
//////////////////////////////////////////////////////////////////////////////////////////////////////
  //hl = log(2)/a;
  //a = log(2)/hl;
  //sigma2_y = vy*(2*a);

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Regression - either constraint for direct cov or adaptive for random cov

//Set up X matrix
  X = design_matrix( N,  1,  a,  T_term, direct_cov,random_cov, Z, n_regimes, n_lineages,
  max_node_length, nodes, nodes_time, t_end, t_beginning,regime_time,regimes_matrix);

//Set up V matix
  V = varcov_model(N,  tij,  tja,  ta,  random_cov, Z,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term, n_regimes);
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
  //matrix[N,N] V_me_final;
  matrix[N,N] Vt_final;
  matrix[N,Z+n_regimes] X_evol;
  vector[Z+n_regimes] beta_evol;

  //a = log(2)/hl;
  hl = log(2)/a;
  vy = sigma2_y/(2*a);
  //sigma2_y = vy*(2*a);
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Calculate V matrix
  Vt_final = varcov_model(N,  tij,  tja,  ta,  random_cov, Z,  sigma2_y,  a, brownian_mean,  sigma_squared_x,  beta,  T_term, n_regimes);
  //V_me_final = varcov_measurement(N, Z, ta, direct_cov, mv_direct_cov, mv_random_cov, sigma_squared_x, beta);
  //V_final = Vt_final + V_me_final + diag_matrix(mv_response);
  V_final=Vt_final;
  
//Calculate evolutionary regression slope
  X_evol = design_matrix( N,  1,  a,  T_term, direct_cov,random_cov, Z, n_regimes, n_lineages,
  max_node_length, nodes, nodes_time, t_end, t_beginning,regime_time,regimes_matrix);

  beta_evol = inverse(X_evol'*inverse(V_final)*X_evol)*X_evol'*inverse(V_final)*Y; //Hansen et al. 2008
//Calculate r2 based on constraint or adaptive regression

  pred_mean = (X_evol*beta_evol);
  grand_mean = ((rep_vector(1,N))' * (V_final') * Y) / sum(inverse(V_final));
  sst = ((Y - grand_mean)' * inverse(V_final) * (Y - grand_mean));
  sse = ((Y - pred_mean)' * inverse(V_final) * (Y - pred_mean));
  r_squared = (sst - sse) / sst;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////



}

//Blouch OU model reprogrammed - 2023
//Combination of regime model with direct effect predictors
//Regime model - for single regime painting and SIMMAPS
//Model accounts for measurement error on X and/or Y variables
//See below for which lines to comment/uncomment based on whether measurement error is present

functions {
  int num_matches(vector x, real y) { //Thanks to Stan Admin Jonah -https://discourse.mc-stan.org/t/how-to-find-the-location-of-a-value-in-a-vector/19768/2
  int n = 0;
  for (i in 1:rows(x))
    if (x[i] == y)
      n += 1;
  return(n);
  }
  
  int[] which_equal(vector x, real y) {
    int match_positions[num_matches(x, y)];
    int pos = 1;
    //for (i in 1:size(x)) {
    for (i in 1:(dims(x)[1])) {  
      if (x[i] == y) {
        match_positions[pos] = i;
        pos += 1;
        }
      }
    return(match_positions);
  }

  vector weight_segments(real a, vector t_beginning, vector t_end, real time, int nodes){//Individual lineage, calculate weights per segment
    vector[nodes] weights = append_row(exp(-a * t_beginning) - exp(-a * t_end),exp(-a * time));
    return(weights);
  }
  
  row_vector weights_regimes(int n_reg, real a, vector t_beginning, vector t_end, real time, vector reg_match, int nodes){//
    //Individual lineage, calculate weights for regimes on each segement
    vector[nodes] weight_seg = weight_segments(a, t_beginning[1:(nodes-1)], t_end[1:(nodes-1)], time, nodes);
    //print(weight_seg);
    vector[n_reg] reg_weights = rep_vector(0,n_reg);
    for(i in 1:n_reg){//reg_match should have values 1,2,3 denoting different regimes
      int ids[num_matches(reg_match, i)] = which_equal(reg_match, i); //Returns indixes of matching regimes in weight_segments vector
      //print(ids);
      //print(weight_seg[ids]);
      reg_weights[i] = sum(weight_seg[ids]);
      //print(reg_weights[i]);
      }
    return(reg_weights');
  }
  
  matrix calc_optima_matrix(int N, int n_reg, real a, matrix t_beginning, matrix t_end, matrix times, matrix reg_match, int[] nodes){
    matrix[N,n_reg] optima_matrix = rep_matrix(0,N,n_reg);
    for(i in 1:N){ //For each tip/lineage, figure out weighting of regimes
      optima_matrix[i,] = weights_regimes(n_reg, a, t_beginning[i,]', t_end[i,]', times[i,1], reg_match[i,]', nodes[i]);
      //print(i);
      //print(optima_matrix[i,]);
      }
    return(optima_matrix);
  }
  matrix calc_direct_V( real a,real sigma2_y,matrix ta, matrix tij) {
        int N = dims(ta)[1];
        matrix[N, N] Vt;
        Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)
        return Vt;
  }
}

data {
  int N; //Numnber of tips
  int n_reg; //Number of regimes
  int Z_direct;
  int Z_X_error;
  int max_node_num; //Max number of nodes in lineage
  vector[N] Y_obs; //Y observed
  matrix[N,Z_direct] X_obs;
  vector[N] Y_error;
  matrix[N,Z_X_error] X_error;
  matrix[N,N] ta; //Time from tip to ancestor
  matrix[N,N] tij;
  matrix[N,N] tja;
  vector[N] T_term;
  matrix[N, max_node_num] t_beginning; //Matrix of times for beginning of segments to node
  matrix[N, max_node_num] t_end; //Matrix of times for end of segments to 
  matrix[N, max_node_num] times; //Matrix of root to node times
  matrix[N, max_node_num] reg_match; //Matrix of 1,2,3 denoting each regime for each node in a lineage. 0 if no node
  int nodes[N]; //Vector of number of nodes per lineage
}

parameters {
  real<lower=0> hl;
  real <lower=0> vy;
  vector[n_reg] optima;
  vector[Z_direct] beta;
  vector[N] Y;
  matrix[N,Z_direct] X;
}

model {
  matrix[N,N] V;
  vector[N] mu;
  matrix[N,N] L_v;
  matrix[N,n_reg+Z_direct] dmX;
  real a = log(2)/hl;
  real sigma2_y = vy*(2*(log(2)/hl));
  matrix[N,n_reg] optima_matrix;
  vector[n_reg+Z_direct] optima_beta = append_row(optima,beta);
  //hl ~ lognormal(log(0.25),0.25);
  target += lognormal_lpdf(hl|log(0.25),0.75);
  //vy ~ exponential(5);
  target += exponential_lpdf(vy|10);
  //optima ~ normal(-1.179507,0.75);
  //beta ~ normal(6.304451,1.5);
  target += normal_lpdf(optima|-1.179507,0.75);
  target += normal_lpdf(beta|6.304451,1.5);

  for(i in 1:(Z_direct)){//Given measurement error in X variable, uncomment this nested statement
    //X[,i] ~ normal(0,1);
    target += normal_lpdf(X[,i]|0,1);
    //X_obs[,i] ~ normal(X[,i], X_error[,i]);
    target += normal_lpdf(X_obs[,i]|X[,i],X_error[,i]);
  }
  optima_matrix = calc_optima_matrix(N, n_reg, a, t_beginning, t_end, times, reg_match, nodes);
  dmX = append_col(optima_matrix,X); //Given measurement error in X variable, uncomment this nested statement
  //dmX = append_col(optima_matrix,X_obs); //Given no measurement error in X variable, uncomment this nested statement
  V = calc_direct_V(a, sigma2_y,ta, tij);
  L_v = cholesky_decompose(V);
  mu = dmX*optima_beta;
  //Y ~ multi_normal_cholesky(mu , L_v);//Given measurement error in Y variable, uncomment this statement
  //Y_obs ~ normal(Y,Y_error);//Given measurement error in Y variable, uncomment this statement
  //Y_obs ~ multi_normal_cholesky(mu , L_v); //Given no measurement error in Y variable, uncomment this statement
  target += multi_normal_cholesky_lpdf(Y | mu , L_v);
  target += normal_lpdf(Y_obs | Y, Y_error);
}
generated quantities {
  matrix[N,N] V;
  matrix[N,N] inv_V;
  vector[N] mu;
  real g_i;
  real sigma_ii;
  real sigma_i;
  real u_i;
  vector[N] log_lik;
  matrix[N,n_reg+Z_direct] dmX;
  matrix[N,n_reg] optima_matrix;
  vector[n_reg+Z_direct] optima_beta = append_row(optima,beta);

  real sigma2_y = vy*(2*(log(2)/hl));
  real a = log(2)/hl;
  vector[N] rho = (1 - (1 - exp(-a * T_term))./(a * T_term)); 
  //LOO-CV for multivariate normal models
  optima_matrix = calc_optima_matrix(N, n_reg, a, t_beginning, t_end, times, reg_match, nodes);
  dmX = append_col(optima_matrix,X); //Given measurement error in X variable, uncomment this nested statement
  //dmX = append_col(optima_matrix,X_obs); //Given no measurement error in X variable, uncomment this nested statement
  V = calc_direct_V(a, sigma2_y,ta, tij);
  inv_V = inverse(V);
  mu = dmX*optima_beta;
  for(i in 1:N){
      g_i = (inv_V*(Y-mu))[i];
      sigma_ii = inv_V[i,i];
      u_i = Y[i]-g_i/sigma_ii;
      sigma_i = 1/sigma_ii;
      
      log_lik[i] = -0.5*log(2*pi()*sigma_i)-0.5*(square(Y[i]-u_i)/sigma_i);
      }
}

//Blouch OU model reprogrammed - 2023
//Combination of regime model with multitrait direct effect model with mesurement error
//Regime model - for single regime painting and SIMMAPS
//Multilevel model - varying effects - varying intercepts and varying slopes across regimes with covariance between parameters

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
  int reg_tips[N];
}

parameters {
}
transformed parameters{
}
model {

}
generated quantities {
  vector[N] Y_sim;
  vector[N] Y_sim_obs;
  matrix[N,Z_direct] X;
  real optima_bar;
  vector[Z_direct] beta_bar;
  corr_matrix[1+Z_direct] Rho;
  
  matrix[n_reg,1+Z_direct] v;
  vector<lower=0>[1+Z_direct] sigma;
  
  vector[n_reg] optima;
  matrix[n_reg,Z_direct] beta;

  matrix[N,N] V;
  vector[N] mu;
  matrix[N,N] L_v;
  matrix[N,Z_direct] pred_X;
  matrix[N,n_reg+Z_direct] dmX;
  matrix[N,n_reg] optima_matrix;
  vector[1+Z_direct] ab_bar;

  real<lower=0> hl = lognormal_rng(log(0.25),0.75);
  real<lower=0> vy = exponential_rng(20);
  real a = log(2)/hl;
  real sigma2_y = vy*(2*(log(2)/hl));
  
  Rho = lkj_corr_rng(2,4);
  optima_bar = normal_rng(-1.179507,0.75);
  for(i in 1:Z_direct){
    beta_bar[i] = normal_rng(6.304451,1.75);  
  }
  
  ab_bar[1] = optima_bar;
  ab_bar[2:(Z_direct+1)] = beta_bar;
  
  
  //sigma ~ exponential(5);
  for (i in 1:(1+Z_direct)){
    sigma[i] = normal_rng(0,1);  
  }
  
for(i in 1:(Z_direct)){//Given measurement error in X variable, uncomment this nested statement
    for(j in 1:N){
      X[j,i] = normal_rng(X_obs[j,i], X_error[j,i]);  
    }
  }
  optima_matrix = calc_optima_matrix(N, n_reg, a, t_beginning, t_end, times, reg_match, nodes); //X data
  
  V = calc_direct_V(a, sigma2_y,ta, tij);
  L_v = cholesky_decompose(V);
  
  for ( i in 1:n_reg ){
    v[i,:] = multi_normal_rng(ab_bar,quad_form_diag(Rho , sigma))';
  }
  
  beta = v[, 2:(Z_direct+1)];
  optima = v[, 1];

  for(i in 1:N){
    mu[i] = optima_matrix[i,]*optima+X[i,]*beta[reg_tips[i],]';//Given measurement error in X variable, uncomment this statement
    }
  Y_sim = multi_normal_cholesky_rng(mu , L_v);//Given measurement error in Y variable, uncomment this statement
  
  for(i in 1:N){
    Y_sim_obs[i] = normal_rng(Y_sim[i],Y_error[i]); //Given measurement error in Y variable, uncomment this statement
  }
}

//Blouch OU model reprogrammed
//Using Hansen (1997)
//Regime model - for single regime painting and SIMMAPS
//Multilevel model - varying intercepts across regimes - non-centered transformation

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

  matrix calc_V(real a,real sigma2_y,matrix ta, matrix tij) {
        int N = dims(ta)[1];
        matrix[N, N] Vt;
        Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)
        return(Vt);
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
}

data {
  int N; //Numnber of tips
  int n_reg; //Number of regimes
  int max_node_num; //Max number of nodes in lineage
  vector[N] Y_obs; //Y observed
  vector[N] Y_error;
  matrix[N,N] ta; //Time from tip to ancestor
  matrix[N,N] tij;
  matrix[N, max_node_num] t_beginning; //Matrix of times for beginning of segments to node
  matrix[N, max_node_num] t_end; //Matrix of times for end of segments to
  matrix[N, max_node_num] times; //Matrix of root to node times
  matrix[N, max_node_num] reg_match; //Matrix of 1,2,3 denoting each regime for each node in a lineage. 0 if no node
  int nodes[N]; //Vector of number of nodes per lineage
  vector[2] hl_prior;
  real vy_prior;
  vector[2] optima_prior;
  vector[2] sigma_prior;

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
  matrix[N,n_reg] dmX;
  vector[N] mu;
  vector[N] Y_sim_obs;
  vector[N] Y_sim;
  real<lower=0> hl = lognormal_rng(hl_prior[1],hl_prior[2]);
  real<lower=0> vy = exponential_rng(vy_prior);
  real<lower=0> sigma = abs(normal_rng(sigma_prior[1],sigma_prior[2]));
  vector[n_reg] z;
  real optima_bar = normal_rng(optima_prior[1],optima_prior[2]);
  for(i in 1:n_reg){
    z[i] = normal_rng(0,0.5);
  }
  vector[n_reg] optima = optima_bar + z*sigma;
  real sigma2_y = vy*(2*(log(2)/hl));
  real a = log(2)/hl;
  V = calc_V(a, sigma2_y,ta, tij);
  L_v = cholesky_decompose(V);
  dmX = calc_optima_matrix(N, n_reg, a, t_beginning, t_end, times, reg_match, nodes);
  mu = dmX*optima;
  Y_sim = multi_normal_cholesky_rng(mu , L_v);
  for(i in 1:N){
    Y_sim_obs[i] = normal_rng(Y_sim[i],Y_error[i]);
  }

}


//Blouch OU model reprogrammed
//Using Hansen (1997) 
//Regime model - for multi-regime painting and SIMMAPS
//Multilevel model partial pooling across trees
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
    for (i in 1:dims(x)[1]) {
    //for (i in 1:size(x)) {
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
    vector[n_reg] reg_weights = rep_vector(0,n_reg);
    //print(weight_seg);
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
  int T; //number of trees
  int max_node_num; //Max number of nodes in lineage
  vector[N] Y_obs; //Y observed
  matrix[N,N] ta; //Time from tip to ancestor
  matrix[N,N] tij;
  matrix[N, max_node_num] t_beginning[T]; //Matrix of times for beginning of segments
  matrix[N, max_node_num] t_end[T]; //Matrix of times for end of segments
  matrix[N, max_node_num] times[T]; //Matrix of root to node times
  matrix[N, max_node_num] reg_match[T]; //Matrix of 1,2,3 denoting each regime for each node in a lineage. 0 if no node
  int nodes[N,T]; //Matrix of number of nodes/segments per lineage
}

parameters {
  matrix[T,2] v; //variance covariance matrix between hl and sigms2y
  //matrix[T,n_reg] o; //V/CV matrix for pooling across regimes
  corr_matrix[2] Rho;
  vector<lower=0>[2] hbar_sbar;
  vector<lower=0>[2] sigma;
  matrix[T,n_reg] optima; //Regime Coefficients
  vector[n_reg] optima_bar; //
  vector<lower=0>[n_reg] optima_sigma; //R
}

transformed parameters{
    vector[T] hl;
    vector[T] sigma2_y;
    hl = v[, 1];
    sigma2_y = v[, 2];
}

model {
  matrix[N,N] V;
  vector[N] mu;
  matrix[N,N] L_v;
  matrix[N,n_reg] dmX;
  real a;
  hbar_sbar ~ normal(0,1);
  optima_bar ~ normal(0,1);
  optima_sigma ~ exponential(1);
  Rho ~ lkj_corr( 4 );
  sigma ~ exponential( 1 );
  for(i in 1:T){
    v[i,:] ~ multi_normal(hbar_sbar , quad_form_diag(Rho , sigma)); //Two features - hl and sigma2y using partial pooling across features and across trees
    a = log(2)/hl[i];
    V = calc_V(a, sigma2_y[i],ta, tij);
    L_v = cholesky_decompose(V);
    dmX = calc_optima_matrix(N, n_reg, a, t_beginning[i], t_end[i], times[i], reg_match[i], nodes[,i]);
    optima[i,]~normal(optima_bar,optima_sigma);
    mu = dmX*optima[i,]';
  }
  Y_obs ~ multi_normal_cholesky(mu , L_v);
}
generated quantities {
  real vy = hlbar_sigma2ybar[2]/(2*(log(2)/hlbar_sigma2ybar[1]));
}

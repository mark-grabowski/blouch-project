num_matches<-function(x, y) { #Thanks to Stan Admin Jonah -https:#discourse.mc-stan.org/t/how-to-find-the-location-of-a-value-in-a-vector/19768/2
    n = 0;
    for (i in 1:rows(x))
      if (x[i] == y)
        n = n+1;
    return(n);
  }
  
which_equal<-function(x, y) {
    match_positions[num_matches(x, y)];
    post<-1
    #for (i in 1:size(x)) {
      for (i in 1:(dims(x)[1])) {  
        if (x[i] == y) {
          match_positions[pos] = i;
          pos<-pos+1;
        }
      }
      return(match_positions);
    }
    
weight_segments<-function( a,  t_beginning,  t_end,  time,  nodes){#Individual lineage, calculate weights per segment
      weights = append_row(exp(-a * t_beginning) - exp(-a * t_end),exp(-a * time));
      return(weights);
    }
    
weights_regimes<-function( n_reg,  a,  t_beginning,  t_end,  time,  reg_match,  nodes){#
        #Individual lineage, calculate weights for regimes on each segement
      weight_seg = weight_segments(a, t_beginning[1:(nodes-1)], t_end[1:(nodes-1)], time, nodes);
      #print(weight_seg);
      reg_weights = rep_vector(0,n_reg);
      for(i in 1:n_reg){#reg_match should have values 1,2,3 denoting different regimes
        ids[num_matches(reg_match, i)] = which_equal(reg_match, i); #Returns indixes of matching regimes in weight_segments vector
        #print(ids);
        #print(weight_seg[ids]);
        reg_weights[i] = sum(weight_seg[ids]);
        #print(reg_weights[i]);
      }
      return(t(reg_weights));
  }
  
calc_optima_matrix<-function( N,  n_reg,  a,  t_beginning,  t_end,  times,  reg_match, nodes){
    optima_matrix = rep_matrix(0,N,n_reg);
    for(i in 1:N){ #For each tip/lineage, figure out weighting of regimes
      optima_matrix[i,] = weights_regimes(n_reg, a, t(t_beginning[i,]), t(t_end[i,]), times[i,1], t(reg_match[i,]), nodes[i]);
    }
    return(optima_matrix);
  }
calc_direct_V<-function(  a, sigma2_y, ta,  tij) {
    Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) %*% exp(-a * tij)); #From Hansen (1997)
    return(Vt);
  }

post.ve.link.fxn<-function(dat){
  #post.ve.link.fxn<-function(N,n_reg,Z_direct,Z_X_error,max_node_num,Y_obs,X_obs,Y_error,X_error,ta,tja,
                             #T_term,t_beginning,t_end,times,reg_match,nodes,reg_tips,hl,vy,optima,beta){
    #Using antler post.ve
  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));
  optima_matrix<-calc_optima_matrix(N, n_reg, a, t_beginning, t_end, times, reg_match, nodes); #X data
  V = calc_direct_V(a, sigma2_y,ta, tij);
  L_v = chol(V);
  mu<-rep(0,N)
  for(i in 1:N){
    mu[i]<-optima_matrix[i,]*optima[1,]+X[i,]*t(beta[reg_tips[i],1])
  }
return(mu)
}

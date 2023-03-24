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

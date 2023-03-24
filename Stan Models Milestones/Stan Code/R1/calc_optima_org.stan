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
      //for(j in 1:max_node_length){ //Commented out in revised version - what is this doing?
        optima_matrix[i,z] = sum(regimes_matrix_sep[i,,z]);
      //}
    }
  }

return(optima_matrix);
}
  
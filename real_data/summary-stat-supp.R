library(igraph)

# Number of hyperedges =========================================================
no_hypedge <- function(hyp_set) {
  return(as.numeric(length(hyp_set)))
  }

# Number of vertices ===========================================================

no_vertex <- function(hyp_set){
  return(as.numeric(length(unique(unlist(hyp_set)))))
}
# Weighted Adjacency Matrix ====================================================
weighted_adj <- function(hyp_set){
  
  m <- as.numeric(length(hyp_set))
  n <- as.numeric(length(unique(unlist(hyp_set))))
  n_max <- max(unique(unlist(hyp_set)))
  
  wt.A <- matrix(0,n_max,n_max)
  for (ii in 1:m) {
    wt.A[hyp_set[[ii]],hyp_set[[ii]]] <- wt.A[hyp_set[[ii]],hyp_set[[ii]]] + 1
  }
  diag(wt.A) <- 0
  
  return(wt.A[sort(unique(unlist(hyp_set))), sort(unique(unlist(hyp_set)))])
}


# Binary Adjacency Matrix ======================================================
binary_adj <- function(hyp_set){
  return(ifelse(weighted_adj(hyp_set) >= 1, 1, 0))
}
# Degree distribution ========================================================== 
total_collaboration <- function(hyp_set){
  wt.A <-  weighted_adj(hyp_set)
  return(rowSums(wt.A))
}

unique_collaborator <- function(hyp_set){
  return(rowSums(binary_adj(hyp_set)))
}

unique_collaboration <- function(hyp_set){
  n <- no_vertex(hyp_set)
  m <- no_hypedge(hyp_set)
  number_hyp_incl <- unlist(lapply(1:n, function(i){
    sum(unlist(lapply(1:m, function(j) i %in% hyp_set[[j]] )))
  }))
  return(number_hyp_incl)
}



# Type 2 two stars =============================================================

type2_twostars <- function(hyp_set, apx_itr){
  # ============================================================================
  get_common_vertices <- function(x, y){
    return(sort(union(x,y)))
  }
  # ============================================================================
  restricted_adjacency_matrix <- function(hyp_set, sub_n){
    wt_A <- matrix(0, nrow = sub_n, ncol = sub_n)
    wt_A[hyp_set, hyp_set] <- 1
    diag(wt_A) <- 0
    return(wt_A)
  }
  # ============================================================================
  twostar_2_func <- function(mat1, mat2){
    denom <- sum(mat1 %*% mat2 - mat1 * mat2, na.rm = T)
    return(c( denom))
  }
  # ============================================================================
  color_twostar_count_2 <- function(hyp_set, m,n, apx_itr){
    twostar_list <- lapply(1:apx_itr, function(apx_itr_i){
      subset_m2 <- sample(1:m, 2, replace = F)
      
      all_vertices <- get_common_vertices(hyp_set[[subset_m2[1]]], hyp_set[[subset_m2[2]]])
      
      sub_n <- length(all_vertices)
      
      wt_A_i  <- restricted_adjacency_matrix(match(hyp_set[[subset_m2[1]]], all_vertices), sub_n)
      wt_A_j  <- restricted_adjacency_matrix(match(hyp_set[[subset_m2[2]]], all_vertices), sub_n)
      
      return(twostar_2_func(wt_A_i, wt_A_j))
    })
    return(sum(unlist(twostar_list)))
  }
  m <- as.numeric(length(hyp_set))
  n <- as.numeric(length(unique(unlist(hyp_set))))
  return(color_twostar_count_2(hyp_set, m,n, apx_itr)/apx_itr)
}

# Type 2 clustering coefficient ================================================
type2_clustering_coefficient <- function(hyp_set, apx_itr){
  
  # ==============================================================================
  get_common_vertices <- function(x, y){
    return(sort(union(x,y)))
  }
  # ==============================================================================
  restricted_adjacency_matrix <- function(hyp_edge, sub_n){
    wt_A <- matrix(0, nrow = sub_n, ncol = sub_n)
    wt_A[hyp_edge, hyp_edge] <- 1
    diag(wt_A) <- 0
    return(wt_A)
  }
  # ==============================================================================
  clus_2_mat_tr <- function(mat1, mat2){
    num <- sum(diag(mat1 %*% mat2 %*% mat2 + mat1 %*% mat1 %*% mat2), na.rm = T) 
    denom <- sum(mat1 %*% mat2 - mat1 * mat2, na.rm = T)
    return(c(num, denom))
  }
  # ==============================================================================
  color_clus_coeff_2 <- function(hyp_set, m,n, apx_itr){
    clus_list <- lapply(1:apx_itr, function(apx_itr_i){
      subset_m2 <- sample(1:m, 2, replace = F)
      
      all_vertices <- get_common_vertices(hyp_set[[subset_m2[1]]], hyp_set[[subset_m2[2]]])
      
      sub_n <- length(all_vertices)
      
      wt_A_i  <- restricted_adjacency_matrix(match(hyp_set[[subset_m2[1]]], all_vertices), sub_n)
      wt_A_j  <- restricted_adjacency_matrix(match(hyp_set[[subset_m2[2]]], all_vertices), sub_n)
      
      return(clus_2_mat_tr(wt_A_i, wt_A_j))
    })
    return( sum(unlist(clus_list)[-2*c(1:apx_itr)])/sum(unlist(clus_list)[2*c(1:apx_itr)]) )
  }
  m <- as.numeric(length(hyp_set))
  n <- as.numeric(length(unique(unlist(hyp_set))))
  return(color_clus_coeff_2(hyp_set, m,n, apx_itr))
}


# Binarized edge density =======================================================
bin_edge_density <- function(hyp_set){
  return(edge_density(graph_from_adjacency_matrix(binary_adj(hyp_set))))
}

# Binarized two stars ==========================================================
bin_twostars_scaled <-  function(hyp_set){
  Adj_mat <-  binary_adj(hyp_set)
  n_v <- nrow(Adj_mat)
  rho <- sum(Adj_mat)/(n_v*(n_v-1))
  deg <- sum(choose(rowSums(Adj_mat), 2))
    
  return((deg/choose(n_v,3))/(rho^2))
}

# Binarized cluserting coeff ===================================================
bin_clustering_coefficient_scaled <- function(hyp_set){
  Adj_mat <-  binary_adj(hyp_set)
  n_v <- nrow(Adj_mat)
  graph_gen <- graph_from_adjacency_matrix(Adj_mat)
  rho <- sum(Adj_mat)/(n_v*(n_v-1))
    # transitivity/sparsity
  return(transitivity(graph_gen)/rho) 
}
# ==============================================================================
# Function to calculate A union B

get_common_vertices <- function(x, y){
  return(sort(union(x,y)))
}

# ==============================================================================
# Given a hyperedge (as a vector) produces an adjacency matrix of a subset of vertices

restricted_adjacency_matrix <- function(hyp_edge, sub_n){
  wt_A <- matrix(0, nrow = sub_n, ncol = sub_n)
  wt_A[hyp_edge, hyp_edge] <- 1
  diag(wt_A) <- 0
  return(wt_A)
}

# ==============================================================================
# Function to compute Type 2 clustering coefficient given two adjacency matrices 
# or a subset
# num = Type 2 triangles
# denom = Type 2 two-stars

clus_2_mat_tr <- function(mat1, mat2){
  num <- sum(diag(mat1 %*% mat2 %*% mat2 + mat1 %*% mat1 %*% mat2), na.rm = T) 
  denom <- sum(mat1 %*% mat2 - mat1 * mat2, na.rm = T)
  return(c(num, denom))
}

# ==============================================================================
# Function to compute the Type 2 clustering coefficient of a hypergraph 
# hyp_set = hypergraph as a list
# m = hyper of hyperedges
# n = number of vertices
# apx_itr = approximate number of iteration for incomplete U-stat 

color_clus_coeff_2 <- function(hyp_set, m,n, apx_itr){
  clus_list <- lapply(1:apx_itr, function(apx_itr_i){
    # sample a pair of hyperedges randomly
    subset_m2 <- sample(1:m, 2, replace = F)
    
    all_vertices <- get_common_vertices(hyp_set[[subset_m2[1]]], hyp_set[[subset_m2[2]]])
    
    sub_n <- length(all_vertices)
    
    wt_A_i  <- restricted_adjacency_matrix(match(hyp_set[[subset_m2[1]]], all_vertices), sub_n)
    wt_A_j  <- restricted_adjacency_matrix(match(hyp_set[[subset_m2[2]]], all_vertices), sub_n)
    # Type 2 triangle and Type 2 two-stars
    return(clus_2_mat_tr(wt_A_i, wt_A_j))
  })
  # compute the Type 2 clustering coefficient for the hypergraph
  return( sum(unlist(clus_list)[-2*c(1:apx_itr)])/sum(unlist(clus_list)[2*c(1:apx_itr)]) )
}

# ==============================================================================
# Main function for subsampling
# n = number of vertices
# m = hyper of hyperedges
# hyp_set = hypergraph as a list
# sub.rep = Number of subsamples
# s.m = subsample size
# apx_itr_sub = approximate number of iteration for incomplete U-stat for subsamples 

get.val2 <- function(n, m, hyp_set, sub.rep, s.m, apx_itr_sub){
  # Subsamplimg iteration for each hyper graph
  sub.clus.ct <- mclapply(1:sub.rep, function(sub_itr){
    # choose sample of size s.m
    samp.hyp <- sample(1:m, s.m, replace = F)
    sub.hyp_set <- sapply(samp.hyp, function(list_i) hyp_set[[list_i]])
    sub.clus <-  color_clus_coeff_2(sub.hyp_set, s.m, n, apx_itr_sub)
    return(sub.clus)
  }, mc.cores = n.core)
  return(unlist(sub.clus.ct))
}

# ==============================================================================
# Function to compute CI for two sample testing
# _num = results from sample 1
# _denom = results from sample 2
# m = number of hyperedges
# s.m = subsample size
# res = subsampling result from "get.val2"
# centering = statistic of on the full data
# This function returns both CI for num/denom and denom/num

type2_clus_coeff_sub_ci <- function(m_num, s.m_num, res_num, centering_num, 
                                    m_denom, s.m_denom, res_denom, centering_denom, 
                                    alpha = 0.95){
  z <- qnorm(1-alpha/2)
  which.finite <- function(x) return(x[is.finite(x)])
  
  sub_val <- sqrt(s.m_num)*unlist(res_num)/unlist(res_denom)
  sub_val <- which.finite(sub_val)
  
  CI1 <- c(centering_num/centering_denom - 
             z *sqrt(mean((sub_val-centering_num/centering_denom )^2))/sqrt(m_num),
           centering_num/centering_denom + 
             z *sqrt(mean((sub_val-centering_num/centering_denom )^2))/sqrt(m_num))
  
  sub_val <- sqrt(s.m_denom)*unlist(res_denom)/unlist(res_num)
  sub_val <- which.finite(sub_val)
  
  CI2 <- c(centering_denom/centering_num - 
             z *sqrt(mean((sub_val-centering_denom/centering_num )^2))/sqrt(m_denom),
           centering_denom/centering_num + 
             z *sqrt(mean((sub_val-centering_denom/centering_num )^2))/sqrt(m_denom))
  
  return(list(CI1, CI2))
}

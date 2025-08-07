library(igraph)
library(rje)
library(parallel)

# Data generation ==============================================================
# n = number of vertices
# m = number of hyperedges
# n.prob = node appearance probability

data.gen <- function(n,m,n.prob) {
  hyp.set <- list()
  
  # two step process:
  # step 1: select the number of vertices
  # step 2: select the set of vertices according 1/j^2 given the hyperedge size
  
  for (ii in 1:m){
    # select the hyperedge size
    n.sam <- sample(2:n, 1, prob = n.prob) 
    # select the hyperedge given the size
    n.nodes <- sort(sample(1:n, n.sam, replace = F, 
                           prob = (1/c(1:n)^2)/sum(1/c(1:n)^2))) 
    
    hyp.set[[ii]] <- n.nodes # store hyperedge as a vector
  }
  return(hyp.set)
}

# ==============================================================================
# Function to calculate A union B union C

get_common_vertices <- function(x, y, z){
  return(sort(union(union(x,y), z)))
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
# Function to compute Type 3 clustering coefficient given two adjacency matrices 
# or a subset
# num = Type 3 triangles
# denom = Type 2 two-stars

clus_3_func <- function(mat1, mat2, mat3){
  num <- sum(diag(mat1 %*% mat2 %*% mat3), na.rm = T) 
  denom <- sum(mat1 %*% mat2 - mat1 * mat2, na.rm = T) +
    sum(mat2 %*% mat3 - mat2 * mat3, na.rm = T) +
    sum(mat1 %*% mat3 - mat1 * mat3, na.rm = T)
  return(c(num, denom))
}

# ==============================================================================
# Function to compute the Type 3 clustering coefficient of a hypergraph 
# hyp_set = hypergraph as a list
# m = hyper of hyperedges
# n = number of vertices
# apx_itr = approximate number of iteration for incomplete U-stat 
# filter_id = degree filtered vertices

color_clus_coeff_3 <- function(hyp_set, m,n, apx_itr, filter_id){
  clus_list <- lapply(1:apx_itr, function(apx_itr_i){
    # sample a triplet of hyperedges randomly
    subset_m3 <- sort(sample(1:m, 3, replace = F))
    
    all_vertices <- get_common_vertices(hyp_set[[subset_m3[1]]], 
                                        hyp_set[[subset_m3[2]]],
                                        hyp_set[[subset_m3[3]]])
    
    sub_n <- length(all_vertices)
    
    wt_A_i  <- restricted_adjacency_matrix(
      match(intersect(hyp_set[[subset_m3[1]]], filter_id),all_vertices),
      sub_n)
    wt_A_j  <- restricted_adjacency_matrix(
      match(intersect(hyp_set[[subset_m3[2]]], filter_id),all_vertices),
      sub_n)
    wt_A_k  <- restricted_adjacency_matrix(
      match(intersect(hyp_set[[subset_m3[3]]], filter_id),all_vertices),
      sub_n)
  
    # Type 3 triangle and Type 2 two-stars
    return(clus_3_func(wt_A_i, wt_A_j, wt_A_k))
  })
  return( sum(unlist(clus_list)[-2*c(1:apx_itr)])/sum(unlist(clus_list)[2*c(1:apx_itr)]) )
}

# Wrapper function =============================================================
# Function generates the data and Type 3 clustering coefficient
# n = number of vertices
# m = number of hyperedges
# n.prob = node appearance probability
# apx_itr = approximate number of iteration for incomplete U-stat 
# d = degree filtering, default : d = 1 (no filtering)

get.val1 <-  function(n,m, n.prob, apx_itr, d = 1){
  hyp_set <- data.gen(n,m,n.prob)
  
  # Form weighted adjacency matrix
  wt.A <- matrix(0, ncol = n, nrow = n)
  
  for (ii in 1:m){
    wt.A[hyp_set[[ii]], hyp_set[[ii]]] <- wt.A[hyp_set[[ii]], hyp_set[[ii]]]+1
  }
  diag(wt.A) <- 0
  filter_id <- c(1:n)[rowSums(wt.A) >= d]
  
  n.cluscoeff <- color_clus_coeff_3(hyp_set, m, n, apx_itr, filter_id)
  return(n.cluscoeff)
}

# ==============================================================================
# Main function for subsampling
# n = number of vertices
# m = hyper of hyperedges
# hyp_set = hypergraph as a list
# sub.rep = Number of subsamples
# s.m = subsample size
# apx_itr = approximate number of iteration for incomplete U-stat 
# apx_itr_sub = approximate number of iteration for incomplete U-stat for subsamples 
# d = degree filtering, default : d = 1 (no filtering)

get.val2 <- function(n, m, n.prob, sub.rep, s.m, apx_itr, apx_itr_sub, d){
  hyp.set <- data.gen(n,m,n.prob) # generate a hyper graph
  
  wt.A <- matrix(0, ncol = n, nrow = n)
  
  for (ii in 1:m){
    wt.A[hyp.set[[ii]], hyp.set[[ii]]] <- wt.A[hyp.set[[ii]], hyp.set[[ii]]]+1
  }
  diag(wt.A) <- 0
  filter_id <- c(1:n)[rowSums(wt.A) >= d]
  
  sample_val <- color_clus_coeff_3(hyp.set,m,n,apx_itr,filter_id)
  
  # Subsamplimg iterations for each hyper graph
  sub.clus.ct <- lapply(1:sub.rep, function(sub_itr){
    # choose sample of size s.m
    samp.hyp <- sort(sample(1:m, s.m, replace = F))
    sub.hyp.set <- sapply(samp.hyp, function(list_i) hyp.set[[list_i]])
    sub.clus <-  color_clus_coeff_3(sub.hyp.set, s.m, n, apx_itr_sub, filter_id)
    return(sub.clus)
  })
  return(c(sample_val,unlist(sub.clus.ct)))
}

# Empirical coverage for CI using subsampling ==================================
# type3.clus.true = result from get.val1
# type3.clus.sub = result from get.val2
# m = number of hyperedges
# s.m = s.m subsample sizes
# sub.MC.rep = number of MC iteration for subsampling

type3.clus.coeff.sub.ci <- function(type3.clus.true, type3.clus.sub,
                                    m, s.m, sub.MC.rep, alpha = 0.05){
  return(mean(unlist(lapply(1:sub.MC.rep, function(ii){
    z <- qnorm(1-alpha/2)
    # True value
    trueval <- mean(unlist(type3.clus.true))
    
    sub_res_list <- unlist(type3.clus.sub[[ii]])
    
    # subsample estimates
    sub_res <- sub_res_list[-1]
    
    # Data value
    sub_est <- sub_res_list[1]
    
    return( trueval <= sub_est +  sqrt(s.m * var (sub_res)/ m) * z &
              trueval >= sub_est -  sqrt(s.m * var (sub_res)/m) * z )
  }))))
}
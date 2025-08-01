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
# Function to calculate A union B

get_common_vertices <- function(x, y){
  return(sort(union(x,y)))
}

# ==============================================================================
# Given a hyperedge (as a vector) produces an adjacency matrix of a subset of vertices

restricted_adjacency_matrix <- function(hyp_edge, n, restricted_vertices){
  wt_A <- matrix(0, nrow = n, ncol = n)
  wt_A[hyp_edge, hyp_edge] <- 1
  diag(wt_A) <- 0
  return(wt_A[restricted_vertices, restricted_vertices])
}

# ==============================================================================
# Function to compute Type 2 two-stars given two adjacency matrices or a subset

twostar_2_func <- function(mat1, mat2){
  return(sum( mat1 %*% mat2 - mat1 * mat2))
}

# ==============================================================================
# Function to compute the Type 2 two-stars of a hypergraph 
# hyp_set = hypergraph as a list
# m = hyper of hyperedges
# n = number of vertices
# apx_itr = approximate number of iteration for incomplete U-stat 
# filter_id = degree filtered vertices

color_twostar_count_2 <- function(hyp_set, m, n, apx_itr, filter_id){
  return(sum(unlist(lapply(1:apx_itr, function(apx_itr_i){
    # sample a pair of hyperedges randomly
    subset_m2 <- sample(1:m, 2, replace = F)
    
    all_vertices <- get_common_vertices(hyp_set[[subset_m2[1]]], hyp_set[[subset_m2[2]]])
    
    deg_fil_vertices <- intersect(all_vertices, filter_id)
    
    wt_A_i  <- restricted_adjacency_matrix(hyp_set[[subset_m2[1]]], n, deg_fil_vertices)
    wt_A_j  <- restricted_adjacency_matrix(hyp_set[[subset_m2[2]]], n, deg_fil_vertices)
    # Type 2 two-stars
    return(twostar_2_func(wt_A_i, wt_A_j))
  }))))
}

# Wrapper function =============================================================
# Function generates the data and count colorless triangles
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
  
  n.twostar <- color_twostar_count_2(hyp_set, m, n, apx_itr, 1:n)
  return(n.twostar)
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
  
  sample_val <- color_twostar_count_2(hyp.set,m,n,apx_itr,filter_id)

  # Subsamplimg iterations for each hyper graph
  sub.twostar.ct <- lapply(1:sub.rep, function(sub_itr){
    # choose sample of size s.m
    samp.hyp <- sample(1:m, s.m, replace = F)
    sub.hyp.set <- sapply(samp.hyp, function(list_i) hyp.set[[list_i]])
    sub.twostar <-  color_twostar_count_2(sub.hyp.set, s.m, n, apx_itr_sub, filter_id)
  })
  return(c(sample_val,unlist(sub.twostar.ct)))
}

# Empirical coverage for CI using subsampling ==================================
# type2.twostar.true = result from get.val1
# type2.twostar.sub = result from get.val2
# m = number of hyperedges
# s.m = s.m subsample sizes
# sub.MC.rep = number of MC iteration for subsampling

type2.twostar.sub.ci <- function(type2.twostar.true, type2.twostar.sub,
                                 m, apx_itr_sub, sub.MC.rep, alpha = 0.05){
  return(mean(unlist(lapply(1:sub.MC.rep, function(ii){
    z <- qnorm(1-alpha/2)
    # True value
    trueval <- mean(unlist(type2.twostar.true)/apx_itr)
    
    sub_res_list <- unlist(type2.twostar.sub[[ii]])
    
    # subsample estimates
    sub_res <- sub_res_list[-1]/apx_itr_sub
    
    # Data value
    sub_est <- sub_res_list[1]/apx_itr
    
    return( trueval <= sub_est +  sqrt(apx_itr_sub * var (sub_res)/ m) * z &
              trueval >= sub_est -  sqrt(apx_itr_sub * var (sub_res)/m) * z )
  }))))
}
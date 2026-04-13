library(igraph)
library(rje)
library(pbmcapply)

# Data generation ==============================================================
# n = number of vertices
# m = number of hyperedges
# n.prob = node appearance probability

data.gen <- function(n,m,n.prob, exponent) {
  hyp.set <- list()
  
  # two step process:
  # step 1: select the number of vertices
  # step 2: select the set of vertices according 1/j^2 given the hyperedge size
  
  for (ii in 1:m){
    # select the hyperedge size
    n.sam <- sample(2:n, 1, prob = n.prob) 
    # select the hyperedge given the size
    n.nodes <- sort(sample(1:n, n.sam, replace = F, 
                           prob = (1/c(1:n)^exponent)/sum(1/c(1:n)^exponent))) 
    
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

restricted_adjacency_matrix <- function(hyp_edge, sub_n){
  wt_A <- matrix(0, nrow = sub_n, ncol = sub_n)
  wt_A[hyp_edge, hyp_edge] <- 1
  diag(wt_A) <- 0
  return(wt_A)
}

# ==============================================================================
# Function to compute Type 2 triangles given two adjacency matrices or a subset

triangles_2_func <- function(mat1, mat2){
  return(sum(c(sum(diag(mat1 %*% mat2 %*% mat2))), sum(diag(mat1 %*% mat1 %*% mat2))))
}

# ==============================================================================
# Function to compute the Type 2 triangles of a hypergraph 
# hyp.set = hypergraph as a list
# m = hyper of hyperedges
# n = number of vertices
# apx_itr = approximate number of iteration for incomplete U-stat 
# filter_id = degree filtered vertices

color_triangles_count_2 <- function(hyp.set, m, n, apx_itr, filter_id){
  return(sum(unlist(lapply(1:apx_itr, function(apx_itr_i){
    # sample a pair of hyperedges randomly
    subset_m2 <- sample(1:m, 2, replace = F)
    
    all_vertices <- get_common_vertices(hyp.set[[subset_m2[1]]], 
                                        hyp.set[[subset_m2[2]]])
    
    sub_n <- length(all_vertices)
    
    wt_A_i  <- restricted_adjacency_matrix(
      match(intersect(hyp.set[[subset_m2[1]]], filter_id),all_vertices),
      sub_n)
    wt_A_j  <- restricted_adjacency_matrix(
      match(intersect(hyp.set[[subset_m2[2]]], filter_id),all_vertices),
      sub_n)
    # Type 2 triangles
    return(triangles_2_func(wt_A_i, wt_A_j))
  }))))
}


# Wrapper function =============================================================
# Function generates the data and count Type 2 triangles
# n = number of vertices
# m = number of hyperedges
# n.prob = node appearance probability
# apx_itr = approximate number of iteration for incomplete U-stat 
# d = degree filtering, default : d = 1 (no filtering)

get.val1 <-  function(n,m, n.prob, exponent, apx_itr, d = 1){
  hyp.set <- data.gen(n,m,n.prob,exponent)
  
  # Form weighted incedence matrix
  Inc_mat <- matrix(0, ncol = n, nrow = m)
  
  for (ii in 1:m){
    Inc_mat[ii,hyp.set[[ii]]] <- 1
  }
  filter_id <- c(1:n)[colSums(Inc_mat) >= d]
  
  n.tri <- color_triangles_count_2(hyp.set, m, n, apx_itr, filter_id)
  return(n.tri)
}

# ==============================================================================
# Main function for subsampling
# n = number of vertices
# m = hyper of hyperedges
# hyp.set = hypergraph as a list
# sub.rep = Number of subsamples
# s.m = subsample size
# apx_itr = approximate number of iteration for incomplete U-stat 
# apx_itr_sub = approximate number of iteration for incomplete U-stat for subsamples 
# d = degree filtering, default : d = 1 (no filtering)

get.val2 <- function(n, m, n.prob, exponent, sub.rep, s.m, apx_itr, apx_itr_sub, d){
  hyp.set <- data.gen(n,m,n.prob, exponent) # generate a hyper graph
  
  # Form weighted incedence matrix
  Inc_mat <- matrix(0, ncol = n, nrow = m)
  
  for (ii in 1:m){
    Inc_mat[ii,hyp.set[[ii]]] <- 1
  }
  filter_id <- c(1:n)[colSums(Inc_mat) >= d]
  
  sample_val <- color_triangles_count_2(hyp.set,m,n,apx_itr,filter_id)
  
  # Subsamplimg iterations for each hyper graph
  sub.triangle.ct <- lapply(1:sub.rep, function(sub_itr){
    # choose sample of size s.m
    samp.hyp <- sample(1:m, s.m, replace = F)
    sub.hyp.set <- sapply(samp.hyp, function(list_i) hyp.set[[list_i]])
    sub.triangle <-  color_triangles_count_2(sub.hyp.set, s.m, n, apx_itr_sub, filter_id)
    return(sub.triangle)
  })
  return(c(sample_val,unlist(sub.triangle.ct)))
}

# Empirical coverage for CI using subsampling ==================================
# type2.triangle.true = result from get.val1
# type2.triangle.sub = result from get.val2
# m = number of hyperedges
# s.m = s.m subsample sizes
# apx_itr = approximate number of iteration for incomplete U-stat 
# apx_itr_sub = approximate number of iteration for incomplete U-stat for subsamples 
# sub.MC.rep = number of MC iteration for subsampling

type2.triangle.sub.ci <- function(type2.triangle.true, type2.triangle.sub,
                                  m, apx_itr, apx_itr_sub, sub.MC.rep, 
                                  alpha = 0.05){
  return(mean(unlist(lapply(1:sub.MC.rep, function(ii){
    z <- qnorm(1-alpha/2)
    # True value
    trueval <- mean(unlist(type2.triangle.true)/apx_itr)
    
    sub_res_list <- unlist(type2.triangle.sub[[ii]])
    
    # subsample estimates
    sub_res <- sub_res_list[-1]/apx_itr_sub
    
    # Data value
    sub_est <- sub_res_list[1]/apx_itr
    
    return( trueval <= sub_est +  sqrt(apx_itr_sub * var (sub_res)/ m) * z &
              trueval >= sub_est -  sqrt(apx_itr_sub * var (sub_res)/m) * z )
  }))))
}
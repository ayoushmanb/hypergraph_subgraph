require(RSpectra)
require(igraph)
#===============================================================================
# Vertex subsampling
# b: subsample size (integer)
# N number of subsamples
vertex_subsampling <- function(Adj_mat, b,N, my_function,output_size,...) {
  
  boot_out  <-  matrix(NA,nrow=N, ncol=output_size)
  
  n_row  <-  nrow(Adj_mat)
  
  for(ii in 1: N){
    print(ii)
    boot_indices <-  sample(n_row, size = b, replace=FALSE)
    boot_subgraph  <-  Adj_mat[boot_indices,boot_indices]
    boot_out[ii,]  <-  my_function(boot_subgraph,...) 
  }
  return(boot_out)
}

#===============================================================================
# Subsampling Function
# subsampling_CI  <- function(x,...) UseMethod("subsampling_CI")

subsampling_CI.numeric  <- function(numeric, est.param, centering.val,
                                  tau_seq = function(x) return(x^0.5), 
                                  n, b, alpha=0.05){
  tau_vec  <- Vectorize(tau_seq, vectorize.args="x")  
  
  quantile_values  <- c(1-alpha/2,alpha/2)
  subsample_normalized  <- tau_vec(b)*(numeric-centering.val)
  quantiles_subsample  <- quantile(subsample_normalized,
                                 probs = quantile_values,na.rm=TRUE)
  ci_out  <- est.param - quantiles_subsample/tau_seq(n)
  return(ci_out)
}

# ==============================================================================
twostars_fun_scaled <- function(Adj_mat){
  n_v <- nrow(Adj_mat)
  rho <- sum(Adj_mat)/(n_v*(n_v-1))
  deg <- sum(choose(rowSums(Adj_mat), 2))
  
  return((deg/choose(n_v,3))/(rho^2))
}

# Binary Adjacency Matrix ======================================================
binary_adj <- function(hyp_set){
  m <- as.numeric(length(hyp_set))
  n <- as.numeric(length(unique(unlist(hyp_set))))
  n_max <- max(unique(unlist(hyp_set)))
  
  wt.A <- matrix(0,n_max,n_max)
  for (ii in 1:m) {
    wt.A[hyp_set[[ii]],hyp_set[[ii]]] <- wt.A[hyp_set[[ii]],hyp_set[[ii]]] + 1
  }
  diag(wt.A) <- 0
  
  return(ifelse(wt.A >= 1, 1, 0))
}
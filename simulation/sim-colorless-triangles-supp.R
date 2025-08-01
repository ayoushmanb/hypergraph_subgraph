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

# Colorless triangles ==========================================================
# n = number of vertices
# m = number of hyperedges
# hyp.set = hyperedge list
# d = degree filtering

colorless.tri.func.df <- function(n, m, hyp.set, d) {
  
  # Form weighted adjacency matrix
  wt.A <- matrix(0, ncol = n, nrow = n)
  
  for (ii in 1:m){
    wt.A[hyp.set[[ii]], hyp.set[[ii]]] <- wt.A[hyp.set[[ii]], hyp.set[[ii]]]+1
  }
  diag(wt.A) <- 0
  filter_id <- rowSums(wt.A) >= d
  
  # count triangles 
  tri.ct <- sum(diag(wt.A[filter_id,filter_id] %*% 
                       wt.A[filter_id,filter_id] %*% 
                       wt.A[filter_id,filter_id]))/6
  return(tri.ct)
}

# Wrapper function =============================================================
# Function generates the data and count colorless triangles
# n = number of vertices
# m = number of hyperedges
# n.prob = node appearance probability
# d = degree filtering, default : d = 1 (no filtering)

get.val1 <-  function(n,m, n.prob,d = 1){
  hyp.set <- data.gen(n,m,n.prob)
  n.tri <- colorless.tri.func.df(n,m, hyp.set,d)
  return(n.tri)
}

# Subsampling ==================================================================
# Function performs subsampling
# n = number of vertices
# m = number of hyperedges
# n.prob = node appearance probability
# sub.rep = number of independent subsamples
# s.m = subsample sizes
# d = degree filtering, default : d = 1 (no filtering)

get.val2 <- function(n, m, n.prob, sub.rep, s.m, d = 1){
  hyp.set <- data.gen(n,m,n.prob) # generate a hyper graph
  
  sample_val <- colorless.tri.func.df(n,m,hyp.set,d)
  
  # Subsamplimg iterations for each hyper graph
  sub.tri.ct <- lapply(1:sub.rep, function(sub_itr){
    # choose sample of size s.m
    samp.hyp <- sample(1:m, s.m, replace = F)
    sub.hyp.set <- sapply(samp.hyp, function(list_i) hyp.set[[list_i]])
    # count triangles for subsamples
    sub.tri <-  colorless.tri.func.df(n, s.m, sub.hyp.set,d)
  })
  return(c(sample_val,unlist(sub.tri.ct)))
}

# Empirical coverage for CI using subsampling ==================================
# colorless.tri.true = result from get.val1
# colorless.tri.sub = result from get.val2
# m = number of hyperedges
# s.m = s.m subsample sizes
# sub.MC.rep = number of MC iteration for subsampling


colorless.tri.sub.ci <- function(colorless.tri.true, colorless.tri.sub,
                                 m, s.m, sub.MC.rep, alpha = 0.05){
  return(mean(unlist(lapply(1:sub.MC.rep, function(ii){
    z <- qnorm(1-alpha/2)
    # True value
    trueval <- mean(unlist(colorless.tri.true)/choose(m,3))
    
    sub_res_list <- unlist(colorless.tri.sub[[ii]])
    
    # subsample estimates
    sub_res <- sub_res_list[-1]/choose(s.m,3)
    
    # Data value
    sub_est <- sub_res_list[1]/choose(m,3) 
    
    return( trueval <= sub_est +  sqrt(s.m*mean((sub_res - sub_est)^2)/m) * z &
              trueval >= sub_est -  sqrt(s.m*mean((sub_res - sub_est)^2)/m) * z )
  }))))
}

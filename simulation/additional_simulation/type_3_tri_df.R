library(igraph)
library(rje)
library(parallel)
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
# Function to compute Type 3 triangles given three adjacency matrices or a subset

triangles_3_func <- function(mat1, mat2, mat3){
  return(sum(diag(mat1 %*% mat2 %*% mat3)))
}

# ==============================================================================
# Function to compute the Type 3 triangles of a hypergraph 
# hyp.set = hypergraph as a list
# m = hyper of hyperedges
# n = number of vertices
# apx_itr = approximate number of iteration for incomplete U-stat 
# filter_id = degree filtered vertices


color_triangles_count_3 <- function(hyp.set, m,n, apx_itr, filter_id){
  sum(unlist(lapply(1:apx_itr, function(apx_itr_i){
    # sample a triplet of hyperedges randomly
    subset_m3 <- sample(1:m, 3, replace = F)
    
    all_vertices <- get_common_vertices(hyp.set[[subset_m3[1]]], 
                                        hyp.set[[subset_m3[2]]],
                                        hyp.set[[subset_m3[3]]])
    
    sub_n <- length(all_vertices)
    
    wt_A_i  <- restricted_adjacency_matrix(
      match(intersect(hyp.set[[subset_m3[1]]], filter_id),all_vertices),
      sub_n)
    wt_A_j  <- restricted_adjacency_matrix(
      match(intersect(hyp.set[[subset_m3[2]]], filter_id),all_vertices),
      sub_n)
    wt_A_k  <- restricted_adjacency_matrix(
      match(intersect(hyp.set[[subset_m3[3]]], filter_id),all_vertices),
      sub_n)
    
    # Type 3 triangles
    return(triangles_3_func(wt_A_i, wt_A_j, wt_A_k))
  })))
}

# Wrapper function =============================================================
# Function generates the data and count Type 3 triangles
# n = number of vertices
# m = number of hyperedges
# n.prob = node appearance probability
# apx_itr = approximate number of iteration for incomplete U-stat 
# d = degree filtering, default : d = 1 (no filtering)


get.val1 <-  function(n,m, n.prob,exponent, apx_itr, d_values = 1){
  hyp.set <- data.gen(n,m,n.prob, exponent)
  
  Inc_mat <- matrix(0, ncol = n, nrow = m)
  for (ii in 1:m){
    Inc_mat[ii,hyp.set[[ii]]] <- 1
  }
  
  n_tri_values <- unlist(lapply(d_values, function(d){
    filter_id <- c(1:n)[colSums(Inc_mat) >= d]
    
    n.tri <- color_triangles_count_3(hyp.set, m, n, apx_itr, filter_id)
    return(n.tri)
  }))
}


m <- 500 # number of hyp edges
n <- 1000 # number of vertices
MC.rep <- 200 # MC iteration to approximate true variance
apx_itr <-  round(m^1.5) # approximate iteration for incomplete U-stat
exponent <- 4
# d <- c(1, ceiling(c(0.5, 1, 5, 10, 20)*m^(1/2-1/exponent))) # degree filtering, default : d = 1 (no filtering)
d_values <- round(c(1, m^(c(0.1, 0.167, 0.25, 0.5, 0.6, 0.75, 0.9))))

n.prob <- dpois(2:n, 6) # Hyperedge sizes follow Poi(6)
n.prob <- n.prob/sum(n.prob) # probability


type_3_values <- pbmclapply(1:MC.rep, function(mc_rep){
  return(get.val1(n,m, n.prob, exponent,apx_itr,d_values))
}, mc.cores = 20, mc.set.seed = F, mc.style = "ETA")

type_3_values_mat <- matrix(unlist(type_3_values), byrow = T, nrow = MC.rep)/apx_itr
type_3_values_mat <- sqrt(m)*(type_3_values_mat- mean(type_3_values_mat[,1]))/sd(sqrt(m)*type_3_values_mat[,1])

# =======

library(ggplot2)
library(tidyr)

# plot_density_matrix <- function(mat, c_vals) {
#   
#   df <- as.data.frame(mat)
#   
#   if (is.null(colnames(df))) {
#     colnames(df) <- paste0("V", seq_len(ncol(df)))
#   }
#   
#   df_long <- pivot_longer(
#     df,
#     cols = everything(),
#     names_to = "variable",
#     values_to = "value"
#   )
#   
#   # Build legend expressions: 0, m^c(val)
#   legend_labels <- c(
#     list(expression(0)),
#     lapply(c_vals, function(x) {
#       bquote(m^.(x))
#     })
#   )
#   
#   ggplot(df_long, aes(x = value, color = variable)) +
#     geom_density(linewidth = 1.5) +
#     scale_color_manual(
#       values = scales::hue_pal()(ncol(mat)),
#       labels = legend_labels
#     ) +
#     labs(
#       x = "",
#       y = "",
#       color = expression(d)
#     ) +
#     theme_bw(base_size = 35) +
#     theme(panel.grid = element_blank())
# }

plot_density_matrix <- function(mat, c_vals) {
  
  df <- as.data.frame(mat)
  
  if (is.null(colnames(df))) {
    colnames(df) <- paste0("V", seq_len(ncol(df)))
  }
  
  df_long <- tidyr::pivot_longer(
    df,
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  )
  
  # Legend labels
  legend_labels <- c(
    list(expression(0)),
    lapply(c_vals, function(x) bquote(m^.(x)))
  )
  
  ggplot(df_long, aes(x = value, color = variable)) +
    geom_density(linewidth = 1.2,key_glyph = "point") +
    scale_color_manual(
      values = scales::hue_pal()(ncol(mat)),
      labels = legend_labels
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = 16,     # solid dot
          size = 7        # dot size
        )
      )
    ) +
    labs(
      x = "",
      y = "",
      color = expression(d)
    ) +
    theme_bw(base_size = 28) +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1)
}


pdf(paste0("type_3_tri_df_",m,"_",n,"_exp",exponent,".pdf"))
plot_density_matrix(type_3_values_mat,
                    c(0.1,0.167, 0.25, 0.5, 0.6, 0.75, 0.9)        )
dev.off()

save.image(paste0("type_3_tri_df_",m,"_",n,"_exp",exponent,".RData"))


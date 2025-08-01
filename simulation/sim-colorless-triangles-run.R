setwd("~/hyperedge")
source("sim-colorless-triangles-supp.R")

# Simulation setting ===========================================================

m <- 100 # number of hypperedges
n <- 100 # number of vertices
MC.rep <- 10 # MC iteration to approximate true variance
d <- 1 # degree filtering, default : d = 1 (no filtering)

const <- 0.6 # c(0.6, 0.8, 1, 1.2, 1.4)
s.m <- round(const*m/log(m)) # Subsample sizes
sub.rep <- 10 # Number of subsamples
sub.MC.rep <- 10 # MC iteration for subsampling

n.core <- 5 # number of cores to use

n.prob <- dpois(2:n, 6) # Hyperedge sizes follow Poi(6)
n.prob <- n.prob/sum(n.prob) # probability

# MC approximation =============================================================

ST <- Sys.time()
set.seed(1234)
colorless.tri.true <- unlist(mclapply(1:MC.rep, function(rep){
  get.val1(n,m,n.prob,d)
}, mc.cores =  n.core))
Sys.time()-ST

# Subsampling ==================================================================

Sys.time()
ST <- Sys.time()
colorless.tri.sub <- mclapply(1:sub.MC.rep, function(sub_rep){
  (get.val2(n, m, n.prob, sub.rep, s.m, d))
}, mc.cores = n.core)
Sys.time()-ST

# save.image("file_name.RData") 

# Results ======================================================================

colorless.tri.sub.ci(colorless.tri.true, colorless.tri.sub,
                     m, s.m, sub.MC.rep)


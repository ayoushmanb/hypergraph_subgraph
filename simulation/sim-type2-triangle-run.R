# setwd("~/hyperedge/sim_1000_1000/")
source("sim-type2-triangle-supp.R")

# Simulation setting ===========================================================

m <- 100 # number of hyp edges
n <- 100 # number of vertices
MC.rep <- 10 # MC iteration to approximate true variance
apx_itr <-  round(m^1.1) # approximate iteration for incomplete U-stat
d <- 1 # degree filtering, default : d = 1 (no filtering)

const <- 4 # c(3.8, 4)
s.m <- round(const*m/log(m)) # Subsample sizes
sub.rep <- 10 # 1000 # Number of subsamples
sub.MC.rep <- 10 # 1000 # MC iteration for subsampling
apx_itr_sub <- round(s.m^1.1) # approximate iteration for incomplete U-stat

n.core <- 5 # number of cores to use

n.prob <- dpois(2:n, 6) # Hyperedge sizes follow Poi(6)
n.prob <- n.prob/sum(n.prob) # probability

# MC approximation =============================================================

ST <- Sys.time()
set.seed(1234)
type2.triangle.true <- unlist(mclapply(1:MC.rep, function(rep){
  get.val1(n,m,n.prob, apx_itr, 1)
}, mc.cores =  n.core))
Sys.time()-ST

# Subsampling ==================================================================

Sys.time()
ST <- Sys.time()
type2.triangle.sub <- mclapply(1:sub.MC.rep, function(sub_rep){
  return(get.val2(n, m, n.prob, sub.rep, s.m, apx_itr, apx_itr_sub, d))
}, mc.cores = n.core)
Sys.time()-ST
Sys.time()

# save.image("file_name.RData") 

# Results ======================================================================

type2.triangle.sub.ci(type2.triangle.true, type2.triangle.sub,
                     m, apx_itr_sub, sub.MC.rep, alpha = 0.05)

setwd("~/hyperedge/real-data-test")
source('~/hyperedge/real-data-test/real-data-type2-clustering-test-supp.R')
library(igraph)
library(rje)
library(parallel)

# Subsampling ==================================================================

m1 <- length(hyp_edge_list1) # number of hyp edges
n1 <- length(unique(unlist(hyp_edge_list1))) # number of vertices

m2 <- length(hyp_edge_list2) # number of hyp edges
n2 <- length(unique(unlist(hyp_edge_list2))) # number of vertices

m3 <- length(hyp_edge_list3) # number of hyp edges
n3 <- length(unique(unlist(hyp_edge_list3))) # number of vertices

m4 <- length(hyp_edge_list4) # number of hyp edges
n4 <- length(unique(unlist(hyp_edge_list4))) # number of vertices

const <- 1.5
s.m1 <- round(const*m1/log(m1)) # Subsample sizes
s.m2 <- round(const*m2/log(m2)) # Subsample sizes
s.m3 <- round(const*m3/log(m3)) # Subsample sizes
s.m4 <- round(const*m4/log(m4)) # Subsample sizes
sub.rep <- 1000 # Number of subsamples

n.core <- 10

# Construct the confidence intervals ===========================================

Sys.time()
ST <- Sys.time()
sub.tri.ct.mc1 <- get.val2(n1, m1, hyp_edge_list1, 1000, s.m1, choose(s.m1, 2)) #jcgs
Sys.time()-ST


Sys.time()
ST <- Sys.time()
sub.tri.ct.mc2 <- get.val2(n2, m2, hyp_edge_list2, 1000, s.m2, choose(s.m2, 2)) #jasa
Sys.time()-ST


Sys.time()
ST <- Sys.time()
sub.tri.ct.mc3 <- get.val2(n3, m3, hyp_edge_list3, 1000, s.m3, choose(s.m3, 2)) #emi
Sys.time()-ST


Sys.time()
ST <- Sys.time()
sub.tri.ct.mc4 <- get.val2(n4, m4, hyp_edge_list4, 1000, s.m4, choose(s.m4, 2)) #movie 
Sys.time()-ST

save.image("real-data-test-clus-coeff-2_v2.RData")

# ==============================================================================

# jcgs vs jasa
type2_clus_coeff_sub_ci(m1, s.m1, sub.tri.ct.mc1, 0.2872,
                        m2, s.m2, sub.tri.ct.mc2, 0.1683)
# [1] 1.612212 1.800741
# [1] 0.5068612 0.6651444


# jcgs vs emi
type2_clus_coeff_sub_ci(m1, s.m1, sub.tri.ct.mc1, 0.2872,
                        m3, s.m3, sub.tri.ct.mc3, 0.3588)
# [1] 0.7718633 0.8290286
# [1] 1.131201 1.367406

#  movie vs jcgs
type2_clus_coeff_sub_ci(m4, s.m4, sub.tri.ct.mc4, 0.0430,
                        m1, s.m1, sub.tri.ct.mc1, 0.2872)
# [1] 0.1209605 0.1784824
# [1] 6.242364 7.115776

# jasa vs emi
type2_clus_coeff_sub_ci(m2, s.m2, sub.tri.ct.mc2, 0.1683,
                        m3, s.m3, sub.tri.ct.mc3, 0.3588)
# [1] 0.4501157 0.4880114
# [1] 2.025166 2.238649

# jasa vs movie
type2_clus_coeff_sub_ci(m4, s.m4, sub.tri.ct.mc4, 0.0430,
                        m2, s.m2, sub.tri.ct.mc2, 0.1683)
# [1] 0.2358825 0.2751098
# [1] 3.651743 4.176164


# movie vs emi
type2_clus_coeff_sub_ci(m4, s.m4, sub.tri.ct.mc4, 0.0430,
                        m3, s.m3, sub.tri.ct.mc3, 0.3588)
# [1] 0.1131185 0.1265693
# [1] 7.884328 8.804044

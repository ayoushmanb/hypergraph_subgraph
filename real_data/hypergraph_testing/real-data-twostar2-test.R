setwd("~/hyperedge/real-data-test")
source('~/hyperedge/real-data-test/real-data-twostar2-test-supp.R')
library(igraph)
library(rje)
library(parallel)

# =============================================================================
# Load the data set

# JCGS
m1 <- length(hyp_edge_list1) # number of hyp edges
n1 <- length(unique(unlist(hyp_edge_list1))) # number of vertices

# JASA
m2 <- length(hyp_edge_list2) # number of hyp edges
n2 <- length(unique(unlist(hyp_edge_list2))) # number of vertices

# EMI
m3 <- length(hyp_edge_list3) # number of hyp edges
n3 <- length(unique(unlist(hyp_edge_list3))) # number of vertices

# Movie
m4 <- length(hyp_edge_list4) # number of hyp edges
n4 <- length(unique(unlist(hyp_edge_list4))) # number of vertices

const <- 1.5
s.m1 <- round(const*m1/log(m1)) # Subsample sizes
s.m2 <- round(const*m2/log(m2)) # Subsample sizes
s.m3 <- round(const*m3/log(m3)) # Subsample sizes
s.m4 <- round(const*m4/log(m4)) # Subsample sizes
sub.rep <- 1000 # Number of subsamples

n.core <- 10 # number of cores to use

# =============================================================================

Sys.time()
ST <- Sys.time()
sub.tri.ct.mc1 <- get.val2(n1, m1, hyp_edge_list1, 1000, s.m1, choose(s.m1, 2)) # JCGS
Sys.time()-ST


Sys.time()
ST <- Sys.time()
sub.tri.ct.mc2 <- get.val2(n2, m2, hyp_edge_list2, 1000, s.m2, choose(s.m2, 2)) # JASA
Sys.time()-ST


Sys.time()
ST <- Sys.time()
sub.tri.ct.mc3 <- get.val2(n3, m3, hyp_edge_list3, 1000, s.m3, choose(s.m3, 2)) # EMI
Sys.time()-ST


Sys.time()
ST <- Sys.time()
sub.tri.ct.mc4 <- get.val2(n4, m4, hyp_edge_list4, 1000, s.m4, choose(s.m4, 2)) # Movie
Sys.time()-ST

save.image("real-dat-test-twostar-2_v2.RData")

# ==============================================================================


# jcgs vs jasa
type2_two_star_sub_ci(m1, s.m1, sub.tri.ct.mc1, 0.0157, choose(s.m1,2),
                        m2, s.m2, sub.tri.ct.mc2, 0.0270, choose(s.m2,2))
# [1] 0.5619886 0.6009744
# [1] 1.656108 1.783383


# jcgs vs emi
type2_two_star_sub_ci(m1, s.m1, sub.tri.ct.mc1, 0.2872, choose(s.m1,2),
                        m3, s.m3, sub.tri.ct.mc3, 0.3588, choose(s.m3,2))
# [1] 0.7989713 0.8019206
# [1] -3.313603  5.812210

# movie vs jcgs
type2_two_star_sub_ci(m4, s.m4, sub.tri.ct.mc4, 0.0430, choose(s.m4,2),
                        m1, s.m1, sub.tri.ct.mc1, 0.2872, choose(s.m1,2))
# [1] -0.9830197  1.2824626
# [1] 6.665637 6.692503

# jasa vs emi
type2_two_star_sub_ci(m2, s.m2, sub.tri.ct.mc2, 0.1683, choose(s.m2, 2),
                        m3, s.m3, sub.tri.ct.mc3, 0.3588, choose(s.m3,2))
# [1] 0.4685817 0.4695454
# [1] -0.1321625  4.3959772

# jasa vs wiki
type2_two_star_sub_ci(m4, s.m4, sub.tri.ct.mc4, 0.0430, choose(s.m4,2),
                        m2, s.m2, sub.tri.ct.mc2, 0.1683, choose(s.m2,2))
# [1] -0.3056004  0.8165927
# [1] 3.908154 3.919753


# movie vs emi
type2_two_star_sub_ci(m4, s.m4, sub.tri.ct.mc4, 0.0430, choose(s.m4, 2),
                        m3, s.m3, sub.tri.ct.mc3, 0.3588, choose(s.m3,2))
# [1] 0.1123036 0.1273843
# [1] 8.241471 8.446902


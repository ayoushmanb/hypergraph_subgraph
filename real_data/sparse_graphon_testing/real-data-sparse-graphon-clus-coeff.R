source("~/R/hyp-real-data/real-data-final/real-data-sparse-graphon-clus-coeff-supp.R")

#===============================================================================
# JASA
# load the data set
bin.A.jasa <- binary_adj(hyp_edge_list2)
n.jasa <- nrow(bin.A.jasa)
b.jasa <- ceiling(n.jasa/log(n.jasa))
rho.jasa <- edge_density(graph_from_adjacency_matrix(bin.A.jasa)) # 0.001291248
jasa.true <- transitivity(graph_from_adjacency_matrix(bin.A.jasa)) # 0.4124067

# subsampling
jasa.cc <- vertex_subsampling(bin.A.jasa, b.jasa, 1000,
                              transitivity_fun_scaled, 1)

# # CI using subsampling
# sub.cc.jasa <- subsampling_CI.numeric(jasa.cc, 
#                                       jasa.true/rho.jasa, jasa.true/rho.jasa, 
#                                       n = n.jasa,
#                                       b = b.jasa)

#===============================================================================
# JCGS
# load the data set
bin.A.jcgs <- binary_adj(hyp_edge_list2)
n.jcgs <- nrow(bin.A.jcgs)
b.jcgs <- ceiling(n.jcgs/log(n.jcgs))
rho.jcgs <- edge_density(graph_from_adjacency_matrix(bin.A.jcgs)) # 0.001350759
jcgs.true <- transitivity(graph_from_adjacency_matrix(bin.A.jcgs)) # 0.685331

# subsampling
jcgs.cc <- vertex_subsampling(bin.A.jcgs, b.jcgs, 1000, transitivity_fun_scaled, 1)

# CI using subsampling
# sub.cc.jcgs <- subsampling_CI.numeric(jcgs.cc,
#                                       jcgs.true/rho.jcgs, jcgs.true/rho.jcgs,
#                                       n = n.jcgs,
#                                       b = b.jcgs)

#===============================================================================
# EMI
# load the data set
bin.A.emi <- binary_adj(hyp_edge_list2)
n.emi <- nrow(bin.A.emi)
b.emi <- ceiling(n.emi/log(n.emi))
rho.emi <- edge_density(graph_from_adjacency_matrix(bin.A.emi)) # 0.001014454
emi.true <- transitivity(graph_from_adjacency_matrix(bin.A.emi)) # 0.6199805

# subsampling
emi.cc <- vertex_subsampling(bin.A.emi, b.emi, 1000, transitivity_fun_scaled, 1)

# CI using subsampling
# sub.cc.emi <- subsampling_CI.numeric(emi.cc,
#                                      emi.true/rho.emi, emi.true/rho.emi, 
#                                      n = n.emi , 
#                                      b = b.emi)


#===============================================================================
# Movie
# load the data set
bin.A.mov <- binary_adj(hyp_edge_list2)
n.mov <- nrow(bin.A.mov)
b.mov <- ceiling(n.mov/log(n.mov))
rho.mov <- edge_density(graph_from_adjacency_matrix(bin.A.mov)) # 0.001812831
mov.true <- transitivity(graph_from_adjacency_matrix(bin.A.mov)) # 0.2718159

# subsampling
mov.cc <- vertex_subsampling(bin.A.mov, b.mov, 1000, transitivity_fun_scaled, 1)

# # CI using subsampling
# sub.cc.mov <- subsampling_CI.numeric(mov.cc,
#                                      mov.true/rho.mov, mov.true/rho.mov, 
#                                      n = n.mov ,
#                                      b = b.mov)




#===============================================================================
# JCGS vs JASA
subsampling_CI.numeric(jcgs.cc/(sqrt(b.jasa/b.jcgs)*jasa.cc),
                       (jcgs.true/rho.jcgs)/(jasa.true/rho.jasa), 
                       (jcgs.true/rho.jcgs)/(sqrt(n.jasa/n.jcgs)*jasa.true/rho.jasa),
                       n = n.jcgs *(n.jcgs/n.jasa),
                       b = b.jcgs) 
# 0.1093683 1.9843195

#===============================================================================
# JCGS vs EMI
subsampling_CI.numeric(jcgs.cc/(sqrt(b.emi/b.jcgs)*emi.cc),
                       (jcgs.true/rho.jcgs)/(emi.true/rho.emi), 
                       (jcgs.true/rho.jcgs)/(sqrt(n.emi/n.jcgs)*emi.true/rho.emi),
                       n = n.jcgs*(n.jcgs/n.emi),
                       b = b.jcgs) 
# 0.5771839 0.9805596

#===============================================================================
# mov vs jcgs
subsampling_CI.numeric(mov.cc/(sqrt(b.jcgs/b.mov)*jcgs.cc),
                       (mov.true/rho.mov)/(jcgs.true/rho.jcgs), 
                       (mov.true/rho.mov)/(sqrt(n.jcgs/n.mov)*jcgs.true/rho.jcgs),
                       n = n.mov*(n.mov/n.jcgs),
                       b = b.mov)
# 0.1729660 0.3431789  

#===============================================================================
# jasa vs emi
subsampling_CI.numeric(jasa.cc/(sqrt(b.emi/b.jasa)*emi.cc),
                       (jasa.true/rho.jasa)/(emi.true/rho.emi), 
                       (jasa.true/rho.jasa)/(sqrt(n.emi/n.jasa)*emi.true/rho.emi),
                       n = n.jasa *(n.jasa/n.emi),
                       b = b.jasa)
# 0.2482379 0.6460953

#===============================================================================
# mov vs jasa
subsampling_CI.numeric(mov.cc/(sqrt(b.jasa/b.mov)*jasa.cc),
                       (mov.true/rho.mov)/(jasa.true/rho.jasa), 
                       (mov.true/rho.mov)/(sqrt(n.jasa/n.mov)*jasa.true/rho.jasa),
                       n = n.mov *(n.mov/n.jasa),
                       b = b.mov)
# 0.1126088 0.5626465  

#===============================================================================
# mov vs emi
subsampling_CI.numeric(mov.cc/(sqrt(b.emi/b.mov)*emi.cc),
                       (mov.true/rho.mov)/(emi.true/rho.emi), 
                       (mov.true/rho.mov)/(sqrt(n.emi/n.mov)*emi.true/rho.emi),
                       n = n.mov*(n.mov/n.emi),
                       b = b.mov) 
# 0.2012835 0.2693520   
#===============================================================================

save.image("real-data-graphon_comparison-clus-coeff.RData")

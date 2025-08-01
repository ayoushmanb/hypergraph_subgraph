source("~/R/hyp-real-data/real-data-final/real-data-sparse-graphon-two-stars-supp.R")

#===============================================================================
# JASA
# load the data set
bin.A.jasa <- binary_adj(hyp_edge_list2)
n.jasa <- nrow(bin.A.jasa)
b.jasa <- ceiling(n.jasa/log(n.jasa))
rho.jasa <- edge_density(graph_from_adjacency_matrix(bin.A.jasa)) # 0.001291248
jasa.true <- twostars_fun_scaled(bin.A.jasa) # 4.664715

# subsampling
jasa.cc <- vertex_subsampling(bin.A.jasa, b.jasa, 1000,
                              twostars_fun_scaled, 1)

# # CI using subsampling
# sub.cc.jasa <- subsampling_CI.numeric(jasa.cc, 
#                                       jasa.true, jasa.true, 
#                                       n = n.jasa,
#                                       b = b.jasa)


#===============================================================================
# JCGS
# load the data set
bin.A.jcgs <- binary_adj(hyp_edge_list2)
n.jcgs <- nrow(bin.A.jcgs)
b.jcgs <- ceiling(n.jcgs/log(n.jcgs))
rho.jcgs <- edge_density(graph_from_adjacency_matrix(bin.A.jcgs)) # 0.001350759
jcgs.true <- twostars_fun_scaled(bin.A.jcgs) # 3.858501

# subsampling
jcgs.cc <- vertex_subsampling(bin.A.jcgs, b.jcgs, 1000, twostars_fun_scaled, 1)

# # CI using subsampling
# sub.cc.jcgs <- subsampling_CI.numeric(jcgs.cc,
#                                       jcgs.true, jcgs.true,
#                                       n = n.jcgs,
#                                       b = b.jcgs)


#===============================================================================
# EMI
# load the data set
bin.A.emi <- binary_adj(hyp_edge_list2)
n.emi <- nrow(bin.A.emi)
b.emi <- ceiling(n.emi/log(n.emi))
rho.emi <- edge_density(graph_from_adjacency_matrix(bin.A.emi)) # 0.001014454
emi.true <- twostars_fun_scaled(bin.A.emi) # 5.325487

# subsampling
emi.cc <- vertex_subsampling(bin.A.emi, b.emi, 1000, twostars_fun_scaled, 1)

# CI using subsampling
# sub.cc.emi <- subsampling_CI.numeric(emi.cc,
#                                      emi.true, emi.true, 
#                                      n = n.emi , 
#                                      b = b.emi)

#===============================================================================
# Movie
# load the data set
bin.A.mov <- binary_adj(hyp_edge_list2)
n.mov <- nrow(bin.A.mov)
b.mov <- ceiling(n.mov/log(n.mov))
rho.mov <- edge_density(graph_from_adjacency_matrix(bin.A.mov)) # 0.001812831
mov.true <- twostars_fun_scaled(bin.A.mov) # 7.288776

# subsampling
mov.cc <- vertex_subsampling(bin.A.mov, b.mov, 1000, twostars_fun_scaled, 1)

# CI using subsampling
# sub.cc.mov <- subsampling_CI.numeric(mov.cc,
#                                      mov.true, mov.true, 
#                                      n = n.mov ,
#                                      b = b.mov)




#===============================================================================
# JCGS vs JASA

# subsampling_CI.numeric(jasa.cc/(sqrt(b.jcgs/b.jasa)*jcgs.cc),
#                        (jasa.true)/(jcgs.true), 
#                        (jasa.true)/(sqrt(n.jcgs/n.jasa)*jcgs.true),
#                        n = n.jasa *(n.jasa/n.jcgs),
#                        b = b.jasa) 
# # 0.2555828 1.4908742

subsampling_CI.numeric(1/(jasa.cc/(sqrt(b.jcgs/b.jasa)*jcgs.cc)),
                       1/((jasa.true)/(jcgs.true)), 
                       1/((jasa.true)/(sqrt(n.jcgs/n.jasa)*jcgs.true)),
                       n = n.jcgs *(n.jcgs/n.jasa),
                       b = b.jcgs) 
# 0.2742979 1.0330615
#===============================================================================
# JCGS vs EMI
subsampling_CI.numeric(jcgs.cc/(sqrt(b.emi/b.jcgs)*emi.cc),
                       (jcgs.true)/(emi.true), 
                       (jcgs.true)/(sqrt(n.emi/n.jcgs)*emi.true),
                       n = n.jcgs*(n.jcgs/n.emi),
                       b = b.jcgs) 
# 0.4537878 0.8816997

#===============================================================================
# mov vs jcgs
subsampling_CI.numeric(mov.cc/(sqrt(b.jcgs/b.mov)*jcgs.cc),
                       (mov.true)/(jcgs.true), 
                       (mov.true)/(sqrt(n.jcgs/n.mov)*jcgs.true),
                       n = n.mov*(n.mov/n.jcgs),
                       b = b.mov)
# 0.9025169 2.1993422 

#===============================================================================
# jasa vs emi
subsampling_CI.numeric(jasa.cc/(sqrt(b.emi/b.jasa)*emi.cc),
                       (jasa.true)/(emi.true), 
                       (jasa.true)/(sqrt(n.emi/n.jasa)*emi.true),
                       n = n.jasa *(n.jasa/n.emi),
                       b = b.jasa)
# 0.585716 1.044672 

#===============================================================================
# mov vs jasa
subsampling_CI.numeric(mov.cc/(sqrt(b.jasa/b.mov)*jasa.cc),
                       (mov.true)/(jasa.true), 
                       (mov.true)/(sqrt(n.jasa/n.mov)*jasa.true),
                       n = n.mov *(n.mov/n.jasa),
                       b = b.mov)
# 0.8442397 1.8056675 

#===============================================================================
# mov vs emi
subsampling_CI.numeric(mov.cc/(sqrt(b.emi/b.mov)*emi.cc),
                       (mov.true)/(emi.true), 
                       (mov.true)/(sqrt(n.emi/n.mov)*emi.true),
                       n = n.mov*(n.mov/n.emi),
                       b = b.mov) 
# 1.226910 1.453314 
#===============================================================================

save.image("real-data-graphon_comparison_twostar.RData")

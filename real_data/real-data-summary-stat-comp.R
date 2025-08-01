source("summary-stat-supp.R")

hyp_edge_list <- hyp_edge_list # hypergraph as list

m <- no_hypedge(hyp_edge_list)
n <- no_vertex(hyp_edge_list)

type2_twostars(hyp_edge_list, choose(m,2))

type2_clustering_coefficient(hyp_edge_list, choose(m,2))

bin_twostars_scaled(hyp_edge_list)

bin_clustering_coefficient_scaled(hyp_edge_list)


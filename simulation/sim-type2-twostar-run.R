# setwd("~/hyperedge/sim_1000_1000/")
source("sim-type2-twostar-supp.R")

# Simulation setting ===========================================================

m <- 100 # number of hyp edges
n <- 100 # number of vertices
MC.rep <- 10 # MC iteration to approximate true variance
apx_itr <-  round(m^1.1) # approximate iteration for incomplete U-stat
d <- 5 # degree filtering, default : d = 1 (no filtering)

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
type2.twostar.true <- unlist(mclapply(1:MC.rep, function(rep){
  get.val1(n,m,n.prob, apx_itr, 1)
}, mc.cores =  n.core))
Sys.time()-ST

# Subsampling ==================================================================

Sys.time()
ST <- Sys.time()
type2.twostar.sub <- mclapply(1:sub.MC.rep, function(sub_rep){
  return(get.val2(n, m, n.prob, sub.rep, s.m, apx_itr, apx_itr_sub, d))
}, mc.cores = n.core)
Sys.time()-ST
Sys.time()

# save.image("file_name.RData") 

# Results ======================================================================

type2.twostar.sub.ci(type2.twostar.true, type2.twostar.sub,
                     m, apx_itr, apx_itr_sub, sub.MC.rep, alpha = 0.05)

# Density plot =================================================================

true_val <- as.vector(scale(unlist(type2.twostar.true)))
sub_res <- as.vector(scale(unlist(type2.twostar.sub)[-seq(from = 1, 
                                                          to = sub.MC.rep*(sub.rep+1), 
                                                          by = (sub.rep+1))]))


data <- data.frame(
  value = c(true_val, sub_res),
  group = rep(c("True", "Subsampling"), times = c(length(true_val), length(sub_res)))
)
# remove(true_val, sub_res)

library(ggplot2)

# pdf(paste0("twostar_500_",const,".pdf"))
ggplot(data, aes(x = value, fill = group, color = group)) +
  geom_density(aes(y = ..density..), alpha = 0.2, position = "identity", 
               size = 1) + 
  stat_function(fun = function(x) dnorm(x, mean = 0, sd = 1),
                color = "black", size = 1) +
  scale_fill_manual(values = c("blue", "red")) +  # Custom fill colors
  scale_color_manual(values = c("blue", "red")) +  # Matching density curve colors
  theme_bw() +
  labs(title = "", x = "Value", y = "Density") +
  theme(axis.text.x = element_text(size = 20),  # Increase x-axis font size
        axis.text.y = element_text(size = 20, angle = 90),  # Increase y-axis font size
        axis.title.x = element_text(size = 20), # Increase x-axis title size
        axis.title.y = element_text(size = 20),
        legend.position = "none") # Increase y-axis title sizelegend.position = "none"
# dev.off()
# remove(data)
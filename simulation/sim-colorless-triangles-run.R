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

# Density plot =================================================================

true_val <- as.vector(scale(unlist(colorless.tri.true)))
sub_res <- as.vector(scale(unlist(colorless.tri.sub)[-seq(from = 1, to = sub.MC.rep*(sub.rep+1), by = (sub.rep+1))]))

data <- data.frame(
  value = c(true_val, sub_res),
  group = rep(c("True", "Subsampling"), times = c(length(true_val), length(sub_res)))
)
# remove(true_val, sub_res)

library(ggplot2)

# pdf("file_name.pdf")
ggplot(data, aes(x = value, fill = group, color = group)) +
  geom_density(aes(y = ..density..), alpha = 0.2, position = "identity") + 
  stat_function(fun = function(x) dnorm(x, mean = 0, sd = 1),
                color = "black", size = 0.5) +
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


### ------------------------------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("viridis")
library("parallel")
library("ggpubr")
library("ggthemes")
library("marmap")
library("FactoMineR")
library("ecodist")
library("vegan")

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### Perform Mantel tests between env. dissimilarity and FD indices to assess how env. dissim correlates with dissimilarity in FD indices
### Can env. dissim. predict FD in marine copepods ? 

setwd("/net/kryo/work/fabioben/GODLY/data") 
fd <- get(load("table_mean_ann_FD_indices_baseline+BCP+env_21.11.23.RData"))

### Requires 2 distance matrices: one based on the div indices (one index or several together) + one based on the PCs issued from the PCA
### Because we do not want to subjectively choose the variables to describe env. dissimilarity, we summarize this dissimilarity through PCs
vars <- colnames(fd)[c(1:16,32,39:56)]; vars
data4pca <- na.omit(fd[,vars])
data4pca$NPPv2 <- log10(data4pca$NPPv2)
pca <- PCA(X = data4pca[,c(17:length(data4pca))], scale.unit = T, graph = F, ncp = 7)
# summary(pca)
#                        Dim.1   Dim.2   Dim.3   Dim.4   Dim.5   Dim.6   Dim.7
# Variance               8.192   2.692   2.244   1.586   0.894   0.865   0.749
# % of var.             43.114  14.167  11.808   8.349   4.705   4.554   3.940
# Cumulative % of var.  43.114  57.281  69.089  77.439  82.144  86.697  90.637

### Let's retain the first 6-7 PCs to account for 87-90% of the dissimilarity in env. covariates
data4pca[,c("PC1","PC2","PC3","PC4","PC5","PC6")] <- pca$ind$coord[,c(1:6)]
rm(pca,vars); gc()

### Perform Mantel tests FD index per index - store results on kryo
# https://stats.oarc.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/ 
# These are the two matrices which the function will be testing for a correlation. The test consists of calculating the correlation of the entries in the matrices, then permuting the matrices and calculating the same test statistic under each permutation and comparing the original test statistic to the distribution of test statistics from the permutations to generate a p-value. The number of permutations defines the precision with which the p-value can be calculated. The function to perform the Mantel test is mantel.rtest and the required arguments are the two distance matrices. The number of permutations can also be specified by the user, but is by default 99.
# Mantel r values can fall within a range between -1 to 1

# Wiht ecodist::mantel(), If only one independent variable is given, the simple Mantel r (r12) is calculated. If more than one independent variable is given, the partial Mantel r (ryx|x1 ...) is calculated by permuting one of the original dissimilarity matrices. The bootstrapping is actually resampling without replacement, because duplication of samples is not useful in a dissimilarity context (the dissimilarity of a sample with itself is zero). Resampling within dissimilarity values is inappropriate, just as for permutation. 

# Create the reference env.dist matrix (will be used in all tets below)
env.dist <- dist(data4pca[,c("PC1","PC2","PC3","PC4","PC5","PC6")], method = "euclidean")
vars <- colnames(data4pca)[c(4:9,11,12,14,16)]; vars
# v <- "SR"

mclapply(vars, function(v) {
    
        message(paste("Performing Mantel test for ",v, sep = ""))
    
        fdiv.dist <- dist(data4pca[,v], method = "euclidean")
    
        mantel.test <- ecodist::mantel(fdiv.dist ~ env.dist, nperm = 99, mrank = T)
        
        setwd("/net/kryo/work/fabioben/GODLY/data/mantel_tests/")
        save(x = mantel.test, file = paste("mantel_test_result_",v,"_23_11_23.RData", sep = ""))
    
    }, mc.cores = length(vars)

) # eo mclapply - v in vars

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
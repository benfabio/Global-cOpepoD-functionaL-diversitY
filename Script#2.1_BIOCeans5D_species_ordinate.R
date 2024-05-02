
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 21/07/23: R script to create ordination plots showing the position of copepod species  © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - Compute Gower distance matrix + PCoA to plot species in multidimensional trait space (the one used for FD indices from dbFD())
# - Compute FAMD (+ Euclidean dist) to plot species in multidimensional trait space (the one from Benedetti et al., 2023)
# - Compute AUC criterion from Mouillot et al. (2021) to evaluate quality of trait space
# - Compute AUC score for increasing nb of max PCoA dimensions and FAMD dimensions (2 plots sensu Mouillot et al. 2021)

### Latest update: 24/07/23

### ------------------------------------------------------------------------------------------------------------------------------------------------------

# install.packages("coRanking")
library("gawdis")
library("FD")
library("ape")
library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("viridis")
library("ggpubr")
library("parallel")
library("xlsx")
library("readxl")
library("lubridate")
library("worrms")
library("clValid")
library("flashClust")

world <- map_data("world") # coastlines for maps

setwd("/net/kryo/work/fabioben/BIOCeans5D/data") # working dir

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 1°) Load and format functional traits table 
traits <- read.csv("traits_table_Benedetti2023.csv", h = T, sep = ";", dec = ",")
# dim(traits); str(traits)
# colnames(traits)
colnames(traits)[8] <- "Body.length" # size vector that we will keep 
# summary(traits)
# unique(traits$Trophic.group)

# Replace "" by NA
library("naniar")
traits <- traits %>% replace_with_na_all(condition = ~.x == "")
# Convert Feeding.mode, Trophic.group and Spawning.mode to factors
traits$Spawning.mode <- as.factor(traits$Spawning.mode)
traits$Trophic.group <- as.factor(traits$Trophic.group)
traits$Feeding.mode <- as.factor(traits$Feeding.mode)
# Count NA
names <- colnames(traits)[c(8:11,16)] ; names
traits$na_count <- apply(traits[,names], 1, function(x) sum(is.na(x)))
# Drop species with missing body length and more than two missing traits
traits_red <- traits[!is.na(traits$Body.length),]
# dim(traits_red) # 384/385 (retained all spp but 1)
traits_red2 <- traits_red[traits_red$na_count < 2,]
# dim(traits_red2) # 355/385 (retained 92.2% of the initial pool of species)
# Convert to 1 and 0 for Gower (as to be factors for FAMD+Eucli though)
traits_red2$Myelination <- as.integer(as.logical(traits_red2$Myelination))
traits_red2$Omnivore <- as.integer(as.logical(traits_red2$Omnivore))
traits_red2$Carnivore <- as.integer(as.logical(traits_red2$Carnivore))
traits_red2$Herbivore <- as.integer(as.logical(traits_red2$Herbivore))
traits_red2$Detritivore <- as.integer(as.logical(traits_red2$Detritivore))
traits_red2$Current <- as.integer(as.logical(traits_red2$Current))
traits_red2$Cruise <- as.integer(as.logical(traits_red2$Cruise))
traits_red2$Ambush <- as.integer(as.logical(traits_red2$Ambush))
traits_red2 <- data.frame(traits_red2)

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 2°) Compute Gower's distance matrix
# Compute Gower's distance matrix 
# ?gawdis
# Why gawdis is better than gowdis: https://cran.r-project.org/web/packages/gawdis/gawdis.pdf 
# gawdis provides a solution to the problem of unequal traits contribution when combining different traits in a multi-trait dissimilarity (without care the
# contribution of some traits can much stronger than others, i.e. the correlation of the dissimilarity of individual trait, with the multi-trait dissimilarity, will be much stronger for some traits, particularly categorical ones)
# Other advantage of gawdis: you can specify groups of variables that specify the same traits and that are encoded in a fuzzy way
gow <- gawdis(x = traits_red2[,c(8,9,10,12:15,17:19)], groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))

# Project in PCoA
# ?pcoa
# pcoa <- pcoa(D = gow, correction = "cailliez")
# str(pcoa)
# head(pcoa$values)
# biplot(pcoa)

# Or with 'vegan'
library("vegan")
pcoa2 <- wcmdscale(d = gow, eig = T)
# summary(pcoa2) 
# str(pcoa2)
# plot(pcoa2)
pcoa.scores <- data.frame(pcoa2$points)
scores1 <- paste0("PCoA 1 (",floor(pcoa2$eig[1]*100)/100,"%)")
scores2 <- paste0("PCoA 2 (",floor(pcoa2$eig[2]*100)/100,"%)")
scores3 <- paste0("PCoA 3 (",floor(pcoa2$eig[3]*100)/100,"%)")
scores4 <- paste0("PCoA 4 (",floor(pcoa2$eig[4]*100)/100,"%)")
scores5 <- paste0("PCoA 5 (",floor(pcoa2$eig[5]*100)/100,"%)")
scores6 <- paste0("PCoA 6 (",floor(pcoa2$eig[6]*100)/100,"%)")

P1 <- ggplot(aes(x = Dim1, y = Dim2), data = pcoa.scores) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..)), fill = "#7fbc41") + 
    geom_point(colour = "#276419", size = 1) + scale_alpha_continuous(range = c(0, 1)) + 
    theme_bw() + xlab(scores1) + ylab(scores2)+ guides(alpha = "none")

P2 <- ggplot(aes(x = Dim3, y = Dim4), data = pcoa.scores) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..)), fill = "#7fbc41") + 
    geom_point(colour = "#276419", size = 1) + scale_alpha_continuous(range = c(0, 1)) + 
    theme_bw() + xlab(scores3) + ylab(scores4)+ guides(alpha = "none")
 
P3 <- ggplot(aes(x = Dim5, y = Dim6), data = pcoa.scores) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..)), fill = "#7fbc41") + 
    geom_point(colour = "#276419", size = 1) + scale_alpha_continuous(range = c(0, 1)) + 
    theme_bw() + xlab(scores5) + ylab(scores6)+ guides(alpha = "none")
    
# Save
panel <- ggarrange(P1,P2,P3, labels = letters, align = "hv", ncol = 2, nrow = 2)

 
### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 3°) Compute FAMD + Euclid distance matrix
library("FactoMineR")
library("missMDA")

# Estimate ncp for fct
npfamd <- estim_ncpFAMD(traits_red2[,c(8,9,10,12:15,17:19)], ncp.min = 0, ncp.max = 6)

# Perform FAMD for data with missing values
compfamd <- imputeFAMD(traits_red2[,c(8,9,10,12:15,17:19)], ncp = npfamd$ncp)
FAMD <- FAMD(traits_red2[,c(8,9,10,12:15,17:19)], tab.disj = compfamd$tab.disj, graph = F, ncp = npfamd$ncp)
# summary(FAMD) 

### Compare both FAMDs by aggregating loadings in one table that you'll add to Appendix S4
eigFAMD <- data.frame(perc = FAMD$eig[,"percentage of variance"], nb = c(1:nrow(FAMD$eig)) )
# sum(eigFAMD$perc)
# % of explained variance
famd1 <- paste0("FAMD 1 (",floor(eigFAMD$perc[1]*100)/100,"%)")
famd2 <- paste0("FAMD 2 (",floor(eigFAMD$perc[2]*100)/100,"%)")
famd3 <- paste0("FAMD 3 (",floor(eigFAMD$perc[3]*100)/100,"%)")
famd4 <- paste0("FAMD 4 (",floor(eigFAMD$perc[4]*100)/100,"%)")

# Extract scores and plot FAMD space
# str(FAMD)
# FAMD$var

# Bind in table and plo
ind.coords <- data.frame(FAMD$ind$coord)
var.coords <- data.frame(FAMD$var$coord)

library("ggrepel")
p1 <- ggplot() + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..), x = Dim.1, y = Dim.2), fill = "#7fbc41", data = ind.coords) + 
    geom_point(aes(x = Dim.1, y = Dim.2), colour = "#276419", size = 1, data = ind.coords) + 
    geom_point(aes(x = Dim.1*3, y = Dim.2*3), data = var.coords, colour = "black", fill = "#c51b7d", pch = 21, size = 2) +    
    geom_text_repel(aes(x = Dim.1*3, y = Dim.2*3), data = var.coords, colour = "#c51b7d", label = rownames(var.coords)) +  
    scale_alpha_continuous(range = c(0,1)) + 
    xlab(famd1) + ylab(famd2) + theme_bw() + guides(alpha = "none")

p2 <- ggplot() + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..), x = Dim.3, y = Dim.4), fill = "#7fbc41", data = ind.coords) + 
    geom_point(aes(x = Dim.3, y = Dim.4), colour = "#276419", size = 1, data = ind.coords) + 
    geom_point(aes(x = Dim.3*3, y = Dim.4*3), data = var.coords, colour = "black", fill = "#c51b7d", pch = 21, size = 2) +    
    geom_text_repel(aes(x = Dim.3*3, y = Dim.4*3), data = var.coords, colour = "#c51b7d", label = rownames(var.coords)) +  
    scale_alpha_continuous(range = c(0,1)) + 
    xlab(famd3) + ylab(famd4) + theme_bw() + guides(alpha = "none")
    
# Plot panle    
ggarrange(p1,p2, labels = letters, align = "hv", ncol = 2, nrow = 1)

# To derive Euclid. distance matrix
famd_dist <- dist(ind.coords, method = "euclidean")
# clust <- flashClust::hclust(d = famd_dist, method = "ward")
# plot(clust)
    
    
### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 4°) Compute AUC criterion - should be above 0.7 according to Mouillot et al., 2021 (Ecol. Lett.)
library("dimRed")
library("coRanking")

### Calculate coranking matrix for the FAMD-based euclid distance matrix
Q = coranking(Xi = gow, X = famd_dist, input_Xi = c("dist"), input_X = c("dist"), use = "R")
coRanking::imageplot(Q)
NX <- coRanking::R_NX(Q)
AUC <- coRanking::AUC_ln_K(NX)
#round(AUC,3) # 0.750 ! Good


### To evaluate whethrr keeping 4-5 PCoA axis retains enough information compared to the full Gower matrix
euclid_dist <- dist(pcoa.scores[,1:5], method = "euclidean")

Q = coranking(Xi = gow, X = euclid_dist, input_Xi = c("dist"), input_X = c("dist"), use = "R")
coRanking::imageplot(Q)
NX <- coRanking::R_NX(Q)
AUC <- coRanking::AUC_ln_K(NX)
#round(AUC,3)
# 0.834 for 5 dimensions
# 0.822 for 4 
# 0.786 for 3
# 0.749 for 2 though
### --> stick to 4 dimensions of PCoA when computing dbFD based on Gower 


### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 5°) Plot AUC score for increasing nb of max PCoA dimensions and FAMD dimensions 
library("dimRed")
library("coRanking")
library("vegan")

### 5.A) With PCoA axes
pcoa2 <- wcmdscale(d = gow, eig = T)
pcoa.scores <- data.frame(pcoa2$points)
# dim(pcoa.scores) # 355 spp x 106 pcoa dimensions

res <- mclapply(c(2:35), function(i) {
    
        euclid_dist <- dist(pcoa.scores[,1:i], method = "euclidean")
        Q <- coranking(Xi = gow, X = euclid_dist, input_Xi = c("dist"), input_X = c("dist"), use = "R")
        #coRanking::imageplot(Q)
        NX <- coRanking::R_NX(Q)
        AUC <- coRanking::AUC_ln_K(NX)
        s <- round(AUC,4)
        # Return i and s
        return( data.frame(Ndim = i, AUC = s) )
    
    }, mc.cores = 20

) # eo mclapply
# Rbind
tab1 <- dplyr::bind_rows(res)
# tab

# ggplot(aes(x = Ndim, y = AUC), data = tab) + geom_point() +
#     xlab("Number of dimensions (PCoA axes)") + 
#     ylab("Quality of species trait space (AUC)") +
#     theme_bw()
### --> 4 or 5 more than enough


### 5.B) With FAMD axes
library("FactoMineR")
library("missMDA")

npfamd <- estim_ncpFAMD(traits_red2[,c(8,9,10,12:15,17:19)], ncp.min = 0, ncp.max = 9)
compfamd <- imputeFAMD(traits_red2[,c(8,9,10,12:15,17:19)], ncp = 9)
FAMD <- FAMD(traits_red2[,c(8,9,10,12:15,17:19)], tab.disj = compfamd$tab.disj, graph = F, ncp = 9)
FAMD.scores <- data.frame(FAMD$ind$coord)
# dim(FAMD.scores)

# Compute AUC for N dims 2:9
res <- mclapply(c(2:9), function(i) {
    
        euclid_dist <- dist(FAMD.scores[,1:i], method = "euclidean")
        Q <- coranking(Xi = gow, X = euclid_dist, input_Xi = c("dist"), input_X = c("dist"), use = "R")
        NX <- coRanking::R_NX(Q)
        AUC <- coRanking::AUC_ln_K(NX)
        s <- round(AUC,4)
        # Return i and s
        return( data.frame(Ndim = i, AUC = s) )
    
    }, mc.cores = 20

) # eo mclapply
# Rbind
tab2 <- dplyr::bind_rows(res)
# tab

# ggplot(aes(x = Ndim, y = AUC), data = tab) + geom_point() +
#     xlab("Number of dimensions (FAMD axes)") + 
#     ylab("Quality of species trait space (AUC)") +
#     theme_bw()

### Rbind both tabs and plot together 
tab1$method <- "PCoA"
tab2$method <- "FAMD"
tabs <- rbind(tab1,tab2)

ggplot(aes(x = Ndim, y = AUC, colour = factor(method)), data = tabs) + geom_point() +
    xlab("Number of dimensions (FAMD axes)") + ylab("Quality of species trait space (AUC)") +
    geom_vline(xintercept = 4, linetype = "dashed") + theme_bw()

### Use PCoA + Euclidean distance based on 4 axes :-)  

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
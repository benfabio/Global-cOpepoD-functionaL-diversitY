
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 21/07/23: R script to compute the functional distinctiveness index (IntDi) for the copepod species © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - compute the functional distinctiveness index (IntDi) of Coulon et al., (2023) - doi:10.1111/geb.13731
# - plot profile of integrated IntDi per species (Fig. 2 of Coulon et al., (2023))
# - explore covariance between IntDi and traits (continuous and categorical)

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
library("gtools")

world <- map_data("world") # coastlines for maps

setwd("/net/kryo/work/fabioben/GODLY/data") # working dir

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 1°) Load and format functional traits table 
traits <- read.csv("traits_table_Benedetti2023.csv", h = T, sep = ";", dec = ",")
colnames(traits)[8] <- "Body.length" # size vector that we will keep 
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
traits_red2 <- traits_red[traits_red$na_count < 2,]
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


### 2°) Compute IntDi based on the R code made available at: https://figshare.com/articles/software/coulon_et_al_2023_IntDi_calculation/22317643/1

# Some documentation on the method below, from Coulon et al., 2023
# The functional distinctiveness index (IntDi) that not only considers one combination of traits (i.e. the combination of all available traits) but all possible combinations of available traits. It ensures that a species has a high distinctiveness value not because of a single extreme trait value but because it has several uncommon trait value. 
# We calculated the dissimilarity matrix (i.e. dissimilarities between species based on their traits) of each possible combination of traits among the 11 traits selected for this study, considering that a minimum of 4 is needed to characterize the difference between species (Petchey & Gaston, 2002). This procedure therefore provides a total of N dissimilarity matrices. We then calculated the IntDi value of each species from the integrated dissimilarity matrix itself computed as the average between the N matrices. Additionally, to prevent any disproportionate contribution of categorical/qualitative traits in dissimilarity matrices and thus on distinctiveness computation, we applied the approach developed by de Bello et al. (2021) using the ‘gawdis’ package. Functional distinctiveness was standardized so that it ranges between 0 and 1. We followed the same procedure for the resulting trait database, the database with no missing data, and the database with missing data predicted with the R package ‘missForest’ (Stekhoven & Buehlmann, 2012).

## Creating the full list of traits combinations
# colnames(traits_red2) # 
traits.names <- colnames(traits_red2)[c(8,9,10,12:15,17:19)] ; traits.names # 10 traits
tab <- traits_red2[,traits.names] # traits table to consider for defining all possible traits combin
rownames(tab) <- traits_red2$Species
# dim(tab) # should be 355 spp x 10 traits dimensions
n_traits_min <- 4 + 1 # You must select at least 5 traits to make combinations of 4 traits
n_traits <- ncol(tab) - 1
combi <- as.data.frame(t(combn(unique(n_traits), n_traits))) # empty ddf where all traits combin will be stored

for(i in n_traits_min:n_traits-1) {
  # i <- 6
  combi <- bind_rows(combi, as.data.frame(combinations(n = n_traits, r = i, v = seq(1, n_traits, 1), repeats.allowed = F)))
} # eo for loop
# nrow(combi) # 382 possible traits computations 

# Get dissimilarity matrices
# list_mat <- list()
# i <- 300 # for testing gawdis

### Run into a mclapply (long otherwise)
list_mat <- mclapply( c(1:nrow(combi)), function(i) {
    
        data <- tab
        codex <- as.vector(t(combi[i,]))
        data <- dplyr::select(data, all_of(na.omit(codex)))
        # dim(data) # subsetted
        gaw <- gawdis(data, w.type = "optimized", opti.maxiter = 300) # decrease opti.maxiter for shorter comput times
        mat_gaw <- as.matrix(gaw)
        # list_mat <- c(list_mat, list(mat_gaw))
        cat(paste("----", Sys.time(), "----", i, "/", nrow(combi), "combination(s) selected ----", sep = " "))
        return(mat_gaw)
        
    }, mc.cores = 25
    
) # eo mclapply
# Save output
# class(list_mat)
# str(list_mat)
# length(list_mat) # should be nrow(combi)
# list_mat[[300]] ; str(list_mat[[300]]) ; class(list_mat[[300]])
save(list_mat, file = "list_gawdis_mat_IntDi_24_07_23.Rdata")

### Version from Coulon et al.: long unparallelized for loop
# for (i in c(1:nrow(combi))) {
#     #data <- dplyr::select(traits, -taxon)
#     data <- tab
#     codex <- as.vector(t(combi[i,]))
#     data <- dplyr::select(data, all_of(na.omit(codex)))
#     # dim(data) # subsetted
#     gaw <- gawdis(data, w.type = "optimized", opti.maxiter = 300) # decrease opti.maxiter for shorter comput times
#     mat_gaw <- as.matrix(gaw)
#     list_mat <- c(list_mat, list(mat_gaw))
#     cat(paste("----", Sys.time(), "----", i, "/", nrow(combi), "combination(s) selected ----", sep = " "))
# }

## Compute average matrix
#load("output/list_mat_with_NA.Rdata")
# ?Reduce 
mean_matrix_gaw <- Reduce("+",list_mat)/length(list_mat)
# dim(mean_matrix_gaw)
# str(mean_matrix_gaw)
# summary(mean_matrix_gaw)

## Compute IntDi
IntDi <- apply(mean_matrix_gaw, 1, sum, na.rm = T) / (nrow(mean_matrix_gaw) + 1)
IntDi <- cbind(as.data.frame(IntDi), traits_red2$Species)
colnames(IntDi)[2] <- "Species"
# Examine
# dim(IntDi); head(IntDi); str(IntDi)
# summary(IntDi)
# save(IntDi, "IntDi_gaw_with_NA.Rdata")
#IntDi[IntDi$IntDi < .2,] # Why 0? 
#IntDi[IntDi$IntDi > .4,] # 

IntDi[order(IntDi$Species, decreasing = F),]

### Plot
IntDi$Species <- factor(IntDi$Species, levels = IntDi[order(IntDi$IntDi),"Species"]) # 
ggplot(IntDi[IntDi$IntDi > 0,], aes(Species, IntDi)) + geom_point() + ylab("Integrated distinctiveness (IntDi)") +
    scale_x_discrete(name = "Species", labels = NULL) + theme_classic() + theme(axis.ticks = element_blank()) 


### 3°) Explore relationships between IntDi and traits

### Supply to 'traits_red2' and examine relationship with traits
traits_red2$IntDi <- IntDi$IntDi

### ~ body length
ggplot() + geom_point(aes(x = Body.length, y = IntDi), data = traits_red2[traits_red2$IntDi > 0,]) +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Maximum body length (mm)") +
    theme_classic()
cor(traits_red2$Body.length, traits_red2$IntDi, method = "spearman")
# 0.0328 # NS 


### ~ FG
ggplot(aes(y = IntDi, x = factor(FG)), data = traits_red2[traits_red2$IntDi > 0 & !is.na(traits_red2$FG),]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("FG (Benedetti et al., 2023)") +
    theme_classic()
### --> FGs that show highest IntDiv: 2,5,1,4...then 3,6.
unique(traits_red2[traits_red2$FG == 2,"Species"]) # Oncaeidae (only detritivores...)
unique(traits_red2[traits_red2$FG == 5,"Species"]) # Very large carnivores
unique(traits_red2[traits_red2$FG == 1,"Species"]) # Clausocalanus because of spawning...
unique(traits_red2[traits_red2$FG == 4,"Species"]) # Small ambush carnivores (Corycaeus)
# Most specialized taxa. Makes sense! Not your average copepod
### --> FGs that show LOWEST IntDiv: 9,10,11
unique(traits_red2[traits_red2$FG == 9,"Species"]) # small-medium sized omnivores
unique(traits_red2[traits_red2$FG == 10,"Species"]) # oithona
unique(traits_red2[traits_red2$FG == 11,"Species"]) # Acartia/Centropages


### ~ Order/Families
ggplot(aes(y = IntDi, x = factor(Family)), data = traits_red2[traits_red2$IntDi > 0,]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Family") + theme_classic() + theme(axis.text.x=element_text(angle = -90, hjust = 0))
# --> Euchaeta, Sapphirina, Oncaea

ggplot(aes(y = IntDi, x = factor(Order)), data = traits_red2[traits_red2$IntDi > 0,]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Order") + theme_classic() + theme(axis.text.x=element_text(angle = -90, hjust = 0))
# --> Calanoida LESS distinct than Cyclopoida (as expected)


### ~ Myelination
ggplot(aes(y = IntDi, x = factor(Myelination)), data = traits_red2[traits_red2$IntDi > 0,]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Myelination (1/0)") + theme_classic() 
# n.s. 


### ~ Spawning strat
ggplot(aes(y = IntDi, x = factor(Spawning.mode)), data = traits_red2[traits_red2$IntDi > 0,]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Spawning mode") + theme_classic() 
# --> sac-spawners are more distinct


### ~ Trophic group
ggplot(aes(y = IntDi, x = factor(Trophic.group)), data = traits_red2[traits_red2$IntDi > 0,]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Trophic group") + theme_classic() 
# --> Detritivore and pure carnivores > distinct


### ~ Feeding.mode
ggplot(aes(y = IntDi, x = factor(Feeding.mode)), data = traits_red2[traits_red2$IntDi > 0,]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Feeding mode") + theme_classic() 
# --> cruise-feeders and more distinct


### 4°) ~ 4 PCoA axes! (corr, plots)
mat.gow <- gawdis(x = traits_red2[,c(8,9,10,12:15,17:19)], groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))
library("vegan")
pcoa2 <- wcmdscale(d = mat.gow, eig = T)
pcoa.scores <- data.frame(pcoa2$points)
scores1 <- paste0("PCoA 1 (",floor(pcoa2$eig[1]*100)/100,"%)")
scores2 <- paste0("PCoA 2 (",floor(pcoa2$eig[2]*100)/100,"%)")
scores3 <- paste0("PCoA 3 (",floor(pcoa2$eig[3]*100)/100,"%)")
scores4 <- paste0("PCoA 4 (",floor(pcoa2$eig[4]*100)/100,"%)")

# dim(pcoa.scores)
IntDi$PCoA1 <- pcoa.scores$Dim1
IntDi$PCoA2 <- pcoa.scores$Dim2
IntDi$PCoA3 <- pcoa.scores$Dim3
IntDi$PCoA4 <- pcoa.scores$Dim4
# Test correlations

round(cor(IntDi[IntDi$IntDi > 0,"IntDi"], IntDi[IntDi$IntDi > 0,"PCoA1"], method = "spearman"),3) # -0.279
round(cor(IntDi[IntDi$IntDi > 0,"IntDi"], IntDi[IntDi$IntDi > 0,"PCoA2"], method = "spearman"),3) # 0.669
round(cor(IntDi[IntDi$IntDi > 0,"IntDi"], IntDi[IntDi$IntDi > 0,"PCoA3"], method = "spearman"),3) # -0.053 N.S.
round(cor(IntDi[IntDi$IntDi > 0,"IntDi"], IntDi[IntDi$IntDi > 0,"PCoA4"], method = "spearman"),3) # 0.261
# PCoA 2 seems to be the main axis of variation

ggplot(aes(x = PCoA1, y = IntDi), data = IntDi[IntDi$IntDi > 0,]) + 
    #geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_point(aes(fill = IntDi), pch = 21, colour = "black") +
    scale_fill_viridis(name = "IntDiv", option = "A", direction = -1) + 
    theme_bw() + xlab(scores1) + ylab("IntDi") 
#
ggplot(aes(x = PCoA2, y = IntDi), data = IntDi[IntDi$IntDi > 0,]) + 
    #geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_point(aes(fill = IntDi), pch = 21, colour = "black") +
    geom_smooth(colour = "black", se = T, method = "lm") + 
    scale_fill_viridis(name = "IntDiv", option = "A", direction = -1) + 
    theme_bw() + xlab(scores2) + ylab("IntDi") 
#
ggplot(aes(x = PCoA3, y = IntDi), data = IntDi[IntDi$IntDi > 0,]) + 
    #geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_point(aes(fill = IntDi), pch = 21, colour = "black") +
    #geom_smooth(colour = "black", se = T, method = "lm") + 
    scale_fill_viridis(name = "IntDiv", option = "A", direction = -1) + 
    theme_bw() + xlab(scores3) + ylab("IntDi") 

# Plot in PCoA space
p1 <- ggplot(aes(x = PCoA1, y = PCoA2), data = IntDi[IntDi$IntDi > 0,]) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..)), fill = "grey") + 
    geom_point(aes(fill = IntDi), pch = 21, colour = "black") +
    scale_fill_viridis(name = "IntDiv", option = "A", direction = -1) + 
    scale_alpha_continuous(range = c(0, 1)) + 
    theme_bw() + xlab(scores1) + ylab(scores2)+ guides(alpha = "none")

p2 <- ggplot(aes(x = PCoA3, y = PCoA4), data = IntDi[IntDi$IntDi > 0,]) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..)), fill = "grey") + 
    geom_point(aes(fill = IntDi), pch = 21, colour = "black") + 
    scale_fill_viridis(name = "IntDiv", option = "A", direction = -1) + 
    scale_alpha_continuous(range = c(0, 1)) + 
    theme_bw() + xlab(scores3) + ylab(scores4)+ guides(alpha = "none")

ggarrange(p1,p2, labels = letters, align = "hv", ncol = 2, nrow = 1)

### Show taxa with highest PCoA2 scores
# IntDi[order(IntDi$PCoA2, decreasing = T),]

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------

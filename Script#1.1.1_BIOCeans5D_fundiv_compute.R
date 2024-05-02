
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 20/07/23: R script to compute the multidimensional FD indices ('FD' package) from the community tables made during the MSc thesis of Jonas Wydler (data from Benedetti et al., 2023 - JBIO) © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - implement functions from the 'FD' R package to compute several FD indices ('dbFD' function)
# - test functions for one SDM (GAM) and one monthly community table (contemp conditions)
# - implement for every SDM (GLM/GAM/ANN) and months and ESMs for future projections; compute uncertainties (stdev)

### Latest update: 27/07/23

### ------------------------------------------------------------------------------------------------------------------------------------------------------

# install.packages("biomod2")
# library("fundiversity")
library("marmap")
library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("viridis")
library("gawdis")
library("FD")
library("parallel")
library("xlsx")
library("readxl")
library("flashClust")
library("naniar")

world <- map_data("world") # coastlines for maps

setwd("/net/kryo/work/fabioben/BIOCeans5D/data") # working dir on kryo

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 1°) Compute multidimensional FD indices 

### 1.A) Load and format functional traits table 
traits <- read.csv("traits_table_Benedetti2023.csv", h = T, sep = ";", dec = ",")
# colnames(traits)
colnames(traits)[8] <- "Body.length" # size vector that we will keep 
# unique(traits$Trophic.group)
# Replace "" by NA
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

### Compute Gower's distance matrix
# ?gawdis
# Why gawdis is better than gowdis: https://cran.r-project.org/web/packages/gawdis/gawdis.pdf 
# gawdis provides a solution to the problem of unequal traits contribution when combining different traits in a multi-trait dissimilarity (without care the
# contribution of some traits can much stronger than others, i.e. the correlation of the dissimilarity of individual trait, with the multi-trait dissimilarity, will be much stronger for some traits, particularly categorical ones)
# Other advantage of gawdis: you can specify groups of variables that specify the same traits and that are encoded in a fuzzy way
gow <- gawdis(x = traits_red2[,c(8,9,10,12:15,17:19)], groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))

# plot(gow)
# str(gow)
# Derive dendrogram
#clust <- flashClust::hclust(d = gow, method = "ward")
#plot(clust)
# Cut tree and check FG composition
#traits_red2$FG2 <- cutree(clust, 10) # 
#unique(traits_red2[traits_red2$FG2 == 6,"Species"])
# Good. Makes sense with gawdis() instead of gowdis as well. Matches groups from Benedetti et al., 2023 when based on Ward
# You can use this Gower dist matrix to project a multdim functional space with a PCoA


### 1.B) Load a standard community table and compute FD indices based on 'FD' functions for each grid cell
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/community_tables_05_01_21/contemp")
comm <- read.table("table_mon_composition_baseline_GLM_mar.txt")
# dim(comm) # 41553   306 (less species that in traits table)
# summary(comm)
### --> Missing values = missing predictors; 0 are to keep
# Subset traits_red2 to spp of interest
spp2keep <- colnames(comm)[c(4:length(comm))] #; spp2keep
# Check if they match in "traits_red2"
commons <- intersect(spp2keep,traits_red2$Species)# ; length(commons) # FD indices based on 284 taxa
# setdiff(spp2keep, traits_red2$Species) # identifies the species that are in the comm table (species modelled) BUT NOT in the traits table
#  [1] "Atrophia_glacialis"          "Clytemnestra_scutellata"
#  [3] "Epilabidocera_amphitrites"   "Goniopsyllus_rostratus"
#  [5] "Jaschnovia_tolli"            "Lubbockia_squillimana"
#  [7] "Miracia_efferata"            "Mormonilla_phasma"
#  [9] "Neomormonilla_minor"         "Parapontella_brevicornis"
# [11] "Parvocalanus_crassirostris"  "Phaenna_spinifera"
# [13] "Pseudoamallothrix_cenotelis" "Pseudoamallothrix_ovata"
# [15] "Scolecithricella_orientalis" "Scolecitrichopsis_ctenopus"
# [17] "Stephos_longipes"            "Temorites_brevis"
# [19] "Temoropia_mayumbaensis"
 
comm_fdiv <- na.omit(comm[,c(commons)])

### Subset both tables based on commons
# colnames(traits_red2)
traits_fdiv <- traits_red2[traits_red2$Species %in% commons,c(8,9,10,12:15,17:19)]
rownames(traits_fdiv) <- traits_red2[traits_red2$Species %in% commons,"Species"]
gow_fdiv <- gawdis(x = traits_fdiv, groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))
# gow2 <- FD::gowdis(x = traits_fdiv)

# Different labels?
# setdiff(labels(gow), commons) # no diff

### AUC criterion test from Mouillot et al., 2021
library("dimRed")
library("coRanking")
library("vegan")
# Calculate coranking matrix for the FAMD-based euclid distance matrix
# To evaluate whether keeping 4-5 PCoA axis retains enough information compared to the full Gower matrix
# Step 1: derive PCoA from full Gower matrix 'gow_fdiv' and retrieve coordinates
# pcoa <- wcmdscale(d = gow_fdiv, eig = T)
# pcoa.scores <- data.frame(pcoa$points) # dim(pcoa.scores) # 88 dimensions...
# # Step 2: derive euclidean dist mat from pcoa
# euclid_dist <- dist(pcoa.scores[,1:4], method = "euclidean")
# # Step 3: compute co-ranking matrix
# Q = coranking(Xi = gow_fdiv, X = euclid_dist, input_Xi = c("dist"), input_X = c("dist"), use = "R")
# coRanking::imageplot(Q)
# NX <- coRanking::R_NX(Q)
# AUC <- coRanking::AUC_ln_K(NX)
# round(AUC,3)
# 0.8332777 for 5 dimensions
# 0.81809 for 4 
# 0.7865 for 3
# 0.743 for 2 though
### --> stick to 4 dimensions of PCoA when computing dbFD based on an euclidean distance matrix (see script 2.1 as well)


### With ape to apply cailliez correction (Villéger & Mouillot etc.)
# pcoa2 <- ape::pcoa(gow_fdiv, correction = "cailliez")
# # str(pcoa2)
# # head(pcoa2$vectors) # use that
# euclid_dist <- dist(pcoa2$vectors[,1:4], method = "euclidean")
# # Step 3: compute co-ranking matrix
# Q = coranking(Xi = gow_fdiv, X = euclid_dist, input_Xi = c("dist"), input_X = c("dist"), use = "R")
# coRanking::imageplot(Q)
# NX <- coRanking::R_NX(Q)
# AUC <- coRanking::AUC_ln_K(NX)
# round(AUC,3)
# 0.832 for 5 dims
# 0.818 for 4
# 0.786 for 3
### --> stick to 4 dimensions of PCoA when computing dbFD based on an euclidean distance matrix 

### Define here the final euclidena distance matrix to use in dbFD()
pcoa2 <- ape::pcoa(gow_fdiv, correction = "cailliez")
euclid_dist <- dist(pcoa2$vectors[,1:4], method = "euclidean") # 4 dimensions from the PCoA based on the Cailliez corr
# is.euclid(euclid_dist)

### Compute FD indices except FRic and FDiv (need to have 1/0 for that)
### Constant error message:
# Erreur dans dbFD(x = euclid_dist, a = dummy, calc.FDiv = F, calc.FRic = F,:
# At least one species does not occur in any community (zero total abundance across all communities).
# Does not like 0...
# Check if they are species thta always have zeroes
#check <- colSums(as.matrix(comm_fdiv), na.rm = F, dims = 1) 
#check[check == 0] # nope
# check if one community has total abundance of zero (no species)
# a <- comm_fdiv[29,]
# abun.sum <- apply(a, 1, sum)
# if ( any(abun.sum == 0) ) stop("At least one community has zero-sum abundances (no species).","\n")
# # check if one species has total abundance of zero (never occurs)
# abun.sum2 <- apply(a, 2, sum)
# if ( any(abun.sum2 == 0) ) stop("At least one species does not occur in any community (zero total abundance across all communities).","\n")
# dummy <- comm_fdiv[333,]
# summary(dummy)
# Replace 0 by very low HSI value (0.000001)
# dummy[dummy == 0] <- 0.000001 # then it works

save(comm_fdiv, file = "comm_fdiv.Rdata")
save(euclid_dist, file = "euclid_dist_postPCoA.Rdata")
save(gow_fdiv, file = "gower_dist.Rdata")


### Tests on one randomly chosen lines
#comm_fdiv_no_zero <- comm_fdiv
#comm_fdiv_no_zero[comm_fdiv_no_zero == 0] <- 0.0001 # dim(comm_fdiv_no_zero)
#system.time(
#    dbFD(x = euclid_dist, a = comm_fdiv_no_zero[1:10,], calc.FDiv = F, calc.FRic = F, w.abun = T, stand.x = F, scale.RaoQ = F, calc.FGR = F, calc.CWM = F)
#) #
# utilisateur     système      écoulé 
#     54.990       9.975      46.268 
### --> 10 rows takes about 50s
### Same if you weight by abundances by the way
### --> Twice longer if your provide the Gower dist matrix instead of the Euclidean one...

### Testing with 1/0 instead of HSI to check if it goes faster
# comm_fdiv_PA <- comm_fdiv
# # ?bm_BinaryTransformation # from biomod2
# #library("biomod2")
# #comm_fdiv_PA <- bm_BinaryTransformation(data = comm_fdiv_PA, threshold = .30)
# # summary(comm_fdiv_PA)
# system.time(
#     dbFD(x = gow_fdiv, a = comm_fdiv_PA[1:10,], m = 4, corr = "cailliez")
# ) #

### Testing with mclapply()
# indices <- dbFD(x = euclid_dist, a = comm_fdiv_no_zero[1:5,], calc.FDiv = F, m = 4, calc.FRic = F, w.abun = T, stand.x = F, scale.RaoQ = F, calc.FGR = F, calc.CWM = F)
# # indices
# n <- nrow(comm_fdiv_no_zero)
# n <- 50 # For testing
#system.time(    
    
# res <- mclapply(c(1:n), function(i) {
#
#         message(paste("Computing FD indices for row i = ",i, sep = ""))
#         indices <- dbFD(x = euclid_dist, a = comm_fdiv_no_zero[i,], calc.FDiv = F, calc.FRic = F, w.abun = T, stand.x = F, scale.RaoQ = F, calc.FGR = F, calc.CWM = F)
#         # Return
#         return(data.frame(index = i, FEve = indices$FEve, FDis = indices$FDis, RaoQ = indices$RaoQ))
#
#     }, mc.cores = 20
#
# ) # eo lapply
#
# # ) # 50 rows --> 36 sec;
# # Rbind
# tab <- bind_rows(res)
# # summary(tab)


### 26/07/23: Run function based on 1000 rows to compute FEve/FDis/Rao'sQ and map to check if it makes sense
# indices <- dbFD(x = euclid_dist, a = comm_fdiv_no_zero, calc.FDiv = F, calc.FRic = F, w.abun = T, stand.x = F, scale.RaoQ = F, calc.FGR = F, calc.CWM = F)
# str(indices) # summary(indices)
# Map to check consistency with Jonas' previous results
# dat <- data.frame(x = na.omit(comm)[1:1000,"x"], y = na.omit(comm)[1:1000,"y"], FEve = indices$FEve, FDis = indices$FDis, RaoQ = indices$RaoQ)
# summary(dat)
# cor(dat$FDis, dat$RaoQ, method = "spearman") # yes, makes sense
# ggplot() + geom_tile(aes(x = x, y = y, fill = FDis), data = dat) + scale_fill_viridis(name = "FDis") +
#      #geom_contour(colour = "grey60", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = FEve), data = dat) +
#      geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
#      coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#          panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
#      scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#      scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

### 26/07/23: Testing with 'fundiversity' (need to convert all traits to numerics though)
# library("fundiversity")
# traits_fdiv2 <- traits_fdiv
# # str(traits_fdiv2)
# # Change class of vectors and convert spawning to 1/0 (1 = sac_spawning; 0 = free spawning)
# traits_fdiv2$sac_spawning <- NA
# traits_fdiv2[traits_fdiv2$Spawning.mode == "Sac","sac_spawning"] <- 1
# traits_fdiv2[traits_fdiv2$Spawning.mode == "Free","sac_spawning"] <- 0
# # Convert other traits to num
# traits_fdiv2$Myelination <- as.numeric(traits_fdiv2$Myelination)
# traits_fdiv2$Omnivore <- as.numeric(traits_fdiv2$Omnivore)
# traits_fdiv2$Carnivore <- as.numeric(traits_fdiv2$Carnivore)
# traits_fdiv2$Herbivore <- as.numeric(traits_fdiv2$Herbivore)
# traits_fdiv2$Detritivore <- as.numeric(traits_fdiv2$Detritivore)
# traits_fdiv2$Current <- as.numeric(traits_fdiv2$Current)
# traits_fdiv2$Cruise <- as.numeric(traits_fdiv2$Cruise)
# traits_fdiv2$Ambush <- as.numeric(traits_fdiv2$Ambush)
#
# # Test the fucntions from 'fundiversity'
# # ?fd_feve
# # na.omit(comm)[,"x"]
# FEVe <- fd_feve(traits = NULL, sp_com = comm_fdiv, dist_matrix = euclid_dist) # returns a data.frame
# # dim(FEVe) ; str(FEVe) ; summary(FEVe)
# FEVe$x <- na.omit(comm)[,"x"]
# FEVe$y <- na.omit(comm)[,"y"]
# Map FEve
# ggplot() + geom_raster(aes(x = x, y = y, fill = FEve), data = FEVe) + scale_fill_viridis(name = "FEve") +
#     geom_contour(colour = "grey60", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = FEve), data = FEVe) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#         panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
# ?fd_fdis
# FDis <- fd_fdis(traits = traits_fdiv2[,c(1,2,4:11)], sp_com = comm_fdiv)
# # Does not take dist matrix and does not take species with NA....
#
# # ?fd_raoq
# RaoQ <- fd_raoq(sp_com = comm_fdiv, dist_matrix = euclid_dist)
# # error

# fd_fric(traits = traits_fdiv2[,c(1,2,4:11)], sp_com = comm_fdiv, stand = FALSE)
### --> removes the species with NA in traits (althouhg they could be put in a PCoA with a Gower distance matrix.....)


### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------

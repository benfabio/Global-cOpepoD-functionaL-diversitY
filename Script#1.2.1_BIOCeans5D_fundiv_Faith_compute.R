
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 26/07/23: R script to compute Faith's index for FRic based on a Gower distance matrix and the community tables made during the MSc thesis of Jonas Wydler (data from Benedetti et al., 2023 - JBIO) © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - compute Faith's index of FRic based on a distance matrix and the global copepod community tables (using pd() from picante)
# - try different distance matrices (Gower, Euclid from Gower, Euclid from FAMD)
# - compute standardized effect sizes (SES) since Faith's index depends on species richness

### Latest update: 27/07/23

### ------------------------------------------------------------------------------------------------------------------------------------------------------

# install.packages("picante")
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
library("picante")

world <- map_data("world") # coastlines for maps

setwd("/net/kryo/work/fabioben/BIOCeans5D/data") # working dir

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 1) Load and format functional traits table 
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
rownames(traits_red2) <- traits_red2$Species
gow <- gawdis(x = traits_red2[,c(8,9,10,12:15,17:19)], groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))
# Derive dendrogram to be used by pd() function 
fit_gow <- flashClust::hclust(gow, method = "average") 
# plot(fit_gow) # plot(flashClust::hclust(gow, method = "ward"))


### 2) Load a standard community table and compute FD indices based on 'FD' functions for each grid cell
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/community_tables_05_01_21/contemp")
comm <- read.table("table_mon_composition_baseline_GAM_apr.txt")
# Subset traits_red2 to spp of interest
spp2keep <- colnames(comm)[c(4:length(comm))] #; spp2keep
commons <- intersect(spp2keep,traits_red2$Species)# ; length(commons) # FD indices based on 284 taxa 
comm_fdiv <- na.omit(comm[,c(commons)])


### 3) Compute Faith's index for 'comm_fdiv' based on phylo object derived from 'fit_gow'
# ?pd
# data(phylocom)
# pd(phylocom$sample, phylocom$phylo)
    
# Compute from HSI table
faith1 <- pd(samp = comm_fdiv, tree = as.phylo(fit_gow), include.root = T)
# summary(faith1)#; dim(faith1)

# Compute from 1/0 table (any HSI > 0 --> 1) to check consistency (it it automatically converts HSI to PA you should not use the HSI tables to compute Faith's index)
comm_fdiv_PA <- comm_fdiv
comm_fdiv_PA[comm_fdiv_PA > 0] <- 1 # dim(comm_fdiv_no_zero)
faith2 <- pd(samp = comm_fdiv_PA, tree = as.phylo(fit_gow), include.root = T)
#summary(faith2)
#cor(faith1$PD, faith2$PD) # if == 1, it means that pd() automatically converts HSI to 1/0 for any HSI > 0. 
### --> NEED TO USE ENSEMBLE OF THRESHOLDS (sensu Benedetti et al., 2021)



### 4) Compute SES PD
# ?ses.mpd
# ses.faith <- ses.pd(samp = comm_fdiv[1:100,], tree = as.phylo(fit_gow), null.model = "taxa.labels", runs = 99)
#summary(ses.faith)
#dim(ses.faith)
# OK good. Write RSCRIPTBATCH code for this

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
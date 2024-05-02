
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 28/07/23: R script to compute functional rarity indices (distinctiveness, uniqueness; Grenié et al., 2017) based on a Gower distance matrix and the community tables made during the MSc thesis of Jonas Wydler (data from Benedetti et al., 2023 - JBIO) © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - compute various functional rarity ('funrar') indices for global copepod communities

### Latest update: 28/07/23

### !!! WAS NOT KEPT FOR MANUSCRIPT IN PREP. AFTER ALL !!! YOU MAY DISREGARD THIS PART OF THE ANALYSES.

### ------------------------------------------------------------------------------------------------------------------------------------------------------

# install.packages("funrar")
library("marmap")
library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("viridis")
library("ggrepel")
library("gawdis")
library("funrar")
library("parallel")
library("xlsx")
library("readxl")
library("naniar")

world <- map_data("world") # coastlines for maps

setwd("/net/kryo/work/fabioben/BIOCeans5D/data") # working dir

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 1) Load and format functional traits table 
traits <- read.csv("traits_table_Benedetti2023.csv", h = T, sep = ";", dec = ",")
colnames(traits)[8] <- "Body.length" # size vector that we will keep 
traits <- traits %>% replace_with_na_all(condition = ~.x == "")
# Convert Feeding.mode, Trophic.group and Spawning.mode to factors
traits$Spawning.mode <- as.factor(traits$Spawning.mode)
traits$Trophic.group <- as.factor(traits$Trophic.group)
traits$Feeding.mode <- as.factor(traits$Feeding.mode)
# Count NA
names <- colnames(traits)[c(8:11,16)] 
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
rownames(traits_red2) <- traits_red2$Species

### Compute functional distance matrix
# ?compute_dist_matrix
dist.mat <- compute_dist_matrix(traits_table = traits_red2[,c(8,9,10,12:15,17:19)], metric = "gower")
# gow <- gawdis(x = traits_red2[,c(8,9,10,12:15,17:19)], groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))


### 2) Load a standard community table and compute FD indices based on 'FD' functions for each grid cell
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/community_tables_05_01_21/contemp")
comm <- read.table("table_mon_composition_baseline_GAM_jul.txt")
spp2keep <- colnames(comm)[c(4:length(comm))] #; spp2keep
commons <- intersect(spp2keep,traits_red2$Species)# ; length(commons) # FD indices based on 284 taxa 
comm_fdiv <- as.matrix(na.omit(comm[,c(commons)]))
# class(comm_fdiv) # needs to be a matrix for funrar()


### 3) Compute Faith's index for 'comm_fdiv' based on phylo object derived from 'fit_gow'
# ?funrar
#  From a site-species matrix and functional distance matrix compute all indices included in the package: functional uniqueness
# (regional, functional), functional distinctiveness (local, functional), geographical restrictedness (regional, extent),
# scarcity (local, abundance). *Note*: scarcity can only be computed if relative abundances are provided in the site-species matrix.
# Test on subset of 'comm_fdiv'
funrar.mat <- funrar(pres_matrix = comm_fdiv, dist_matrix = dist.mat, rel_abund = T)
#str(funrar.mat)
# class(funrar.mat) # returns a list

### What are the different indices computed?
# Ui : vector containing uniqueness values per species (Functional Uniqueness represents how "isolated" is a species in the global species
# pool, it is the functional distance to the nearest neighbor of the species of interest (see Details section for the formula))

# Di : site-species matrix with functional distinctiveness values, per species per site (The Functional Distinctiveness of a species is the average functional distance from a species to all the other in the given community)

# Ri : vector containing geographical restrictedness values per species (Geographical restrictedness is an index related to the extent of a species in a given dataset, it is close to 1 when the species is present in only a single site of the dataset (restricted) and close to 0 when the species is present at all sites. It estimates the geographical extent of a species in a dataset.)

# Si : site-species matrix with scarcity values per species per site. ONLY WHEN rel_abund = T. (Scarcity is close to 1 when a species is rare in a community and close to 0 when it is abundant)
 
### OK, so check uniqueness and restrictedness of the species 
# uniqueness
Uni <- funrar.mat[[1]]
# Uni[order(Uni$Ui, decreasing = T),]

Uni$species <- factor(Uni$species, levels = Uni[order(Uni$Ui),"species"]) # 
ggplot(Uni, aes(species, Ui)) + geom_point() + ylab("Uniqueness") +
    scale_x_discrete(name = "Species", labels = NULL) +
    theme_classic() + theme(axis.ticks = element_blank()) 
# Most species are not unique

# restrictedness
Rest <- funrar.mat[[3]]
# head(Rest)
#Rest[order(Rest$Ri, decreasing = T),]
Rest$species <- factor(Rest$species, levels = Rest[order(Rest$Ri),"species"]) # 
ggplot(Rest, aes(species, Ri)) + geom_point() + ylab("Restrictedness") +
    scale_x_discrete(name = "Species", labels = NULL) +
    theme_classic() + theme(axis.ticks = element_blank()) 
# Species seem to be more restricted than functionally unique (redundancy)
 
# Examine bivariate plot:
Rest$Ui <- Uni$Ui
rownames(Rest) <- Rest$species
ggplot(Rest, aes(Ri, Ui)) + geom_point() + xlab("Restrictedness") + ylab("Uniqueness") +
    geom_text_repel(aes(Ri, Ui, label = species), data = Rest[Rest$Ui > 0.025,]) + 
    theme_classic() + theme(axis.ticks = element_blank()) 


### Map mean functional distinctiveness and mean scarcity per sites
diss <- as.data.frame(funrar.mat[[2]])
# dim(diss); class(diss) # is a matrix
# str(diss)
# Compute mean distinctiveness with rowMeans
diss$distinctiveness <- base::rowMeans(x = as.matrix(diss), na.rm = T)
# summary(diss$distinctiveness)
diss$x <- na.omit(comm)[,"x"]
diss$y <- na.omit(comm)[,"y"]
# Map
ggplot() + geom_raster(aes(x = x, y = y, fill = distinctiveness), data = diss) + scale_fill_viridis(name = "Distinctiveness", option = "A") +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = distinctiveness), data = diss) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

### Scarcity
scar <- as.data.frame(funrar.mat[[4]])
# Compute mean distinctiveness with rowMeans
scar$scarcity <- base::rowMeans(x = as.matrix(scar), na.rm = T)
# summary(scar$scarcity)
scar$x <- na.omit(comm)[,"x"]
scar$y <- na.omit(comm)[,"y"]
# Map
ggplot() + geom_raster(aes(x = x, y = y, fill = scarcity), data = scar) + scale_fill_viridis(name = "Scarcity", option = "A") +
    geom_contour(colour = "grey30", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = scarcity), data = scar) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

### Works well. Run in R CMD BATCH

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
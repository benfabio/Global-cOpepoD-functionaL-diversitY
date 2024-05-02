
### ================================================================================================================

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
library("vegan")
library("biomod2")

setwd("/net/kryo/work/fabioben/GODLY/data") 

### ================================================================================================================

# Load funct traits data and reformat
traits <- read.csv("traits_table_Benedetti2023.csv", h = T, sep = ";", dec = ",")
colnames(traits)[8] <- "Body.length" # size vector that we will keep 
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

# To compute Gower's distance matrix with proper trait encoding
# gow <- gawdis(x = traits_red2[,c(8,9,10,12:15,17:19)], groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))

# The indices we are going to compute rely on presence-absence (1/0) data, not HSI (automatic conversion into 1/0 for any HSI > 0).
# Need to retrieve the TSS probability cutoffs from the scores tables
setwd("/net/kryo/work/fabioben/GODLY/data/sdm/evaluation_scores_05_01_21")
score.files <- dir()#; files
scores <- lapply(score.files, function(f) {d<-get(load(f)); return(d)})
tab.scores <- bind_rows(scores)
rm(score.files); gc()
tab.scores$SDM <- NA
tab.scores[grepl("GLM",rownames(tab.scores)),"SDM"] <- "GLM"
tab.scores[grepl("ANN",rownames(tab.scores)),"SDM"] <- "ANN"
tab.scores[grepl("GAM",rownames(tab.scores)),"SDM"] <- "GAM"
cutoffs <- data.frame( tab.scores %>% group_by(species,SDM) %>% summarise(cutoff = mean(Cutoff_TSS, na.rm = T)/1000) )
# cutoffs; to apply on a species x SDM level

# Load list of community tables
setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/contemp")
files <- dir()[grep("composition_baseline",dir())] #; files
# f <- files[36]
# f <- files[27]

mclapply(files, function(f) {
    
        message(paste("Calculating FD indices using dbFD() for: ", f, sep = ""))
        setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/contemp")
        
        comm <- read.table(f)
        rownames(comm) <- comm$cell_id
        
        # Check if there are columns with only NA?
        empties <- comm %>% keep(~all(is.na(.x))) %>% names
        
        if( length(empties) > 0 ) {
            comm <- comm %>% select(-all_of(empties))
        } # eo if loop
       
        # Subset traits table to spp of interest
        spp2keep <- colnames(comm)[c(4:length(comm))] #; spp2keep
        commons <- intersect(spp2keep,traits_red2$Species)
        comm_fdiv <- na.omit(comm[,c("cell_id","x","y",commons)])
        rm(comm); gc()

        # Subset both tables based on commons
        traits_fdiv <- traits_red2[traits_red2$Species %in% commons,c(8,9,10,12:15,17:19)]
        rownames(traits_fdiv) <- traits_red2[traits_red2$Species %in% commons,"Species"]
        
        # Compute the basic Gower distance matrix, with proper encoding of the traits selected 
        gow_fdiv <- gawdis(x = traits_fdiv, groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))
    
        # For computing Faith's index, need to use 1/0 community data (ensemble of thresholds)
        if( grepl("GLM",f) ) {
            sdm.cutoff <- cutoffs[cutoffs$SDM == "GLM",]
        } else if( grepl("GAM",f) ) {
            sdm.cutoff <- cutoffs[cutoffs$SDM == "GAM",]
        } else if( grepl("ANN",f) ) {
            sdm.cutoff <- cutoffs[cutoffs$SDM == "ANN",]
        } # eo if else loop 
            
        # Use those cutoffs to derive a PA (1/0) comm table
        comm_fdiv_PA <- comm_fdiv
        # For each column ('commons'), get the species' HSI threshold 't' and convert to 1/0 column by column
        # c <- "Paracalanus_parvus"
        for(c in commons) {
            t <- sdm.cutoff[sdm.cutoff$species == c,"cutoff"]
            comm_fdiv_PA[,c] <- bm_BinaryTransformation(data = comm_fdiv_PA[,c], threshold = t)
        } # eo for loop - c in commons
    
        # Check if there are species (i.e., columns) with only 0; retrive name of given species and remove from 'comm_fdiv_PA'
        cs <- colSums(comm_fdiv_PA[,c(4:length(comm_fdiv_PA))]) # for species with only 0
        rs <- rowSums(comm_fdiv_PA[,c(4:length(comm_fdiv_PA))]) # for rows/assemblages with only 0
        
        if( 0 %in% rs ) {
            
            cells2rm <- names(rs[rs == 0]) # finds the grid cell/assemblage to remove
            comm_fdiv_PA <- comm_fdiv_PA[!(row.names(comm_fdiv_PA) %in% cells2rm),]
            # dim(comm_fdiv_PA)
            
        } # eo if loop 0 in 'rs'
        
        if( 0 %in% cs ) {
            
            name2rm <- names(cs[cs == 0]) # finds the species/column to remove
            
            message(paste("!!! Removing ",name2rm," from community because only 0 !!!", sep = ""))
            comm_fdiv_PA <- comm_fdiv_PA %>% select(-all_of(name2rm))
            
            # then need to re-define 'commons' and re-compute 'gow_fdiv'
            commons <- intersect(colnames(comm_fdiv_PA[,c(4:length(comm_fdiv_PA))]), traits_red2$Species)
            
            traits_fdiv <- traits_red2[traits_red2$Species %in% commons,c(8,9,10,12:15,17:19)]
            rownames(traits_fdiv) <- traits_red2[traits_red2$Species %in% commons,"Species"]
            
            # re-compute the basic Gower distance matrix
            gow_fdiv <- gawdis(x = traits_fdiv, groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))
            pcoa <- wcmdscale(d = gow_fdiv, eig = T)
            pcoa.scores <- data.frame(pcoa$points)
            euclid.dist <- dist(pcoa.scores, method = "euclidean")
        
        } else {
            
            pcoa <- wcmdscale(d = gow_fdiv, eig = T)
            pcoa.scores <- data.frame(pcoa$points)
            euclid.dist <- dist(pcoa.scores, method = "euclidean")
            
        } # eo if else loop
    
        # Compute indices
        # indices <- dbFD(x = gow_fdiv, a = comm_fdiv_PA, calc.FDiv = T, calc.FRic = T, w.abun = F, m = 4, corr = "cailliez", scale.RaoQ = T, stand.x = F, calc.FGR = F, calc.CWM = F)
            
        indices <- dbFD(x = euclid.dist, a = comm_fdiv_PA[,c(4:length(comm_fdiv_PA))], calc.FDiv = T, calc.FRic = T, stand.FRic = T, w.abun = F, m = 4, corr = "cailliez", scale.RaoQ = T, stand.x = F, calc.FGR = F, calc.CWM = F)
        
        # Concatenate into a ddf 
        if( length(comm_fdiv_PA[,"x"]) == length(indices$FEve) ) {
            
            message(paste("Saving FD indices for: ", f,"\n", sep = ""))
         
            dat <- data.frame(x = comm_fdiv_PA[,"x"], y = comm_fdiv_PA[,"y"], SR = rowSums(comm_fdiv_PA[,c(4:length(comm_fdiv_PA))]),
                    FEve = indices$FEve, FDis = indices$FDis, RaoQ = indices$RaoQ, FDiv = indices$FDiv, FRic = indices$FRic)
            
            setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/db_FD/PA_Gawdis_PCoA_Euclid/Stand.Fric")
            
            filename <- str_replace(f,"txt","Rdata")
            filename <- str_replace(filename,"composition","FDindices")
            
            save(dat, file = filename)
            
        } else {
            
            message(paste("Failed to calculate FD indices for: ", f,"\n", sep = ""))
            
        } # eo if else loop     
    
        rm(indices,euclid_dist,pcoa2gow_fdiv,traits_fdiv,comm_fdiv)
        gc()
    
    }, mc.cores = 15
    
) # eo mclapply

rm(traits,traits_red,traits_red2)
gc()

### ================================================================================================================


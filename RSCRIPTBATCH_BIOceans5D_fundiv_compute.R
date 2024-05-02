
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

# Load list of community tables
setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/contemp")
files <- dir() #; files
# f <- files[13]

mclapply(files, function(f) {
    
        message(paste("Calculating FD indices using dbFD() for: ", f, sep = ""))
        
        setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/contemp")
        comm <- read.table(f)
        
        # Check if there are columns with only NA?
        empties <- comm %>% keep(~all(is.na(.x))) %>% names
        
        if( length(empties) > 0 ) {
            comm <- comm %>% select(-all_of(empties))
        } # eo if loop
       
        # Subset traits table to spp of interest
        spp2keep <- colnames(comm)[c(4:length(comm))] #; spp2keep
        commons <- intersect(spp2keep,traits_red2$Species)
        comm_fdiv <- na.omit(comm[,c(commons)])
        
        # Keep coordinates to concatenate with FD indices
        long <- na.omit(comm)[,"x"]
        lats <- na.omit(comm)[,"y"]
        rm(comm); gc()

        # Subset both tables based on commons
        traits_fdiv <- traits_red2[traits_red2$Species %in% commons,c(8,9,10,12:15,17:19)]
        rownames(traits_fdiv) <- traits_red2[traits_red2$Species %in% commons,"Species"]
        
        # Compute the basic Gower distance matrix, with proper encoding of the traits selected 
        gow_fdiv <- gawdis(x = traits_fdiv, groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))
    
        # Create the final euclidean distance matrix to use in dbFD()
        #pcoa2 <- ape::pcoa(gow_fdiv, correction = "cailliez")
        #euclid_dist <- dist(pcoa2$vectors[,1:4], method = "euclidean") 
        # For testing
        # comm_fdiv_no_zero <- comm_fdiv
        # comm_fdiv_no_zero[comm_fdiv_no_zero == 0] <- 0.000001 # dim(comm_fdiv_no_zero)
        # indices <- dbFD(x = gow_fdiv, a = comm_fdiv_no_zero[1:5,], calc.FRic = T, calc.FDiv = T, w.abun = T, m = 4, corr = "cailliez", stand.x = F, scale.RaoQ = F, calc.FGR = F, calc.CWM = F)
        ### NOTE: FDiv: Cannot be computed when 'calc.FRic' is FALSE
    
        # Compute indices  
        indices <- dbFD(x = gow_fdiv, a = comm_fdiv, calc.FRic = T, calc.FDiv = T, w.abun = T, m = 4, corr = "cailliez", stand.x = F, scale.RaoQ = F, calc.FGR = F, calc.CWM = F)
        
        # Concatenate into a ddf 
        if( length(long) == length(indices$FEve) ) {
            
            message(paste("Saving FD indices for: ", f,"\n", sep = ""))
            setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid")
         
            dat <- data.frame(x = long, y = lats, FEve = indices$FEve, FDis = indices$FDis, RaoQ = indices$RaoQ, FDiv = indices$FDiv, FRic = indices$FRic)
            filename <- str_replace(f,"txt","Rdata")
            filename <- str_replace(filename,"composition","FDindices")
            
            save(dat, file = filename)
            
        } else {
            
            message(paste("Failed to calculate FD indices for: ", f,"\n", sep = ""))
            
        } # eo if else loop     
    
        rm(indices,euclid_dist,pcoa2gow_fdiv,traits_fdiv,comm_fdiv)
        
        gc()
    
    }, mc.cores = 19
    
) # eo mclapply

rm(traits,traits_red,traits_red2)
gc()

### ================================================================================================================

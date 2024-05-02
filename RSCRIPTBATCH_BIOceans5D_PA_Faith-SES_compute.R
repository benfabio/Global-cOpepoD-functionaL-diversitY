
### ================================================================================================================

library("marmap")
library("tidyverse")
library("reshape2")
library("purrr")
library("RColorBrewer")
library("viridis")
library("gawdis")
library("parallel")
library("xlsx")
library("readxl")
library("flashClust")
library("naniar")
library("vegan")
library("picante")
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

# Faith's index works on PA (1/0) data, not HSI (automatic conversion into 1/0 for any HSI > 0).
# Therefore, need to retrieve the TSS probability cutoffs from the scores tables
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
files <- dir()
# f <- files[25] # for testing

mclapply(files, function(f) {
    
        message(paste("Calculating Faith indices of FRic using pd() for: ", f, sep = ""))
        
        setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/contemp")
        comm <- read.table(f)
        
        # Check if there are columns with only NA?
        empties <- comm %>% keep(~all(is.na(.x))) %>% names
        
        if( length(empties) > 0 ) {
            comm <- comm %>% select(-all_of(empties))
        } # eo if loop
        
        # Subset traits table to spp of interest
        spp2keep <- colnames(comm)[c(4:length(comm))]
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
        # Convert to hclust
        fit_gow <- flashClust::hclust(gow_fdiv, method = "average") 
    
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

        for(c in commons) {
            # message(c)
            t <- sdm.cutoff[sdm.cutoff$species == c,"cutoff"]
            comm_fdiv_PA[,c] <- bm_BinaryTransformation(data = comm_fdiv_PA[,c], threshold = t)
        } # eo for loop - c in commons
            
        # Reach out and touch
        faith <- pd(samp = comm_fdiv_PA, tree = as.phylo(fit_gow), include.root = T)

        # Compute SES as well
        ses <- ses.pd(samp = comm_fdiv_PA, tree = as.phylo(fit_gow), null.model = "taxa.labels", runs = 99)
        
        # Concatenate into a ddf 
        if( length(long) == nrow(faith) & length(long) == nrow(ses) ) {
            
                message(paste("Saving Faith indices and associated SES for: ", f,"\n", sep = ""))
         
                setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/Faith/Faith")
                dat1 <- data.frame(x = long, y = lats, Faith = faith$PD, SR = faith$SR)
                filename <- str_replace(f,"composition","Faith")
                filename <- str_replace(filename,"txt","Rdata")
                # Save
                save(dat1, file = filename)
            
                # ggplot() + geom_raster(aes(x = x, y = y, fill = Faith), data = dat1) + scale_fill_viridis(name = "Faith") +
#                      geom_contour(colour = "grey30", binwidth = 0.2, size = 0.25, aes(x = x, y = y, z = Faith), data = dat1) +
#                      geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#                      coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                          panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
#                      scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#                      scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
    
                # Same for SES
                setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/Faith/SES")
                dat2 <- data.frame(ses)
                dat2$x <- long
                dat2$y <- lats
                filename <- str_replace(f,"composition","SES")
                filename <- str_replace(filename,"txt","Rdata")
                # Save
                save(dat2, file = filename)
            
        } else {
            
                message(paste("Failed to calculate Faith index for: ", f,"\n", sep = ""))
            
        } # eo if else loop     
    
        rm(ses,faith,comm_fdiv_PA,filename,dat1,dat2,sdm.cutoff)
        gc()
            
    }, mc.cores = 25
    
) # eo mclapply

rm(traits,traits_red,traits_red2)
gc()


### ================================================================================================================

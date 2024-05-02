
### ================================================================================================================

library("marmap")
library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("viridis")
library("gawdis")
library("funrar")
library("parallel")
library("xlsx")
library("readxl")
library("flashClust")
library("naniar")
library("vegan")

world <- map_data("world") # coastlines for maps

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
rownames(traits_red2) <- traits_red2$Species

# Load list of community tables
setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/contemp")
files <- dir()[grep("composition",dir())] #; files
# f <- files[13]

lapply(files, function(f) {
    
        message(paste("\n","Calculating 'funrar' indices for: ", f, sep = ""))
        
        setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/contemp")
        comm <- read.table(f)
       
        # Check if there are columns with only NA?
        empties <- comm %>% keep(~all(is.na(.x))) %>% names
        
        if( length(empties) > 0 ) {
            comm <- comm %>% select(-all_of(empties))
        } # eo if loop
       
        spp2keep <- colnames(comm)[c(4:length(comm))] #; spp2keep
        commons <- intersect(spp2keep,traits_red2$Species)# ; length(commons) # FD indices based on 284 taxa 
        comm_fdiv <- as.matrix(na.omit(comm[,c(commons)]))
        
        # Keep coordinates to concatenate with FD indices
        long <- na.omit(comm)[,"x"]
        lats <- na.omit(comm)[,"y"]
        rm(comm); gc()

        # Subset both tables based on commons
        traits_fdiv <- traits_red2[traits_red2$Species %in% commons,c(8,9,10,12:15,17:19)]
        rownames(traits_fdiv) <- traits_red2[traits_red2$Species %in% commons,"Species"]
        
        # Compute the basic Gower distance matrix, with proper encoding of the traits selected 
        dist.mat <- compute_dist_matrix(traits_table = traits_red2[,c(8,9,10,12:15,17:19)], metric = "gower")
        
        # Compute indices
        indices <- funrar(pres_matrix = comm_fdiv, dist_matrix = dist.mat, rel_abund = T)
        # ?funrar
        # functional distinctiveness = the average functional distance from a species to all the other in the given community.
        # scarcity = a measure of the *local* rarity in terms of abundances/HSI. If S_i is close to 1 the species has a very
        # low abundances while if it's close to 0, it is quite abundant in the community ; 
        # is close to 1 when a species is rare in a community and close to 0 when it is abundant.
        
        # Concatenate into a ddf 
        diss <- as.data.frame(indices[[2]])
        scar <- as.data.frame(indices[[4]])
        # head(diss)
        
        if( length(long) == nrow(diss) & length(long) == nrow(scar) ) {
            
            message(paste("Saving funrar indices for: ", f,"\n", sep = ""))
         
            dat <- data.frame(x = long, y = lats, distinct = base::rowMeans(x = as.matrix(diss), na.rm = T), scarcity = base::rowMeans(x = as.matrix(scar), na.rm = T))
            
            setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/funrar/")
            filename <- str_replace(f,"txt","Rdata")
            filename <- str_replace(filename,"composition","funrar")
            
            save(dat, file = filename)
                     
            # ggplot() + geom_raster(aes(x = x, y = y, fill = scarcity), data = dat) + scale_fill_viridis(name = "D") +
 #                 geom_contour(colour = "grey30", binwidth = 0.2, size = 0.25, aes(x = x, y = y, z = scarcity), data = dat) +
 #                 geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
 #                 coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 #                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
 #                scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
 #                scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
            
        } else {
            
            message(paste("Failed to calculate funrar indices for: ", f,"\n", sep = ""))
            
        } # eo if else loop     
    
        rm(dist.mat,indices,traits_fdiv,comm_fdiv,dat,diss,scar)
        gc()
    
    }
    
) # eo mclapply


### Same but as above, but tis time retrieve Ui & Ri indices
# Ui : vector containing uniqueness values per species (Functional Uniqueness represents how "isolated" is a species in the global species
# pool, it is the functional distance to the nearest neighbor of the species of interest (see Details section for the formula))

# Ri : vector containing geographical restrictedness values per species (geographical restrictedness is an index related to the extent of a species in a given dataset, it is close to 1 when the species is present in only a single site of the dataset (restricted) and close to 0 when the species is present at all sites. It estimates the geographical extent of a species in a dataset.)


lapply(files, function(f) {
    
        message(paste("\n","Calculating funrar indices for: ", f, sep = ""))
        
        setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/contemp")
        comm <- read.table(f)
       
        # Check if there are columns with only NA?
        empties <- comm %>% keep(~all(is.na(.x))) %>% names
        
        if( length(empties) > 0 ) {
            comm <- comm %>% select(-all_of(empties))
        } # eo if loop
       
        spp2keep <- colnames(comm)[c(4:length(comm))] #; spp2keep
        commons <- intersect(spp2keep,traits_red2$Species)# ; length(commons) # FD indices based on 284 taxa 
        comm_fdiv <- as.matrix(na.omit(comm[,c(commons)]))
        
        # Keep coordinates to concatenate with FD indices
        long <- na.omit(comm)[,"x"]
        lats <- na.omit(comm)[,"y"]
        rm(comm); gc()

        # Subset both tables based on commons
        traits_fdiv <- traits_red2[traits_red2$Species %in% commons,c(8,9,10,12:15,17:19)]
        rownames(traits_fdiv) <- traits_red2[traits_red2$Species %in% commons,"Species"]
        
        # Compute the basic Gower distance matrix, with proper encoding of the traits selected 
        dist.mat <- compute_dist_matrix(traits_table = traits_red2[,c(8,9,10,12:15,17:19)], metric = "gower")
        
        # Compute indices
        indices <- funrar(pres_matrix = comm_fdiv, dist_matrix = dist.mat, rel_abund = T)
        # str(indices)
        
        # Retrieve indices
        uniqs <- as.data.frame(indices[[1]])
        restr <- as.data.frame(indices[[3]])
        
        if( nrow(uniqs) == nrow(restr) ) {
            
            message(paste("Saving 'funrar' indices for: ", f,"\n", sep = ""))
            dat <- uniqs
            dat$Ri <- restr$Ri
            # dat[order(dat$Ui, decreasing = T),]
            # Basic test plot
            # ggplot(aes(x = Ri, y = Ui), data = dat) + geom_point(size = 1) +
 #                #stat_density_2d(geom = "polygon", aes(alpha = (..level..)), fill = "#7fbc41") +
 #                #scale_alpha_continuous(range = c(0,1)) +
 #                xlab("Restrictedness") + ylab("Functional uniqueness") + theme_bw()   
            
            setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/funrar/")
            filename <- str_replace(f,"txt","Rdata")
            filename <- str_replace(filename,"composition","species_funrar_Ui+Ri")
            
            save(dat, file = filename)
            
        } else {
            
            message(paste("Failed to calculate funrar indices for: ", f,"\n", sep = ""))
            
        } # eo if else loop     
    
        rm(dist.mat,indices,traits_fdiv,comm_fdiv,dat,restr,uniqs)
        gc()
    
    }
    
) # eo mclapply

rm(traits,traits_red,traits_red2)
gc()

### ================================================================================================================


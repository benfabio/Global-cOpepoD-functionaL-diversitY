
### ================================================================================================================

# install.packages("phyloregion")
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
library("biomod2")
library("phyloregion")
library("ape")
library("picante")
library("Matrix")

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
setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/future")
files <- dir()[grep("composition",dir())] #; files

### In case function belows crashes because of memory usage issue son KRYO:
setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/beta.div")
files2 <- dir()[grep("2100-2000",dir())] #; files
files2 <- str_replace(files2,"Rdata","txt")
files2 <- str_replace(files2,"beta.fun.div","composition")
files2run <- setdiff(files,files2)
# f <- files2run[51]

mclapply(files2run, function(f) {
    
        message(paste("\n","Calculating beta fun div indices for: ", f, "\n", sep = ""))
        
        setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/future")
        
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
        # Convert to tree for BAT
        funct.tree <- hclust(gow_fdiv, method = "average")
        # Convert to object 'phylo' for phyloregion package's functions
        phylo.tree <- as.phylo(x = funct.tree)
        
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
        
        # Compute beta functional div
        # ?phylobeta_core
        # ?phylobeta
        # NOTE: x needs to be a 'sparse matrix', and phy needs to be a 'phylo' object
        beta.div <- phyloregion::phylobeta(x = as(as.matrix(comm_fdiv_PA[,c(4:length(comm_fdiv_PA))]),"sparseMatrix"), phy = phylo.tree, index.family = "jaccard")
        
        # Function aboves returns 3 distance matrices: dissimilarity (total+2 components) of each grid cell to another --> compute mean Jac/Jne/Jtu per grid cell
        beta.jac <- beta.div[[3]]
        beta.jtu <- beta.div[[1]]
        beta.jne <- beta.div[[2]]
        
        # Compute average dissim per grid cell
        mat.beta.jac <- as.matrix(beta.jac)
        mat.beta.jtu <- as.matrix(beta.jtu)
        mat.beta.jne <- as.matrix(beta.jne)
        
        # Replace 0 by NA
        mat.beta.jac[mat.beta.jac == 0.00000000] <- NA
        mat.beta.jtu[mat.beta.jtu == 0.00000000] <- NA
        mat.beta.jne[mat.beta.jne == 0.00000000] <- NA
        
        # Compute mean dissim with rowMeans
        mean.beta.jac <- rowMeans(mat.beta.jac, na.rm = T)
        mean.beta.jtu <- rowMeans(mat.beta.jtu, na.rm = T)
        mean.beta.jne <- rowMeans(mat.beta.jne, na.rm = T)
        # summary(mean.beta.jtu / mean.beta.jac) # beta ratio
        
        if( length(mean.beta.jac) == length(comm_fdiv_PA[,"x"]) ) {
            
            message(paste("\n", "Saving beta fun div for: ", f,"\n", sep = ""))
         
         
            dat <- data.frame(x = comm_fdiv_PA[,"x"], y = comm_fdiv_PA[,"y"], SR = rs,
                    beta.jac = mean.beta.jac, beta.jtu = mean.beta.jtu, beta.jne = mean.beta.jne)
            
            # Need to rotate x because future communities are on the 0-360° projection
            dat$x2 <- dat$x
            dat[which(dat$x > 180),"x2"] <- (dat[which(dat$x > 180),"x"]) - 360
            # Relocate
            dat <- dat %>% relocate(x2, .after = x)
         
            setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/beta.div/")
            
            filename <- str_replace(f,"txt","Rdata")
            filename <- str_replace(filename,"composition","beta.fun.div")
            
            save(dat, file = filename)
                     
            # ggplot() + geom_raster(aes(x = x2, y = y, fill = beta.jtu), data = dat) + scale_fill_viridis(name = "ß jtu") +
#                        geom_contour(colour = "grey30", binwidth = .1, size = 0.25, aes(x = x2, y = y, z = beta.jtu), data = dat) +
#                        geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#                        coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                           panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#                       scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#                       scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
    
        } else {
            
            message(paste("\n", "Failed to calculate beta fun div indices for: ", f,"\n", sep = ""))
            
        } # eo if else loop     
    
        rm(dat,traits_fdiv,comm_fdiv,dat,beta.div,beta.jac,beta.jtu,beta.jne,mat.beta.jac,mat.beta.jtu,mat.beta.jne,mean.beta.jac,mean.beta.jtu,mean.beta.jne)
            
        gc()
    
    }, mc.cores = 10
    
) # eo mclapply

rm(traits,traits_red,traits_red2)
gc()


### ================================================================================================================

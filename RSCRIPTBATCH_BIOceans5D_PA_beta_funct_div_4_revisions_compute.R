### ================================================================================================================

### 19/11/24: ADRESSING REVIEWER 2'S COMMENT FOR THE MINOR REVISIONS:
###           NEED TO ASSESS BIASES IN BETA-FD LINKED TO UNEVEN REPRESENTATION OF REGIONS
###           WITH HIGH SPECIES RICHNESS COMPARED TO REGIONS OF LOWER RICHNESS

### Re-run function above on GAM-based monthly communities and ONLY for a spatial subset of the whole ocean cell grid
### by making sure that all SR levels are equally represented

### ================================================================================================================

# install.packages("phyloregion")
library("marmap")
library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("viridis")
library("pals")
library("gawdis")
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

world <- map_data("world")

### ================================================================================================================

# Load funct traits data and reformat
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/")
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
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/evaluation_scores_05_01_21/")
score.files <- dir()#; files
scores <- lapply(score.files, function(f) {d<-get(load(f)); return(d)})
tab.scores <- bind_rows(scores)
rm(score.files); gc()
tab.scores$SDM <- NA
tab.scores[grepl("GLM",rownames(tab.scores)),"SDM"] <- "GLM"
tab.scores[grepl("ANN",rownames(tab.scores)),"SDM"] <- "ANN"
tab.scores[grepl("GAM",rownames(tab.scores)),"SDM"] <- "GAM"
cutoffs <- data.frame( tab.scores %>% group_by(species,SDM) %>% summarise(cutoff = mean(Cutoff_TSS, na.rm = T)/1000) )



# Load list of community tables
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/community_tables_05_01_21/")
files2 <- dir()[grep("composition",dir())]#; files2

# For adjusting FUN to spatial subset
# f <- files2[1]

mclapply(files2, function(f) {
    
        message(paste("\n","Calculating beta fun div indices for: ", f, "\n", sep = ""))
                
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
        
        ### 'rs' contains the SR estimates. Use this vector to subset 1000 grid cells based on their SR and make sure 
        ### high SR cells are not over-represented compared to low SR cells (10 bins of equal density).
        ### Re-compute beta-FD based on this subset.
        ### R code below made with some help from chatGPT3
        diversity_values <- rs

        # Number of bins (e.g., divide into 10 groups for even distribution)
        num_bins <- 10

        # Create bins
        bin_breaks <- seq(min(diversity_values), max(diversity_values), length.out = num_bins+1)
        bins <- cut(diversity_values, breaks = bin_breaks, include.lowest = TRUE)

        # Sample equally from each bin
        samples_per_bin <- 1000/num_bins
        if( samples_per_bin != round(samples_per_bin) ) {
              stop("1000 sites cannot be evenly divided across bins. Adjust the number of bins or sample size.")
        } # eo if loop

        sampled_indices <- unlist(lapply(levels(bins), function(bin) {
                    which_in_bin <- which(bins == bin)
                    sample(which_in_bin, min(samples_per_bin, length(which_in_bin)))
                } # eo FUN
            ) # eo lapply
        ) # eo unlist

        # Get the subset
        sampled_sites <- diversity_values[sampled_indices]

        # Check the distribution of sampled sites
        # hist(sampled_sites, breaks = 15, main = "Distribution of Sampled Diversity Values", xlab = "Diversity Value", col = "grey")
        # good
             
        cells2keep <- names(sampled_sites)  
        
        rm(sampled_sites,sampled_indices,samples_per_bin,bins,bin_breaks,num_bins,diversity_values)
        gc()   
        
        # Compute beta functional div
        # ?phylobeta_core
        # ?phylobeta
        # NOTE: x needs to be a 'sparse matrix', and phy needs to be a 'phylo' object
        beta.div <- phyloregion::phylobeta(x = as(as.matrix(comm_fdiv_PA[cells2keep,c(4:length(comm_fdiv_PA))]),"sparseMatrix"),
                    phy = phylo.tree, index.family = "jaccard")
        
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
        
        if( length(mean.beta.jac) == length(comm_fdiv_PA[cells2keep,"x"]) ) {
            
            message(paste("\n", "Saving beta fun div for: ", f,"\n", sep = ""))
         
            dat <- data.frame(x = comm_fdiv_PA[cells2keep,"x"], y = comm_fdiv_PA[cells2keep,"y"], SR = rs[cells2keep],
                    beta.jac = mean.beta.jac, beta.jtu = mean.beta.jtu, beta.jne = mean.beta.jne)
                        
            filename <- str_replace(f,"txt","Rdata")
            filename <- str_replace(filename,"composition","beta.fun.div")
            filename <- str_replace(filename,"table","table_minor_revisions")
            
            save(dat, file = filename)
            #write.table(x = dat, file = filename, sep = "\t")
    
        } else {
            
            message(paste("\n", "Failed to calculate beta fun div indices for: ", f,"\n", sep = ""))
            
        } # eo if else loop     
    
        rm(dat,traits_fdiv,comm_fdiv,dat,beta.div,beta.jac,beta.jtu,beta.jne,mat.beta.jac,mat.beta.jtu,mat.beta.jne,mean.beta.jac,mean.beta.jtu,mean.beta.jne)
            
        gc()
    
    }, mc.cores = 2
    
) # eo mclapply


rm(traits,traits_red,traits_red2)
gc()

### ================================================================================================================

### Load the beta-FD indices obtained based on the code above. Binf and compute average across months to get 
### mean annual estimates of beta-FD indices. Compare mean trait dissimilarity to SR based on this spayial subsets.

setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/community_tables_05_01_21/")
files2 <- dir()[grep("table_minor_revisions",dir())]
# f <- files2[1]
res <- lapply(files2, function(f) {
            d <- get(load(f))
            filename <- str_replace_all(f,".Rdata","")
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[7]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[8]
            return(d)
        }
) # eo lapply
tab <- bind_rows(res)
rm(res); gc()
#dim(tab); head(tab); summary(tab)
# Add cell id
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) # length(unique(tab$cell_id))

### Examine distributions
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
tab$month <- factor(tab$month,months)
tab$ratio <- tab$beta.jtu/tab$beta.jac

### Compute mean annual indices
ann <- data.frame(
    tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
        rich = mean(SR, na.rm = T),
        jac = mean(beta.jac, na.rm = T),
        jtu = mean(beta.jtu, na.rm = T),
        jne = mean(beta.jne, na.rm = T),
        beta.ratio = mean(ratio, na.rm = T)
    ) 
) # eo ddf
#dim(ann)
#summary(ann)

# jac <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac), data = ann) +
#     scale_fill_gradientn(name = "Trait dissimilarity\n(Jaccard index)", colours = parula(100), guide = "colourbar") +
#     geom_contour(colour = "black", binwidth = .05, size = .4, aes(x = x, y = y, z = jac), data = ann) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
# jtu <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu), data = ann) +
#     scale_fill_gradientn(name = "Trait turnover\n(Jtu)", colours = parula(100), guide = "colourbar") +
#     geom_contour(colour = "black", binwidth = .05, size = .4, aes(x = x, y = y, z = jtu), data = ann) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
# jne <- ggplot() + geom_raster(aes(x = x, y = y, fill = jne), data = ann) +
#     scale_fill_gradientn(name = "Trait nestedness\n(Jne)", colours = parula(100), guide = "colourbar") +
#     geom_contour(colour = "black", binwidth = .02, size = .4, aes(x = x, y = y, z = jne), data = ann) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
# ratio <- ggplot() + geom_raster(aes(x = x, y = y, fill = beta.ratio), data = ann) +
#     scale_fill_gradientn(name = "Ratio (Jtu/Jac)", colours = parula(100), guide = "colourbar") +
#     geom_contour(colour = "black", binwidth = .05, size = .4, aes(x = x, y = y, z = beta.ratio), data = ann) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
# # Save
# require("ggpubr")
# fig2 <- ggarrange(jac,jtu,jne,ratio, align = 'hv', ncol = 2, nrow = 2, labels = letters[1:4])

### Looks similar as before.

### Plot beta FD ~ SR

library("ggpmisc")
formula1 <- y ~ x
formula2 <- y ~ poly(x,2)
formula3 <- y ~ poly(x,3)

p1 <- ggplot(ann, aes(x = rich, y = jac, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula2) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula2, label.y = "top", label.x = "right", size = 3) + 
  xlab("Species richness") + ylab("Trait dissimilarity (Jaccard)") + theme_bw()

p2 <- ggplot(ann, aes(x = rich, y = jtu, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula2) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula2, label.y = "top", label.x = "right", size = 3) + 
  xlab("Species richness") + ylab("Trait turnover") + theme_bw()

p3 <- ggplot(ann, aes(x = rich, y = jne, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula2) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula2, label.y = "top", label.x = "right", size = 3) + 
  xlab("Species richness") + ylab("Trait nestedness") + theme_bw()

panel <- ggarrange(p1,p2,p3, align = 'hv', ncol = 3, nrow = 1, labels = letters, common.legend = T)
ggsave(plot = panel, filename = "Fig.X_betaxSR_rarefied_revised_19.11.24.jpg", dpi = 300, width = 12, height = 4.5)


### ================================================================================================================
### ================================================================================================================
### ================================================================================================================


##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 15/08/23: R script to extract and compute the ∆%values of the chosen FD indices based on the monthly contemporary and future communities © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - load and compute the monthly (NOT annual for once) ∆%values of the FD indices chosen 
# - compute monthly ∆FD for each indices to evaluate impact of CC on zooplankton functional diversity
# - integrate all monthly ∆FD in SDM-specific and ESM-specific annual means; compute ensemble annual mean ∆FD
# - save in common data.frame

### Latest update: 18/09/23

### ------------------------------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("ncdf4")
library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("viridis")
library("parallel")
library("ggpubr")
library("ggthemes")
library("marmap")

world <- map_data("world") # coastlines for maps
world2 <- map_data("world2") # coastlines for maps

setwd("/net/kryo/work/fabioben/BIOCeans5D/data") # working dir

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 1) Choose and retrieve the mean annual values of the chosen indices of copepod FD
# - SR (taxonomic richness) based on PA 
# - Faith's index (functional richness based on length of branches of functional dist matrix) based on PA data
# - SES of Faith
# - Functional beta.div based on PA data 
# - FEve (functional eveness) based on HSI data
# - FDis (functional dispersion) based on HSI data
# - Rao'Q (dispersion too, through quadratic entropy) based on HSI data
# - FDiv (divergence) based on based on HSI data

# Define vector of FD indices for which you have to compute ∆ between monthly future and monthly contemporary values
indices <- c("SR","FR","SES.FR","Jac","Jne","Jtu","FEve","FDis","RaoQ","FDiv")
# i <- "SES.FR"

for(i in indices) {
    
    message(paste("\n","Calculating ∆ values for ",i,"\n", sep = ""))
    
    if(i %in% c("SR","FR")) {
        setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/Faith/Faith")
    } else if (i == "SES.FR") {
        setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/Faith/SES")
    } else if (i %in% c("Jac","Jne","Jtu")) {
        setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/beta.div")
    } else {
        setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid")    
    } # eo 1st if else loop
    
    # For each 180 "future" fields (12*3*5)
    files <- dir()[grep("2100-2000",dir())] #; files # should be 180 files because: 12 months x 3SDMs x 5 ESMs
    baseline.files <- dir()[!grepl("2100-2000", dir())]
    baseline.files <- baseline.files[!grepl("_diff_", baseline.files)]
    
    # f <- files[36] # for testing
    # f <- "table_mon_SES_2100-2000_CESM-BEC_GLM_sep.Rdata"
    mclapply(files, function(f) {
        
        #for(f in files) {
        
            message(paste(f, sep = ""))
            fut <- get(load(f))
            filename <- str_replace_all(f,".Rdata","")
            ESM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
            month <- unlist(strsplit(x = filename, split = "_", fixed = T))[7]
            # ESM; SDM; month
            
            baseline.file <- baseline.files[grepl(month, baseline.files)]
            baseline.file <- baseline.file[grepl(SDM, baseline.file)]
            base <- get(load(baseline.file))
            # dim(base); head(base)
            
            if(i == "SR") {
                
                base <- base[,c("x","y","SR")]
                base$cell_id <- factor(paste(base$x, base$y, sep = "_"))
                base <- base[order(base$cell_id),]
                base <- base %>% relocate(cell_id, .before = x)
                # Same with 'fut'
                fut <- fut[,c("x","y","SR")]
                fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
                fut <- fut[order(fut$cell_id),]
                fut <- fut %>% relocate(cell_id, .before = x)
                
            } else if(i == "FR") {
                
                base <- base[,c("x","y","Faith")]
                base$cell_id <- factor(paste(base$x, base$y, sep = "_"))
                base <- base[order(base$cell_id),]
                base <- base %>% relocate(cell_id, .before = x)
                # Same with 'fut'
                fut <- fut[,c("x","y","Faith")]
                fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
                fut <- fut[order(fut$cell_id),]
                fut <- fut %>% relocate(cell_id, .before = x)  
            
            } else if(i == "SES.FR") {
                
                base <- base[,c("x","y","pd.obs.z")]
                base$cell_id <- factor(paste(base$x, base$y, sep = "_"))
                base <- base[order(base$cell_id),]
                base <- base %>% relocate(cell_id, .before = x)
                # Same with 'fut'
                fut <- fut[,c("x","y","pd.obs.z")]
                fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
                fut <- fut[order(fut$cell_id),]
                fut <- fut %>% relocate(cell_id, .before = x)    
                
            } else if(i == "Jac") {
                
                base <- base[,c("x","y","beta.jac")]
                base$cell_id <- factor(paste(base$x, base$y, sep = "_"))
                base <- base[order(base$cell_id),]
                base <- base %>% relocate(cell_id, .before = x)
                # Same with 'fut'
                fut <- fut[,c("x","y","beta.jac")]
                fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
                fut <- fut[order(fut$cell_id),]
                fut <- fut %>% relocate(cell_id, .before = x)
                
            } else if(i == "Jtu") {
                
                base <- base[,c("x","y","beta.jtu")]
                base$cell_id <- factor(paste(base$x, base$y, sep = "_"))
                base <- base[order(base$cell_id),]
                base <- base %>% relocate(cell_id, .before = x)
                # Same with 'fut'
                fut <- fut[,c("x","y","beta.jtu")]
                fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
                fut <- fut[order(fut$cell_id),]
                fut <- fut %>% relocate(cell_id, .before = x)
                
            } else if(i == "Jne") {
                
                base <- base[,c("x","y","beta.jne")]
                base$cell_id <- factor(paste(base$x, base$y, sep = "_"))
                base <- base[order(base$cell_id),]
                base <- base %>% relocate(cell_id, .before = x)
                # Same with 'fut'
                fut <- fut[,c("x","y","beta.jne")]
                fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
                fut <- fut[order(fut$cell_id),]
                fut <- fut %>% relocate(cell_id, .before = x)
                
            } else if(i == "FEve") {
                
                base <- base[,c("x","y","FEve")]
                base$cell_id <- factor(paste(base$x, base$y, sep = "_"))
                base <- base[order(base$cell_id),]
                base <- base %>% relocate(cell_id, .before = x)
                # Same with 'fut'
                fut <- fut[,c("x","y","FEve")]
                fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
                fut <- fut[order(fut$cell_id),]
                fut <- fut %>% relocate(cell_id, .before = x)
                
            } else if(i == "FDis") {
                
                base <- base[,c("x","y","FDis")]
                base$cell_id <- factor(paste(base$x, base$y, sep = "_"))
                base <- base[order(base$cell_id),]
                base <- base %>% relocate(cell_id, .before = x)
                # Same with 'fut'
                fut <- fut[,c("x","y","FDis")]
                fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
                fut <- fut[order(fut$cell_id),]
                fut <- fut %>% relocate(cell_id, .before = x)
                
            } else if(i == "RaoQ") {
                
                base <- base[,c("x","y","RaoQ")]
                base$cell_id <- factor(paste(base$x, base$y, sep = "_"))
                base <- base[order(base$cell_id),]
                base <- base %>% relocate(cell_id, .before = x)
                # Re-scale Rao's Q
                base$RaoQ <- base$RaoQ/max(base$RaoQ, na.rm = T)
                # Same with 'fut'
                fut <- fut[,c("x","y","RaoQ")]
                fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
                fut <- fut[order(fut$cell_id),]
                fut <- fut %>% relocate(cell_id, .before = x)
                fut$RaoQ <- fut$RaoQ/max(fut$RaoQ, na.rm = T)
                
            } else if(i == "FDiv") {
                
                base <- base[,c("x","y","FDiv")]
                base$cell_id <- factor(paste(base$x, base$y, sep = "_"))
                base <- base[order(base$cell_id),]
                base <- base %>% relocate(cell_id, .before = x)
                # Same with 'fut'
                fut <- fut[,c("x","y","FDiv")]
                fut$cell_id <- factor(paste(fut$x, fut$y, sep = "_"))
                fut <- fut[order(fut$cell_id),]
                fut <- fut %>% relocate(cell_id, .before = x)
                
            } # eo 2nd if else loop
            
            # Sanity checks
            if( exists("base") & exists("fut") ) {
                
                # OK, files defined. But 'fut' are defined on the 0-360° longitudes
                fut$x2 <- fut$x
                fut[fut$x > 179.5,"x2"] <- (fut[fut$x > 179.5,"x"]) - 360
                # summary(base$x) ; summary(fut$x2)
                fut$cell_id <- factor(paste(fut$x2, fut$y, sep = "_"))
                fut <- fut[order(fut$cell_id),]
                
                # Adjust colnames
                colnames(base)[4] <- "base.val"
                colnames(fut)[4] <- "fut.val"
                
                # Compute diff for commons
                commons <- intersect(unique(fut$cell_id), unique(base$cell_id)) # length(commons)
                # Restrict to common cells
                fut2 <- fut[fut$cell_id %in% commons,]
                base2 <- base[base$cell_id %in% commons,]
                
                if( nrow(fut2) == nrow(base2) ) {
                    
                    # Compute diffs and %
                    
                    ### WARNING! Sometimes, 'base' value == 0 --> leads to -Inf in 'perc'
                    # nrow(tab[tab$base.val == 0,])/nrow(tab)
                    ### But usually correspo ds to VERY few instances. Replace 0 by smalles positive values 
                    
                    if( 0 %in% unique(base2$base.val) ) {
                        
                        base2$fut.val <- fut2$fut.val
                        
                        smallest <- min(base2[base2$base.val > 0,"base.val"], na.rm = T)
                        base2[base2$base.val == 0 & !is.na(base2$base.val),"base.val"] <- smallest
                        
                        base2$diff <- (base2$fut.val) - (base2$base.val)
                        base2$perc <- base2$diff/base2$base.val
                        # summary(base2)
                        
                    } else {
                        
                        base2$fut.val <- fut2$fut.val
                        base2$diff <- (base2$fut.val) - (base2$base.val)
                        base2$perc <- base2$diff/base2$base.val
                        # summary(base2)
    
                    } # eo if else loop
                    
                } else {
                    
                    message("!!! ISSUE !!! BASE AND FUT OBJECTS DO NOT HAVE SAME DIMENSIONS !!!")
                    
                } # eo 2nd sanity check
                
                # Save 'base2' as .RData in new dir
                save(base2, file = paste("table_mon_diff_",i,"_",ESM,"_",SDM,"_",month,".Rdata", sep = ""))
                                
            } else {
                
                message("!!! ISSUE !!! BASE AND FUT OBJECTS NOT DEFINED !!!")
                
            } # eo 1st sanity check
        
        }, mc.cores = 25
        
    ) # eo mclapply - f in FILES
    
} # eo FOR LOOP - i in INDICES


### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 2) Calculate mean ensemble values for each FD indices and assess variability per SDm and ESM
indices <- c("SR","FR","SES.FR","Jac","Jne","Jtu","FEve","FDis","RaoQ","FDiv")
#i <- "FEve" # for testing

for(i in indices) {
    
    message(paste("\n","Calculating ∆ values for ",i,"\n", sep = ""))
    
    if(i %in% c("SR","FR")) {
        setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/Faith/Faith")
    } else if (i == "SES.FR") {
        setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/Faith/SES")
    } else if (i %in% c("Jac","Jne","Jtu")) {
        setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/beta.div")
    } else {
        setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid")    
    } # eo 1st if else loop

    # Vector files of interest
    files <- dir()[grep(i,dir())] # f <- files[1]
    
    res <- mclapply(files, function(f) {
                d <- get(load(f))
                filename <- str_replace_all(f,".Rdata","")
                # unlist(strsplit(x = filename, split = "_", fixed = T))
                d$ESM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
                d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
                d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[7]
                return(d)
            }, mc.cores = 20
    ) # eo mclapply
    # Rbind
    tab <- bind_rows(res)
    rownames(tab) <- NULL
    rm(res); gc()
    # dim(tab); head(tab)
    # summary(tab)
        
    # Compute ensemble mean for 'perc'
    ens <- data.frame(tab %>% group_by(cell_id) %>% 
            summarize(x = unique(x), y = unique(y),
                mean = mean(perc, na.rm = T),
                std = sd(perc, na.rm = T)
        )
    ) # eo ddf
    # dim(ens); summary(ens)
    
    # Rotate to 0-360° if needed (Niki usually prefers it)
    ens$x2 <- ens$x
    ens[ens$x < 0,"x2"] <- (ens[ens$x > 179.5,"x"])+360
    # summary(ens$x2)

    ggplot() + geom_raster(aes(x = x, y = y, fill = mean), data = ens) +
         scale_fill_gradient2(name = paste("Mean ∆",i, sep = ""), low = "#3288bd", mid = "white", high = "#b2182b") +
         geom_contour(colour = "grey30", binwidth = .2, size = 0.25, aes(x = x, y = y, z = mean), data = ens) +
         geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
         coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
         scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
         scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
    
    # Map for stdev
    ggplot() + geom_raster(aes(x = x, y = y, fill = std), data = ens[ens$std < .5,]) +
         geom_raster(aes(x = x, y = y), data = ens[ens$std >= .5,], fill = "#800026") +
         scale_fill_distiller(name = paste("Stdev ∆",i, sep = ""), palette = "YlOrRd", direction = 1) +
         #geom_contour(colour = "grey30", binwidth = .1, size = 0.25, aes(x = x, y = y, z = std), data = ens) +
         geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
         coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
         scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
         scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)


} # eo FOR LOOP - i in INDICES

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
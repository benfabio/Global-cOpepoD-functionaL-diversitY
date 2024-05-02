
### ================================================================================================================

library("tidyverse")
library("reshape2")
library("parallel")

setwd("/net/kryo/work/fabioben/GODLY/data") 

### ================================================================================================================

### Script to compute the species-level spatial shifts between contemporary and future distributions (based on HSI).
### Compute per SDM and ESM, on a mean annual scale.

# Define vectors
SDMs <- c("GAM","GLM","ANN")
ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO") 

### For each SDM and ESM combin., compute mean annual HSI at the species-level and compute shift in centroid (weighted by HSI)
# For testing: 
# sdm <- "GAM"
# esm <- "IPSL-PISCES"

for(sdm in SDMs) {
    
    message(paste("\n","Computing position of the mean annual centroids for ",sdm,"\n", sep = ""))
    
    setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/contemp")
    
    files <- dir()[grep(sdm,dir())]
    # Compute species' mean annual HSI
    res <- mclapply(files, function(f) {
            d <- read.table(f, h = T)
            m.d <- melt(d, id.vars = c('cell_id','x','y'))
            rm(d); gc()
            colnames(m.d)[c(4,5)] <- c("species","HSI")
            return(m.d)
        }, mc.cores = length(months)
    ) # eo mclapply
    # Rbind
    tab <- dplyr::bind_rows(res)
    rm(res); gc()
    
    # Compute mean annual HSI for each species 
    ann <- data.frame( tab %>% group_by(species,cell_id) %>% summarize(x = unique(x), y = unique(y), HSI = mean(HSI, na.rm = T)) )
    # str(ann); unique(ann$species); head(ann); summary(ann$HSI)
    ann2 <- na.omit(ann)
    
    # Compute various possible centroids 
    centroids.base <- data.frame(ann2 %>% group_by(species) %>% summarize(lon = weighted.mean(x = x, w = HSI), lat = weighted.mean(x = abs(y), w = HSI)))
    # N Hemisphere only
    centroids.base.NH <- data.frame(ann2[ann2$y >= 0,] %>% group_by(species) %>% summarize(lon = weighted.mean(x = x, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
    # S Hemisphere only
    centroids.base.SH <- data.frame(ann2[ann2$y <= 0,] %>% group_by(species) %>% summarize(lon = weighted.mean(x = x, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
    
    # Rbind and save in dir for later analyses
    centroids.base$region <- "Global"
    centroids.base.NH$region <- "NH"
    centroids.base.SH$region <- "SH"
    centroids <- rbind(centroids.base, centroids.base.NH, centroids.base.SH)
    centroids$SDM <- sdm
    setwd("/net/kryo/work/fabioben/GODLY/data/shift_metrics")
    save(x = centroids, file = paste("table_contemp_mean_ann_centroids_species_",sdm,"_20_11_23.RData", sep = ""))
    
    rm(ann,ann2,tab,centroids.base.SH,centroids.base.NH,centroids.base); gc()
    
    ### Compute the centroids based on the future mean annual HSI
    for(esm in ESMs) {
        
         message(paste("\n","Computing position of the mean annual centroids for ",sdm,"x",esm,"\n", sep = ""))
        
         setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/future")
         files <- dir()[grep(sdm,dir())]  
         files <- files[grepl(esm,files)]     
         
         res <- mclapply(files, function(f) {
                 d <- read.table(f, h = T)
                 # Move longitudes back to -180°/180° # unique(m.base.zoo$x)
                 d$x2 <- d$x
                 d[which(d$x > 180),c("x2")] <- d[which(d$x > 180),"x"] - 360
                 # Drop 'x'
                 d <- subset(d, select = -c(x))
                 m.d <- melt(d, id.vars = c('cell_id','x2','y'))
                 rm(d); gc()
                 colnames(m.d)[c(4,5)] <- c("species","HSI")
                 return(m.d)
             }, mc.cores = length(months)
         ) # eo mclapply
         # Rbind
         tab <- dplyr::bind_rows(res)
         rm(res); gc()  
         
         # Compute mean annual HSI for each species 
         ann <- data.frame( tab %>% group_by(species,cell_id) %>% summarize(x = unique(x2), y = unique(y), HSI = mean(HSI, na.rm = T)) )
         # str(ann); unique(ann$species); head(ann); summary(ann$HSI)
         ann2 <- na.omit(ann)
    
         # Compute various possible centroids 
         centroids.fut <- data.frame(ann2 %>% group_by(species) %>% summarize(lon = weighted.mean(x = x, w = HSI), lat = weighted.mean(x = abs(y), w = HSI)))
         # N Hemisphere only
         centroids.fut.NH <- data.frame(ann2[ann2$y >= 0,] %>% group_by(species) %>% summarize(lon = weighted.mean(x = x, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
         # S Hemisphere only
         centroids.fut.SH <- data.frame(ann2[ann2$y <= 0,] %>% group_by(species) %>% summarize(lon = weighted.mean(x = x, w = HSI), lat = weighted.mean(x = y, w = HSI) ) )
         
         # Combine and save 
         centroids.fut$region <- "Global"
         centroids.fut.NH$region <- "NH"
         centroids.fut.SH$region <- "SH"
         centroids2 <- rbind(centroids.fut, centroids.fut.NH, centroids.fut.SH)
         centroids2$SDM <- sdm
         centroids2$ESM <- esm
         
         setwd("/net/kryo/work/fabioben/GODLY/data/shift_metrics")
         save(x = centroids2, file = paste("table_fut_mean_ann_centroids_species_",sdm,"_",esm,"_20_11_23.RData", sep = ""))
         
         rm(ann2,tab,centroids.fut,centroids.fut.NH,centroids.fut.SH); gc()
        
    } # eo for loop - esm in ESMs
    
} # eo for loop - sdm in SDMs


### ================================================================================================================
### ================================================================================================================
### ================================================================================================================

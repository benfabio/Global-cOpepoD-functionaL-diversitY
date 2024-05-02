
### ================================================================================================================

library("tidyverse")
library("reshape2")
library("parallel")

setwd("/net/kryo/work/fabioben/GODLY/data") 

### ================================================================================================================

### Script to compute the species-level changes in mean HSI between contemporary and future distributions.
### Compute per SDM and ESM, on a mean annual scale.

# Define vectors
SDMs <- c("GAM","GLM","ANN")
ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO") 

### For each SDM and ESM combin., compute mean annual HSI at the species-level and compute difference sin mean HSI
# For testing: 
# sdm <- "GAM"
# esm <- "GFDL-TOPAZ"

for(sdm in SDMs) {
    
    message(paste("\n","Computing mean annual HSI for ",sdm,"\n", sep = ""))
    
    setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/contemp")
    
    files <- dir()[grep(sdm,dir())]
    # Compute species' mean annual HSI
    res <- mclapply(files, function(f) {
            d <- read.table(f, h = T)
            m.d <- melt(d, id.vars = c('cell_id','x','y'))
            rm(d); gc()
            colnames(m.d)[c(4,5)] <- c("species","HSI")
            return(m.d)
        }, mc.cores = length(files)
    ) # eo mclapply
    # Rbind
    tab <- dplyr::bind_rows(res)
    rm(res); gc()
    
    # Compute mean annual HSI for each species 
    ann <- data.frame( tab %>% group_by(species,cell_id) %>% summarize(x = unique(x), y = unique(y), HSI = mean(HSI, na.rm = T)))
    # str(ann); unique(ann$species); head(ann); summary(ann$HSI)
    ann2 <- na.omit(ann)
    
    # Compute 3 possible spatial averages of HSI: globa, NH and SH
    avg.hsi.global <- data.frame(ann2 %>% group_by(species) %>% summarize(mean = mean(HSI, na.rm = T), sd = sd(HSI, na.rm = T)))
    # N Hemisphere only
    avg.hsi.NH <- data.frame(ann2[ann2$y >= 0,] %>% group_by(species) %>% summarize(mean = mean(HSI, na.rm = T), sd = sd(HSI, na.rm = T)))
    # S Hemisphere only
    avg.hsi.SH <- data.frame(ann2[ann2$y <= 0,] %>% group_by(species) %>% summarize(mean = mean(HSI, na.rm = T), sd = sd(HSI, na.rm = T)))
    
    # Rbind and save in dir for later analyses
    avg.hsi.global$region <- "Global"
    avg.hsi.NH$region <- "NH"
    avg.hsi.SH$region <- "SH"
    avg.hsi.all <- rbind(avg.hsi.global, avg.hsi.NH, avg.hsi.SH)
    avg.hsi.all$SDM <- sdm
    
    setwd("/net/kryo/work/fabioben/GODLY/data/shift_metrics")
    save(x = avg.hsi.all, file = paste("table_contemp_mean_ann_hsi_species_",sdm,"_30_11_23.RData", sep = ""))
    
    rm(ann,ann2,tab,avg.hsi.global,avg.hsi.NH,avg.hsi.SH,avg.hsi.all); gc()
    
    
    ### Compute the average HSI based on the future community table
    for(esm in ESMs) {
        
         message(paste("\n","Computing mean annual HSI for ",sdm,"x",esm,"\n", sep = ""))
        
         setwd("/net/kryo/work/fabioben/GODLY/data/community_tables_05_01_21/future")
         files <- dir()[grep(sdm,dir())]  
         files <- files[grepl(esm,files)]     
         
         res <- mclapply(files, function(f) {
                 d <- read.table(f, h = T)
                 d$x2 <- d$x
                 d[which(d$x > 180),c("x2")] <- d[which(d$x > 180),"x"] - 360
                 d <- subset(d, select = -c(x))
                 m.d <- melt(d, id.vars = c('cell_id','x2','y'))
                 rm(d); gc()
                 colnames(m.d)[c(4,5)] <- c("species","HSI")
                 return(m.d)
             }, mc.cores = length(files)
         ) # eo mclapply
         tab <- dplyr::bind_rows(res)
         rm(res); gc()  
         
         # Compute mean annual HSI for each species 
         ann <- data.frame( tab %>% group_by(species,cell_id) %>% summarize(x = unique(x2), y = unique(y), HSI = mean(HSI, na.rm = T)) )
         ann2 <- na.omit(ann)
    
         # Same as for baseline conditions
         avg.hsi.global <- data.frame(ann2 %>% group_by(species) %>% summarize(mean = mean(HSI, na.rm = T), sd = sd(HSI, na.rm = T)))
         avg.hsi.NH <- data.frame(ann2[ann2$y >= 0,] %>% group_by(species) %>% summarize(mean = mean(HSI, na.rm = T), sd = sd(HSI, na.rm = T)))
         avg.hsi.SH <- data.frame(ann2[ann2$y <= 0,] %>% group_by(species) %>% summarize(mean = mean(HSI, na.rm = T), sd = sd(HSI, na.rm = T)))
    
         # Rbind and save in dir for later analyses
         avg.hsi.global$region <- "Global"
         avg.hsi.NH$region <- "NH"
         avg.hsi.SH$region <- "SH"
         avg.hsi.all <- rbind(avg.hsi.global, avg.hsi.NH, avg.hsi.SH)
         avg.hsi.all$SDM <- sdm
         avg.hsi.all$ESM <- esm
         
         setwd("/net/kryo/work/fabioben/GODLY/data/shift_metrics")
         save(x = avg.hsi.all, file = paste("table_fut_mean_ann_hsi_species_",sdm,"_",esm,"_30_11_23.RData", sep = ""))
         
         rm(ann2,tab,avg.hsi.all,avg.hsi.SH,avg.hsi.NH,avg.hsi.global); gc()
        
    } # eo for loop - esm in ESMs
    
} # eo for loop - sdm in SDMs


### ================================================================================================================
### ================================================================================================================
### ================================================================================================================

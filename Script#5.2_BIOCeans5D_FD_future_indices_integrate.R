
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 15/08/23: R script to extract and integrate the chosen FD indices (those indices based on dbFD, Faith's index, Jaccard's based on either HSI or PA future compositional data data) together © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - load and combine in a single data.frame the mean annual values of the FD indices chosen for the future ocean (SR,Faith,SES.Faith,FEve,FDiv,FDis,RaoQ)

### Latest update: 15/09/23

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

# Faith's and SR
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/Faith")
files <- dir()[grep("2100-2000",dir())] # ; files # should be 180 files because: 12 months x 3SDMs x 5 ESMs
# f <- files[5]
res <- mclapply(files, function(f) {
            d <- get(load(f))
            # extract month and sdm from filename 
            filename <- str_replace_all(f,".Rdata","")
            # unlist(strsplit(x = filename, split = "_", fixed = T))
            d$ESM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[7]
            return(d)
        }, mc.cores = 5
) # eo mclapply
tab <- bind_rows(res)
rm(res); gc()
# Add cell id
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) # length(unique(tab$cell_id))

### Compute mean annual indices
ann.sr.faith <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.faith = mean(Faith, na.rm = T), sd.faith = sd(Faith, na.rm = T),
            mean.rich = mean(SR, na.rm = T), sd.rich = sd(SR, na.rm = T))
)
# dim(ann.sr.faith); summary(ann.sr.faith)

# Rotate to -179.5/+179.5
ann.sr.faith$x2 <- ann.sr.faith$x
ann.sr.faith[ann.sr.faith$x > 179.5,"x2"] <- (ann.sr.faith[ann.sr.faith$x > 179.5,"x"])-360
# summary(ann.sr.faith$x2)

ggplot() + geom_raster(aes(x = x2, y = y, fill = mean.faith), data = ann.sr.faith) + scale_fill_viridis(name = "Mean annual\nFR") +
      geom_contour(colour = "grey30", binwidth = .1, size = 0.25, aes(x = x2, y = y, z = mean.faith), data = ann.sr.faith) +
      geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
      coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
      panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
      scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
      scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
# Same but with SR
ggplot() + geom_raster(aes(x = x2, y = y, fill = mean.rich), data = ann.sr.faith) + scale_fill_viridis(name = "Mean annual\nSR") +
      geom_contour(colour = "grey30", binwidth = 20, size = 0.25, aes(x = x2, y = y, z = mean.rich), data = ann.sr.faith) +
      geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
      coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
      panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
      scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
      scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)


# SES Faith
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/SES")
files <- dir()[grep("2100-2000",dir())] #; files
# f <- files[2]
res <- mclapply(files, function(f) {
            d <- get(load(f))
            # extract month and sdm from filename 
            filename <- str_replace_all(f,".Rdata","")
            # unlist(strsplit(x = filename, split = "_", fixed = T))
            d$ESM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[7]
            return(d)
        }, mc.cores = 5
) # eo mclapply
tab <- bind_rows(res)
rm(res); gc()
# dim(tab); head(tab); summary(tab)
# Add cell id
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) # length(unique(tab$cell_id))
colnames(tab)[1] <- "SR"
colnames(tab)[2] <- "FR"

### --> value to plot = pd.obs.z (and pd.obs.p for pval)

### Compute annual means
ann.ses.faith <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.ses = mean(pd.obs.z, na.rm = T),
            sd.ses = sd(pd.obs.z, na.rm = T),
            freq.pval = sum(pd.obs.p[pd.obs.p < .05])*100
    )
)
# dim(ann.ses.faith); summary(ann.ses.faith); head(ann.ses.faith)
# Rotate to -179.5/+179.5
ann.ses.faith$x2 <- ann.ses.faith$x
ann.ses.faith[ann.ses.faith$x > 179.5,"x2"] <- (ann.ses.faith[ann.ses.faith$x > 179.5,"x"])-360
# summary(ann.sr.faith$x2)

ggplot() + geom_raster(aes(x = x2, y = y, fill = mean.ses), data = ann.ses.faith) +
     scale_fill_gradient2(name = "Mean annual\nSES", low = "#3288bd", mid = "white", high = "#b2182b") +
     geom_contour(colour = "grey30", binwidth = 1, size = 0.25, aes(x = x2, y = y, z = mean.ses), data = ann.ses.faith) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)


# B div indices (Jac/Jtu/Jne) 
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/beta.div") # dir()
files <- dir()[grep("2100-2000",dir())] #; files
res <- mclapply(files, function(f) {
            d <- get(load(f))
            # extract month and sdm from filename 
            filename <- str_replace_all(f,".Rdata","")
            # unlist(strsplit(x = filename, split = "_", fixed = T))
            d$ESM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[7]
            return(d)
        }, mc.cores = 5
) # eo mclapply
tab <- bind_rows(res)
rm(res); gc()
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) # length(unique(tab$cell_id))

# Compute beta ratio for Jtu
tab$ratio <- tab$beta.jtu/tab$beta.jac

### Compute mean annual indices
ann.beta <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
        jac = mean(beta.jac, na.rm = T), jac.std = sd(beta.jac, na.rm = T),
        jtu = mean(beta.jtu, na.rm = T), jtu.std = sd(beta.jtu, na.rm = T),
        jne = mean(beta.jne, na.rm = T), jne.std = sd(beta.jne, na.rm = T),
        beta.ratio = mean(ratio, na.rm = T), beta.ratio.std = sd(ratio, na.rm = T)
    ) 
) # eo ddf
#dim(ann.beta); head(ann.beta)
#summary(ann.beta) 

# Rotate to -179.5/+179.5
ann.beta$x2 <- ann.beta$x
ann.beta[ann.beta$x > 179.5,"x2"] <- (ann.beta[ann.beta$x > 179.5,"x"])-360

ggplot() + geom_raster(aes(x = x2, y = y, fill = beta.ratio), data = ann.beta) + scale_fill_viridis(name = "Mean annual\nßratio") +
     geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x2, y = y, z = beta.ratio), data = ann.beta) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)


### Remaining indices
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid/")
files <- dir()[grep("2100-2000",dir())] #; files
res <- lapply(files, function(f) {
            d <- get(load(f))
            # extract month and sdm from filename 
            filename <- str_replace_all(f,".Rdata","")
            # unlist(strsplit(x = filename, split = "_", fixed = T))
            d$ESM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[7]
            return(d)
        }
) # eo lapply
tab <- bind_rows(res)
rm(res); gc()
tab$cell_id <- factor(paste(tab$x,tab$y,sep = "_")) 

# Compute mean annual indices
ann.dbFD.prob <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            FEve.avg = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
            FDis.avg = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
            RaoQ.avg = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
            FDiv.avg = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T) )
)
#dim(ann.dbFD.prob)
#summary(ann.dbFD.prob)

# Rotate to -179.5/+179.5
ann.dbFD.prob$x2 <- ann.dbFD.prob$x
ann.dbFD.prob[ann.dbFD.prob$x > 179.5,"x2"] <- (ann.dbFD.prob[ann.dbFD.prob$x > 179.5,"x"])-360

# Maps
ggplot() + geom_raster(aes(x = x2, y = y, fill = FEve.avg), data = ann.dbFD.prob) + scale_fill_viridis(name = "Mean annual\nFEve") +
     geom_contour(colour = "grey30", binwidth = .05, size = 0.25, aes(x = x2, y = y, z = FEve.avg), data = ann.dbFD.prob) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

ggplot() + geom_raster(aes(x = x2, y = y, fill = FDis.avg), data = ann.dbFD.prob) + scale_fill_viridis(name = "Mean annual\nFDis") +
     geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x2, y = y, z = FDis.avg), data = ann.dbFD.prob) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

ggplot() + geom_raster(aes(x = x2, y = y, fill = RaoQ.avg), data = ann.dbFD.prob) + scale_fill_viridis(name = "Mean annual\nRaoQ") +
     geom_contour(colour = "grey30", binwidth = .05, size = 0.25, aes(x = x2, y = y, z = RaoQ.avg), data = ann.dbFD.prob) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

ggplot() + geom_raster(aes(x = x2, y = y, fill = FDiv.avg), data = ann.dbFD.prob) + scale_fill_viridis(name = "Mean annual\nFDiv") +
     geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x2, y = y, z = FDiv.avg), data = ann.dbFD.prob) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

### All FD indices check out for now - compute ∆FD indices per month and ESM and SDM (see Script#5.3)

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
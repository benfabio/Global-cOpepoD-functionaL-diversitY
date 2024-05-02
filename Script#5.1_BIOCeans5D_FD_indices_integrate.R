
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 07/08/23: R script to extract and integrate the chosen FD indices (those indices based on dbFD, Faith's index, Jaccard's and funrar, based on either HSI or PA data) together © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - load and combine in a single data.frame the mean annual values of the FD indices chosen (SR,Faith,SES Faith,FEve,FDiv,FDis,RaoQ)
# - load and combine with the BCP variables (NPP, C export, PSD, e ratio etc.)
# - combine with mesozooplankton biomass fields from MAREDAT (tuned SDMs) made by Nielja K. and Corentin C.

### Latest update: 22/09/23 (adding Corentin's mesozooplankton biomass fields)

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
# x SR (taxonomic richness) based on PA 
# x Faith's index (functional richness based on length of branches of functional dist matrix) based on PA data
# x SES of Faith
# x Functional beta.div based on PA data 
# x FRic (functional hypervolume) based on PA data
# x FEve (functional eveness) based on HSI data
# x FDis (functional dispersion) based on HSI data
# x Rao'Q (dispersion too, through quadratic entropy) based on HSI data
# x FDiv (divergence) based on based on HSI data

# Faith's and SR
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/Faith")
files <- dir()[grep("baseline",dir())] ; files
res <- mclapply(files, function(f) {
            d <- get(load(f))
            # extract month and sdm from filename 
            filename <- str_replace_all(f,".Rdata","")
            # unlist(strsplit(x = filename, split = "_", fixed = T))
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
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
# ggplot() + geom_raster(aes(x = x, y = y, fill = mean.rich), data = ann.sr.faith) + scale_fill_viridis(name = "Mean annual\nSR") +
#      geom_contour(colour = "grey30", binwidth = 10, size = 0.25, aes(x = x, y = y, z = mean.rich), data = ann.sr.faith) +
#      geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#      coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#      panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#      scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#      scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)



# SES Faith
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/SES")
files <- dir()[grep("baseline",dir())] ; files
# f <- files[2]
res <- mclapply(files, function(f) {
            d <- get(load(f))
            # extract month and sdm from filename 
            filename <- str_replace_all(f,".Rdata","")
            # unlist(strsplit(x = filename, split = "_", fixed = T))
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
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
# ggplot() + geom_raster(aes(x = x, y = y, fill = mean.ses), data = ann.ses.faith) +
#     scale_fill_gradient2(name = "Mean annual\nSES", low = "#3288bd", mid = "white", high = "#b2182b") +
#     geom_contour(colour = "grey30", binwidth = 1, size = 0.25, aes(x = x, y = y, z = mean.ses), data = ann.ses.faith) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)


# B div indices (Jac/Jtu/Jne) 
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/beta.div")
files <- dir()[grep("baseline",dir())]
res <- mclapply(files, function(f) {
            d <- get(load(f))
            # extract month and sdm from filename 
            filename <- str_replace_all(f,".Rdata","")
            # unlist(strsplit(x = filename, split = "_", fixed = T))
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
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
# dim(ann.beta); head(ann.beta)
# summary(ann.beta)
# ggplot() + geom_raster(aes(x = x, y = y, fill = beta.ratio), data = ann.beta) + scale_fill_viridis(name = "Mean annual\nßratio") +
#     geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = beta.ratio), data = ann.beta) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#


# FRic (non essential)
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/PA_Gawdis_PCoA_Euclid/") # dir() # should be of length 36
files <- dir()[grep("FDindices_baseline",dir())] ; files
res <- mclapply(files, function(f) {
            d <- get(load(f))
            # extract month and sdm from filename 
            filename <- str_replace_all(f,".Rdata","")
            # unlist(strsplit(x = filename, split = "_", fixed = T))
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
            return(d)
        }, mc.cores = 5
) # eo mclapply
tab <- bind_rows(res)
rm(res); gc()
# Add cell id
tab$cell_id <- factor(paste(tab$x,tab$y,sep = "_")) 
# Compute mean annual indices
ann.dbFD <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            FRic.avg = mean(FRic, na.rm = T), FRic.std = sd(FRic, na.rm = T),
            FEve.avg = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
            FDis.avg = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
            RaoQ.avg = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
            FDiv.avg = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T) )
)
# dim(ann.dbFD); head(ann.dbFD)
# summary(ann.dbFD)
# ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.avg), data = ann.dbFD) + scale_fill_viridis(name = "Mean annual\nFRic") +
#     geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = FRic.avg), data = ann.dbFD) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)


# Remaining indices
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid/")
files <- dir()[grep("_baseline",dir())]; files
res <- lapply(files, function(f) {
            d <- get(load(f))
            # extract month and sdm from filename 
            filename <- str_replace_all(f,".Rdata","")
            # unlist(strsplit(x = filename, split = "_", fixed = T))
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
            return(d)
        }
) # eo lapply
tab <- bind_rows(res)
rm(res); gc()
tab$cell_id <- factor(paste(tab$x,tab$y,sep = "_")) 
# Compute mean annual indices
ann.dbFD.prob <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            FRic.avg = mean(FRic, na.rm = T), FRic.std = sd(FRic, na.rm = T),
            FEve.avg = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
            FDis.avg = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
            RaoQ.avg = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
            FDiv.avg = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T) )
)
# dim(ann.dbFD.prob)
# summary(ann.dbFD.prob)
# ggplot() + geom_raster(aes(x = x, y = y, fill = FEve.avg), data = ann.dbFD.prob) + scale_fill_viridis(name = "Mean annual\nFEve") +
#     geom_contour(colour = "grey30", binwidth = .05, size = 0.25, aes(x = x, y = y, z = FEve.avg), data = ann.dbFD.prob) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)


### Cbind all indices together 
# dim(ann.sr.faith); dim(ann.ses.faith); dim(ann.beta); dim(ann.dbFD); dim(ann.dbFD.prob) # all same dimensions 
# head(ann.sr.faith$cell_id); head(ann.ses.faith$cell_id); head(ann.beta$cell_id); head(ann.dbFD$cell_id); head(ann.dbFD.prob$cell_id)

all.fd <- data.frame(
    cell_id = ann.sr.faith$cell_id, x = ann.sr.faith$x, y = ann.sr.faith$y,
    SR = ann.sr.faith$mean.rich, Faith = ann.sr.faith$mean.faith, SES.Faith = ann.ses.faith$mean.ses,
    Jac = ann.beta$jac, Jtu = ann.beta$jtu, Jne = ann.beta$jne, Bratio = ann.beta$beta.ratio,
    FEve = ann.dbFD.prob$FEve.avg, FDis = ann.dbFD.prob$FDis.avg, RaoQ = ann.dbFD.prob$RaoQ.avg, FDiv = ann.dbFD.prob$FDiv.avg, 
    FRic =  ann.dbFD$FRic.avg
) # eo ddf
# dim(all.fd) ; summary(all.fd)
# Re-scale Raos'Q
all.fd$RaoQ.scaled <- all.fd$RaoQ/max(all.fd$RaoQ)

### Perform a PCA and map components to assess main spatial patterns
require("FactoMineR")
pca <- PCA(X = all.fd[,c(4:14)], scale.unit = T, ncp = 5, graph = F)
# summary(pca)
#                       Dim.1   Dim.2   Dim.3   Dim.4   
# Variance               7.040   2.230   0.908   0.417   
# % of var.             64.003  20.276   8.252   3.794  
# Cumulative % of var.  64.003  84.279  92.532  96.325
# Basic plots
#plot(pca, axes = c(1,2), choix = "var")
#plot(pca, axes = c(2,3), choix = "var")
#plot(pca, axes = c(3,4), choix = "var")
# str(pca$ind$coord)

# Provide PC scores (1-3) to ddf
all.fd$PC1 <- pca$ind$coord[,1]
all.fd$PC2 <- pca$ind$coord[,2]
all.fd$PC3 <- pca$ind$coord[,3]

# Map
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC1), data = all.fd) +
     scale_fill_gradient2(name = "PC1 (64.0%)", low = "#3288bd", mid = "white", high = "#b2182b") +
     geom_contour(colour = "grey30", binwidth = 1, size = 0.25, aes(x = x, y = y, z = PC1), data = all.fd) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC2), data = all.fd) +
     scale_fill_gradient2(name = "PC2 (20.3%)", low = "#3288bd", mid = "white", high = "#b2182b") +
     geom_contour(colour = "grey30", binwidth = 1.5, size = 0.25, aes(x = x, y = y, z = PC2), data = all.fd) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC3), data = all.fd) +
     scale_fill_gradient2(name = "PC3 (8.3%)", low = "#3288bd", mid = "white", high = "#b2182b") +
     geom_contour(colour = "grey30", binwidth = 1, size = 0.25, aes(x = x, y = y, z = PC3), data = all.fd) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

### Rotate x to also have it as 0°-360°
all.fd$x2 <- all.fd$x
all.fd[which(all.fd$x < 0),"x2"] <- (all.fd[which(all.fd$x < 0),"x"]) + 360
# summary(all.fd$x2)

### Save all indices
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/")
save(x = all.fd, file = "table_mean_ann_FD_indices_baseline_04.09.23.RData")

### Save PC maps
library("ggpubr")
panel <- ggarrange(map1,map2,map3, align = "hv", ncol = 1)
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
ggsave(plot = panel, filename = "map_PC1-3_FD_indices.pdf", dpi = 300, height = 10, width = 8)


### ----------------------------------------------------------------

### 2) Retrieve the mean annual values of the chosen biological carbon pump variables to examine together with the FD indices

### To re-load the mean ensemble FD indices
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/")
all.fd <- get(load("table_mean_ann_FD_indices_baseline_04.09.23.RData"))

# x CHL-A --> globcolour_monthly_100km_CHL_REP
# x NPP v1 --> NPP computed from standard VGPM algorithm (updated on the 05.09.23)
# x NPP v2 - DeVries & Weber 2017
# x PSD slope v2 - Clements et al. 2023 (UVP particles data)
# x POC export below euphotic zone (EZ) v1 - DeVries & Weber 2017
# x POC export below euphotic zone (EZ) v2 - Clements et al. 2023
# x e ratio (POC exp/NPP) v1 - DeVries & Weber 2017
# x PSD slope v1 - Kostadinov et al. 2016

### For satellite CHLA (1997-2022) + contribution of phytoplankton PFTs - out puts from CHLA-NPP climatologies preparer
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/CHLA_PFTs_globcolour_100km_2002-2022/")
files <- dir()[grep("clim_ann",dir())] ; files
# f <- files[3]
res <- lapply(files, function(f) {
    
        d <- get(load(f))
        ff <- str_replace_all(f,".Rdata","")
        
        # Change colname #4 to var of interest
	    if(f == files[1]) {
		    colnames(d)[4] <- substring(ff,43,47) 
	    } else {
		    colnames(d)[4] <- substring(ff,43,48) 
	    } # eo if else loop
        
        return(d)
    }
    
) # eo LAPPLY
# Cbind
sat <- bind_cols(res)
colnames(sat)[c(1:3)] <- c("id","x","y")
sat <- sat[,c(1:4,8,12,16,20,24,28,32,36,40)]
### Add those mean annual values to 'all.df'
# v <- "CHL"
for(v in colnames(sat)[c(4:length(sat))] ) {
    
    message(paste("Adding ",v, sep = ""))
    spg <- sat[,c("x","y",v)]
    colnames(spg)[1] <- "x"
    coordinates(spg) <- ~ x + y
    gridded(spg) <- TRUE
    ras <- raster(spg)
    crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
    all.fd[,v] <- raster::extract(x = ras, y = all.fd[,c("x","y")])
    
}
# Check
# summary(all.fd)
# Guet.

### For satellite NPP (monthly standard VGPM from 2003 to 2022, sourced from http://orca.science.oregonstate.edu/1080.by.2160.monthly.hdf.vgpm.m.chl.m.sst.php)
# Climatologies are stored in: 
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/NPP_VGPM.standard_2002-2023")
npp <- get(load("table_clim_ann_NPP_standard_VGPM_2003-2022.RData"))
# summary(npp)
spg <- npp[,c("x","y","NPP")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # plot(log(ras))
all.fd[,"NPP"] <- raster::extract(x = ras, y = all.fd[,c("x","y")])
rm(sat,npp);gc()

### DeVries & Weber (2017) --> extract NPP + POC export below the EZ --> compute e-ratio

### NOTE: Beware that units here are in mmolC/m^2/YEAR ! NOT mgC/m2/day
### NOTE: Beware that cell grid is 2°x2°

setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/DeVries&Weber_2017")
nc <- nc_open("Cexp_deVries_2017.nc")
lat <- raster::stack("Cexp_deVries_2017.nc", varname = "LAT")
lon <- raster::stack("Cexp_deVries_2017.nc", varname = "LON")
x <- as.data.frame(lon, xy = T)
y <- as.data.frame(lat, xy = T)
nc_close(nc)

vars <- c("NPP","FPOCex","FPOC100m","POCflux") 
# - Net primary production (mmolC/m^2/yr)
# - Sinking POC export at base of euphotic zone (mmolC/m^2/yr)
# - Sinking POC export at 100 m (mmolC/m^2/yr)
# - Sinking POC flux (mmolC/m^2/yr)

# v <- "NPP"
clims <- mclapply(vars, function(v) {
    
				# Get data from the nc file
				ras <- raster::stack("Cexp_deVries_2017.nc", varname = v)
				d <- as.data.frame(ras, xy = T)
				d$x <- x[,3]
				d$y <- y[,3]
                # Add cell ID and compute average 
                d$id <- paste(d$x, d$y, sep = "_") # length( unique(d$id) )
                # Melt and compute clim
                m.d <- melt(d, id.vars = c("id","x","y"))
                # 0 should be NA
                m.d[,"value"] <- na_if(m.d[,"value"], 0) # summary(m.d)
				# Compute multi-model average (12 models)
				clim <- data.frame( m.d %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), mean = mean(value, na.rm = T)) )
                clim <- clim[order(clim$id),]
                # summary(clim)
                # Change colnames 
				colnames(clim)[length(clim)] <- v
                return(clim) 
        
		}, mc.cores = length(vars)
        
) # eo lapply
# Cbind
poc <- bind_cols(clims)
poc <- poc[,c(1:4,8,12,16)]
colnames(poc)[c(1:3)] <- c("id","x","y")
poc$eratio <- (poc$FPOCex)/(poc$NPP)
# summary(poc)
# Remove values < 0 (irrealistic)
poc <- poc[which(poc$POCflux > 0 & poc$FPOCex > 0 & poc$FPOC100m > 0 & poc$eratio > 0 & poc$eratio < .5),]


# Next, convert to raster and combine with FD indices
spg <- poc[,c("x","y","NPP")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
all.fd$NPPv2 <- raster::extract(x = ras, y = all.fd[,c("x2","y")])

# Same with FPOCex below euphotic zone
spg <- poc[,c("x","y","FPOCex")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
all.fd$FPOCex <- raster::extract(x = ras, y = all.fd[,c("x2","y")])

# Same with POCflux below euphotic zone
spg <- poc[,c("x","y","POCflux")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
all.fd$POCflux <- raster::extract(x = ras, y = all.fd[,c("x2","y")])

### And e-ratio
spg <- poc[,c("x","y","eratio")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
all.fd$eratio <- raster::extract(x = ras, y = all.fd[,c("x2","y")])

### Convert "NPP","FPOCex","POCflux" and "POCflux" (which are in mmolesC/m2/year) to mgC/m2.day  
### No need for e-ratio since it is a ratio
all.fd$NPPv2 <- ((all.fd$NPPv2)*12.0107)/365
all.fd$POCflux <- ((all.fd$POCflux)*12.0107)/365
all.fd$FPOCex <- ((all.fd$FPOCex)*12.0107)/365

### Check
# summary(all.fd)
# sanity check for NPP data
#ggplot(aes(x = NPP, y = NPPv2), data = all.fd) + geom_point() + theme_classic()
#ggplot(aes(x = log10(NPP), y = log10(NPPv2)), data = all.fd) + geom_point() + theme_classic()
### Makes sense. Keep log of NPPv2 since less NaN

### For more POC stuff: Clements et al. 2023 (UVP particles data) --> extract POC export below the EZ and PSD slope + compute e ratio from updated NPP (in mgC/m2/d)
# Biovolume = Pred_BV = particle biovolume (ppm)
# Slope = Pred_slope = slope of particle size distribution (unitless)
# Export = carbon export flux from the seasonal euphotic zone (mgC/m2/d)

setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/Clements&al.2023")
clems <- get(load("table_clim_mon_part_BV_slope_carbon_Clements&al.2023_07_08_23.RData"))
colnames(clems)[1] <- "id"
clim <- data.frame( clems %>% group_by(id,variable) %>% summarise(x = unique(x), y = unique(y), annual = mean(mean, na.rm = T)) )
# Dcast to have all three variables as columns
d.clim <- dcast(data = clim, id+x+y ~ variable, value.var = "annual")
# head(d.clim); dim(d.clim)

spg <- d.clim[,c("x","y","Biovolume")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
all.fd$Biovolume <- raster::extract(x = ras, y = all.fd[,c("x","y")])
summary(all.fd$Biovolume)

spg <- d.clim[,c("x","y","Slope")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
all.fd$Slope <- raster::extract(x = ras, y = all.fd[,c("x","y")])
summary(all.fd$Slope)

spg <- d.clim[,c("x","y","Export")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
all.fd$Export <- raster::extract(x = ras, y = all.fd[,c("x","y")])
# summary(all.fd$Export)


### For PSD slope: Kostadinov et al. (2009)
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/Kostadinov&al._2009")
psd <- read.table("table_annual_PSDslope_Kostadinov_1d_29_01_20.txt", h = T, sep = "\t")
# Check if coords are consistent with those in 'div'
psd$x2 <- (psd$x2)+0.5
spg <- psd[,c("x2","y","slope")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
all.fd$Slope2 <- raster::extract(x = ras, y = all.fd[,c("x2","y")])

### 2nd sanity check for slopes this time: 
# ggplot(aes(x = Slope, y = Slope2), data = all.fd) + geom_point(alpha = .1) + geom_abline(slope = 1, intercept = 0, colour = "red") + theme_classic()
### Note: slopes from UVP data are way steeper than based on satellite data. Let's examine their spatial pattern.

### 3rd sanity check: CHL ~ NPP
#ggplot(aes(x = log10(CHL), y = log10(NPP)), data = all.fd) + geom_point(alpha = .1) + theme_classic()

# ggplot() + geom_raster(aes(x = x, y = y, fill = Slope), data = all.fd) +
#     geom_contour(colour = "grey70", binwidth = .25, size = 0.25, aes(x = x, y = y, z = Slope), data = all.fd) +
#     scale_fill_viridis(name = "PSD slope") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
# ggplot() + geom_raster(aes(x = x, y = y, fill = Slope2), data = all.fd) +
#     geom_contour(colour = "grey70", binwidth = .25, size = 0.25, aes(x = x, y = y, z = Slope2), data = all.fd) +
#     scale_fill_viridis(name = "PSD slope (v2)") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

### CCL: Data from Clements et al. (2023) has way more gaps (14% of cells have NAs). So use it for sensitivity analyses.

### Save 
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/")
save(x = all.fd, file = "table_mean_ann_FD_indices_baseline+BCP_06.09.23.RData")


### 22/11/23: Adding global mesozooplankton biomass fields from Nielja and Corentin
all.fd <- get(load("table_mean_ann_FD_indices_baseline+BCP_06.09.23.RData"))
mesozoo <- read.csv("MAREDAT_MESOZOOPLANKTON_TUNED_SDMs_CORENTIN_22_09_23.csv", h = T)
# dim(mesozoo); str(mesozoo)
# summary(mesozoo) # right coordinates

# Convert mesozoo$target back to original dimensions (reverse log10 transform) - values already in mgC/m3
mesozoo$biomass <- (10^(mesozoo$target)-1)
# summary(mesozoo$biomass) # max is 1.37900 mgC/m3

# Compute annual mean
mesozoo$id <- factor(paste(mesozoo$Longitude, mesozoo$Latitude, sep = "_")) # unique(mesozoo$id)
ann.meso <- data.frame(mesozoo %>% group_by(id) %>% summarise(x = unique(Longitude), y = unique(Latitude), biomass = mean(biomass, na.rm = T)))
# summary(ann.meso)

# Quick map to check
# ggplot() + geom_raster(aes(x = x, y = y, fill = biomass), data = ann.meso) +
#     geom_contour(colour = "grey70", binwidth = 1.5, size = 0.25, aes(x = x, y = y, z = biomass), data = ann.meso) +
#     scale_fill_viridis(name = "Mean annual mesozooplankton\nbiomass (mgC/m3)") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

# Pattern looks good. SDM pipeline + annual averaging means values ly range between 0 and 10 mgC/m3
# Add to 'all.fd'
spg <- ann.meso[,c("x","y","biomass")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg) # plot(ras)
all.fd$MESOZOO <- raster::extract(x = ras, y = all.fd[,c("x","y")])
# Save
save(x = all.fd, file = "table_mean_ann_FD_indices_baseline+BCP+biom_22.11.23.RData")

### Examine covariance quickly (code below)
ggplot(data = all.fd, aes(x = SR, y = MESOZOO)) + geom_point(alpha = .1) +
    xlab("Mean annual copepod richness") + ylab("Mean annual mesozooplankton biomass") +
    theme_bw()

ggplot(data = all.fd, aes(x = Jtu, y = MESOZOO)) + geom_point(alpha = .1) +
    xlab("Mean annual turnover") + ylab("Mean annual mesozooplankton biomass") +
    theme_bw()

ggplot(data = all.fd, aes(x = SR, y = MESOZOO, colour = Jtu)) + geom_point(alpha = .1) +
    xlab("Mean annual copepod richness") + ylab("Mean annual mesozooplankton biomass") +
    theme_bw()

ggplot(data = all.fd, aes(x = Bratio, y = MESOZOO)) + geom_point(alpha = .1) +
    xlab("Mean annual Bratio") + ylab("Mean annual mesozooplankton biomass") +
    theme_bw()

ggplot(data = all.fd, aes(x = FEve, y = MESOZOO, colour = Jac)) + geom_point(alpha = .1) +
    scale_colour_viridis(name = "Trait dissimilarity", option = "C") + 
    xlab("Mean annual FEve") + ylab("Mean annual mesozooplankton biomass") +
    theme_bw()

ggplot(data = all.fd, aes(x = FDis, y = MESOZOO)) + geom_point(alpha = .1) +
    xlab("Mean annual FDis") + ylab("Mean annual mesozooplankton biomass") +
    theme_bw()

ggplot(data = all.fd, aes(x = RaoQ.scaled, y = MESOZOO)) + geom_point(alpha = .1) +
    xlab("Mean annual Rao's Q") + ylab("Mean annual mesozooplankton biomass") +
    theme_bw()
#
ggplot(data = all.fd, aes(x = FDiv, y = MESOZOO)) + geom_point(alpha = .1) +
    xlab("Mean annual FDiv") + ylab("Mean annual mesozooplankton biomass") +
    theme_bw()

### Evidence an emergent negatove relationship (likely significant) between copepod diversity and biomass! 

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 3) Examine covariance of BCP variables in a PCA - like you did for FD indices
#colnames(all.fd)
data4pca <- na.omit(all.fd[,c(1:3,21:35,39,40)])
# log10 transform all vars except slope and eratio; colnames(data4pca)
data4pca[,c(4:17,20)] <- log10(data4pca[,c(4:17,20)])

### Perform a PCA and map components to assess main spatial patterns
require("FactoMineR")
pca <- PCA(X = data4pca[,c(4:length(data4pca))], scale.unit = T, ncp = 4, graph = F)
# summary(pca)
#                       Dim.1   Dim.2   Dim.3   Dim.4   Dim.5
# Variance              9.404   5.446   1.040   0.432   0.260
# % of var.             55.315  32.035   6.117   2.542  1.529 
# Cumulative % of var.  55.315  87.350  93.467  96.009  97.53

# Basic plots
#plot(pca, axes = c(1,2), choix = "var")
#plot(pca, axes = c(2,3), choix = "var")
#plot(pca, axes = c(3,4), choix = "var")
# str(pca$ind$coord)

# Provide PC scores (1-3) to ddf
data4pca[,c("PC1","PC2","PC3")] <- pca$ind$coord[,c(1:3)]

# Map
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC1), data = data4pca) +
     scale_fill_gradient2(name = "PC1 (55%)", low = "#3288bd", mid = "white", high = "#b2182b") +
     geom_contour(colour = "grey30", binwidth = 2, size = 0.25, aes(x = x, y = y, z = PC1), data = data4pca) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC2), data = data4pca) +
     scale_fill_gradient2(name = "PC2 (32%)", low = "#3288bd", mid = "white", high = "#b2182b") +
     geom_contour(colour = "grey30", binwidth = 2, size = 0.25, aes(x = x, y = y, z = PC2), data = data4pca) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC3), data = data4pca) +
     scale_fill_gradient2(name = "PC3 (6%)", low = "#3288bd", mid = "white", high = "#b2182b") +
     geom_contour(colour = "grey30", binwidth = 1, size = 0.25, aes(x = x, y = y, z = PC3), data = data4pca) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

panel <- ggarrange(map1,map2,map3, align = "hv", ncol = 1)
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
ggsave(plot = panel, filename = "map_PC1-3_BCP_indices_22_11_23.pdf", dpi = 300, height = 10, width = 8)


### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
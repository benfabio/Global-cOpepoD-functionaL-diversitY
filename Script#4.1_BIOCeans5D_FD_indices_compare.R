
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 10/08/23: R script to compare FD indices based on HSI to those based on PA data © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - load and combine in a single data.frame the mean annual values of the various FD indices
# - examine covariance/overlap

### Latest update: 11/08/23

### ------------------------------------------------------------------------------------------------------------------------------------------------------

library("marmap")
library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("viridis")
library("parallel")
library("ggpubr")
library("ggthemes")

world <- map_data("world") # coastlines for maps

setwd("/net/kryo/work/fabioben/GODLY/data") # working dir

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### Choose and retrieve the mean annual values of the chosen indices to compare:
# - Standardized FRic VS. Faith's index (for functional richness)
# - PA_based FEve/RaoQ/FDis/FDiv VS. HSI-based FEve/RaoQ/FDis/FDiv
# - PA_based funrar VS. HSI-based funrar

### --------------------------------------------------------------------

### 1) Standardized FRic VS. Faith's index
setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/Faith/Faith")
files <- dir()[grep("_baseline_",dir())]
#f <- files[1]
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
# dim(tab); head(tab); summary(tab)
# Add cell id
tab$cell_id <- factor(paste(tab$x,tab$y,sep = "_")) # length(unique(tab$cell_id))
ann.faith <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.faith = mean(Faith, na.rm = T), sd.faith = sd(Faith, na.rm = T),
            mean.rich = mean(SR, na.rm = T), sd.rich = sd(SR, na.rm = T))
)
# dim(ann.faith) ; summary(ann.faith)
# Draw maps
# ggplot() + geom_raster(aes(x = x, y = y, fill = mean.faith), data = ann.faith) + scale_fill_viridis(name = "Mean annual\nFaith's index") +
#     geom_contour(colour = "grey30", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = mean.faith), data = ann.faith) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)


setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/db_FD/PA_Gawdis_PCoA_Euclid/Stand.Fric") 
files <- dir()[grep("FDindices_baseline",dir())]
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
# head(tab); summary(tab)
# Add cell id
tab$cell_id <- factor(paste(tab$x,tab$y,sep = "_")) # length(unique(tab$cell_id))

### Compute mean annual indices
ann.ind.pa <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            FRic.avg = mean(FRic, na.rm = T), FRic.std = sd(FRic, na.rm = T),
            FEve.avg = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
            FDis.avg = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
            RaoQ.avg = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
            FDiv.avg = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T) )
)
# dim(ann.ind.pa); summary(ann.ind.pa) # same nrow as ann.faith
# ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.avg), data = ann.ind.pa) + scale_fill_viridis(name = "Mean annual\nFRic") +
#     geom_contour(colour = "grey30", binwidth = .05, size = 0.25, aes(x = x, y = y, z = FRic.avg), data = ann.ind.pa) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

ann.faith$FRic <- ann.ind.pa$FRic.avg
cor(ann.faith$mean.faith, ann.faith$FRic, method = "spearman") # 0.77
# Standardize them by their max values
ann.faith$Faith.std <- (ann.faith$mean.faith)/max(ann.faith$mean.faith, na.rm = T)
ann.faith$FRic.std <- (ann.faith$FRic)/max(ann.faith$FRic, na.rm = T) # won't change the correlation coeff
# Bi-plot
ggplot(data = ann.faith, aes(x = Faith.std, y = FRic.std, colour = abs(y))) + geom_point() + scale_colour_viridis() + 
    geom_smooth(colour = "red", method = "lm") + xlab("Mean annual Faith's index") + ylab("Mean annual FRic") + 
    theme_classic()
    
### Fairly good correspondence
summary(lm(FRic.std ~ Faith.std, data = ann.faith))
# Multiple R-squared:  0.604,    Adjusted R-squared:  0.604
# F-statistic: 5.441e+04 on 1 and 35678 DF,  p-value: < 2.2e-16

### Per biome
# Tropics
cor(ann.faith[which(abs(ann.faith$y) < 30),"Faith.std"], ann.faith[which(abs(ann.faith$y) < 30),"FRic.std"], method = "spearman") # 0.94
summary(lm(FRic.std ~ Faith.std, data = ann.faith[which(abs(ann.faith$y) < 30),])) # R-squared: 0.899

# High latitudes
cor(ann.faith[which(abs(ann.faith$y) > 60),"Faith.std"], ann.faith[which(abs(ann.faith$y) > 60),"FRic.std"], method = "spearman") # 0.875
summary(lm(FRic.std ~ Faith.std, data = ann.faith[which(abs(ann.faith$y) > 60),])) # R-squared: 0.75

# Temperate latitudes
cor(ann.faith[which(abs(ann.faith$y) < 60 & abs(ann.faith$y) > 30),"Faith.std"], ann.faith[which(abs(ann.faith$y) < 60 & abs(ann.faith$y) > 30),"FRic.std"], method = "spearman")
# 0.76
summary(lm(FRic.std ~ Faith.std, data = ann.faith[which(abs(ann.faith$y) < 60 & abs(ann.faith$y) > 30),])) # R-squared: 0.628

# Departure from correlation outside the tropics because the GLM_bias in the high latitude FRic


### --------------------------------------------------------------------

### 2) PA_based FEve/RaoQ/FDis/FDiv VS. HSI-based FEve/RaoQ/FDis/FDiv
setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid")
files <- dir()[grep("FDindices_baseline",dir())]
# f <- files[2]
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
# dim(tab); head(tab); summary(tab)
tab$cell_id <- factor(paste(tab$x,tab$y,sep = "_")) # length(unique(tab$cell_id))
# Compute mean annual average
ann.ind.hsi <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            FRic.hsi = mean(FRic, na.rm = T), FRic.std = sd(FRic, na.rm = T),
            FEve.hsi = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
            FDis.hsi = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
            RaoQ.hsi = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
            FDiv.hsi = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T) )
) # eo summarise 
# dim(ann.ind.hsi)
ann.ind.hsi$FRic.pa <- ann.ind.pa$FRic.avg
ann.ind.hsi$FEve.pa <- ann.ind.pa$FEve.avg
ann.ind.hsi$FDis.pa <- ann.ind.pa$FDis.avg
ann.ind.hsi$RaoQ.pa <- ann.ind.pa$RaoQ.avg
ann.ind.hsi$FDiv.pa <- ann.ind.pa$FDiv.avg

# Compute correlations coeff
cor(ann.ind.hsi$FRic.pa, ann.ind.hsi$FRic.hsi, method = "spearman") # 0.30
summary(lm(FRic.pa ~ FRic.hsi, data = ann.ind.hsi)) # R-squared:  0.075

cor(ann.ind.hsi$FEve.pa, ann.ind.hsi$FEve.hsi, method = "spearman") # 0.61
summary(lm(FEve.pa ~ FEve.hsi, data = ann.ind.hsi)) # R-squared:  0.298

cor(ann.ind.hsi$FDis.pa, ann.ind.hsi$FDis.hsi, method = "spearman") # 0.87
summary(lm(FDis.pa ~ FDis.hsi, data = ann.ind.hsi)) # R-squared:  0.655

cor(ann.ind.hsi$RaoQ.pa, ann.ind.hsi$RaoQ.hsi, method = "spearman") # 0.77
summary(lm(RaoQ.pa ~ RaoQ.hsi, data = ann.ind.hsi)) # R-squared:  0.579

cor(ann.ind.hsi$FDiv.pa, ann.ind.hsi$FDiv.hsi, method = "spearman") # 0.28
summary(lm(FDiv.pa ~ FDiv.hsi, data = ann.ind.hsi)) # R-squared:  0.104

# Bi-plots
p1 <- ggplot(data = ann.ind.hsi, aes(x = FRic.hsi, y = FRic.pa, colour = abs(y))) + geom_point() + 
    geom_smooth(colour = "red", method = "lm") + scale_colour_viridis(name = "Latitude", direction = -1) + theme_classic()

p2 <- ggplot(data = ann.ind.hsi, aes(x = FEve.hsi, y = FEve.pa, colour = abs(y))) + geom_point() + 
    geom_smooth(colour = "red", method = "lm") + scale_colour_viridis(name = "Latitude", direction = -1) + theme_classic()

p3 <- ggplot(data = ann.ind.hsi, aes(x = FDis.hsi, y = FDis.pa, colour = abs(y))) + geom_point() + 
    geom_smooth(colour = "red", method = "lm") + scale_colour_viridis(name = "Latitude", direction = -1) + theme_classic()

p4 <- ggplot(data = ann.ind.hsi, aes(x = RaoQ.hsi, y = RaoQ.pa, colour = abs(y))) + geom_point() + 
    geom_smooth(colour = "red", method = "lm") + scale_colour_viridis(name = "Latitude", direction = -1) + theme_classic()

p5 <- ggplot(data = ann.ind.hsi, aes(x = FDiv.hsi, y = FDiv.pa, colour = abs(y))) + geom_point() + 
    geom_smooth(colour = "red", method = "lm") + scale_colour_viridis(name = "Latitude", direction = -1) + theme_classic()

ggarrange(p1,p2,p3,p4,p5, align = "hv", ncol = 3, nrow = 2)


# Check FEve without and with tropics
cor(ann.ind.hsi[which(abs(ann.ind.hsi$y) < 20),"FEve.hsi"], ann.ind.hsi[which(abs(ann.ind.hsi$y) < 20),"FEve.pa"], method = "spearman") # 0.46
cor(ann.ind.hsi[which(abs(ann.ind.hsi$y) < 30),"FEve.hsi"], ann.ind.hsi[which(abs(ann.ind.hsi$y) < 30),"FEve.pa"], method = "spearman") # 0.47
cor(ann.ind.hsi[which(abs(ann.ind.hsi$y) < 50),"FEve.hsi"], ann.ind.hsi[which(abs(ann.ind.hsi$y) < 50),"FEve.pa"], method = "spearman") # 0.65
cor(ann.ind.hsi[which(abs(ann.ind.hsi$y) < 60),"FEve.hsi"], ann.ind.hsi[which(abs(ann.ind.hsi$y) < 60),"FEve.pa"], method = "spearman") # 0.67
cor(ann.ind.hsi[which(abs(ann.ind.hsi$y) < 70),"FEve.hsi"], ann.ind.hsi[which(abs(ann.ind.hsi$y) < 70),"FEve.pa"], method = "spearman") # 0.64
# So, better if you remove latitudes lower than 60° but not that much

### --------------------------------------------------------------------

### 3) 'funrar' 
setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/funrar")
files <- dir()[grep("funrar_baseline_",dir())]
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
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) 

### Compute mean annual indices
ann.funrar.hsi <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.dist = mean(distinct, na.rm = T), sd.dist = sd(distinct, na.rm = T),
            mean.scar = mean(scarcity, na.rm = T), sd.scar = sd(scarcity, na.rm = T))
)
# dim(ann.funrar.hsi); summary(ann.funrar.hsi)

setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/funrar/PA") # dir() # should be of length 36
files <- dir()[grep("funrar_baseline_",dir())]
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
# Add cell id
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) 

### Compute mean annual indices
ann.funrar.pa <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.dist = mean(distinct, na.rm = T), sd.dist = sd(distinct, na.rm = T)
    ) # eo summarise
)
# dim(ann.funrar.pa); summary(ann.funrar.pa)

# combine
ann.funrar.hsi$mean.dist.pa <- ann.funrar.pa$mean.dist

# check corr
cor(ann.funrar.hsi$mean.dist.pa, ann.funrar.hsi$mean.dist, method = "spearman") # 0.97
summary(lm(mean.dist.pa ~ mean.dist, data = ann.funrar.hsi)) # R-squared: 0.93

ggplot(data = ann.funrar.hsi, aes(x = mean.dist, y = mean.dist.pa, colour = abs(y))) + geom_point() + 
    geom_smooth(colour = "red", method = "lm") + scale_colour_viridis(name = "Latitude", direction = -1) + 
    theme_classic()

### Same with HSI or PA

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
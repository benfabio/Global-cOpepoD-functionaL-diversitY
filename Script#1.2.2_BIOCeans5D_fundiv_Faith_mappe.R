
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 03/07/23: R script to map Faith's index computed based on a Gower distance matrix and the community tables made during the MSc thesis of Jonas Wydler (data from Benedetti et al., 2023 - JBIO) © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - load the tables containing the Faith indices values + the tables of SES Faith
# - compute and map the mean (ensemble) annual indices, the associated stdev
# - compute and map the monthly averages, and the SDM-specific predictions

### Latest update: 03/07/23

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

### 1) Load the tables of Faith's indices and compute ensemble averages and uncertainties around it
setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/Faith/Faith") # dir() # should be of length 36
files <- dir()
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

### Compute mean annual indices
ann <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.faith = mean(Faith, na.rm = T), sd.faith = sd(Faith, na.rm = T),
            mean.rich = mean(SR, na.rm = T), sd.rich = sd(SR, na.rm = T))
)
# dim(ann); summary(ann)

### Draw maps
p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.faith), data = ann) + scale_fill_viridis(name = "Mean annual\nFaith's index") +
    geom_contour(colour = "grey30", binwidth = 0.2, size = 0.25, aes(x = x, y = y, z = mean.faith), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.faith), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = sd.faith), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.rich), data = ann) + scale_fill_viridis(name = "Mean annual\nrichness") +
    geom_contour(colour = "grey30", binwidth = 15, size = 0.25, aes(x = x, y = y, z = mean.rich), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.rich), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 10, size = 0.25, aes(x = x, y = y, z = sd.rich), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

panel.maps.ann <- ggarrange(p3,p4,p1,p2, align = 'hv', ncol = 2, nrow = 2)    

# To save individual maps 
ggsave(plot = p1, filename = "map_mean_ann_Faith.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = p3, filename = "map_mean_ann_SR.jpg", dpi = 300, height = 4, width = 7)


### Same as above but monthly
mon <- data.frame(tab %>% group_by(month,cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.faith = mean(Faith, na.rm = T), sd.faith = sd(Faith, na.rm = T),
            mean.rich = mean(SR, na.rm = T), sd.rich = sd(SR, na.rm = T))
)
# dim(mon); summary(mon)
# Re-order months correctly
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
mon$month <- factor(mon$month, months)

p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.faith), data = mon) + scale_fill_viridis(name = "Monthly Faith's\nindex") +
    geom_contour(colour = "grey30", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = mean.faith), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.rich), data = mon) + scale_fill_viridis(name = "Monthly\nFaith's index") +
    geom_contour(colour = "grey30", binwidth = 25, size = 0.25, aes(x = x, y = y, z = mean.rich), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)


### Same as above but per SDM
sdms <- data.frame(tab %>% group_by(SDM,cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.faith = mean(Faith, na.rm = T), sd.faith = sd(Faith, na.rm = T),
            mean.rich = mean(SR, na.rm = T), sd.rich = sd(SR, na.rm = T))
)
# dim(sdms); summary(sdms)
sdms$SDM <- factor(sdms$SDM, c("GLM","GAM","ANN"))

p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.faith), data = sdms) + scale_fill_viridis(name = "Mean annual\nFaith's index") +
    geom_contour(colour = "grey30", binwidth = 0.2, size = 0.25, aes(x = x, y = y, z = mean.faith), data = sdms) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(SDM), nrow = 3, ncol = 1)

p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.rich), data = sdms) + scale_fill_viridis(name = "Mean annual\nrichness") +
    geom_contour(colour = "grey30", binwidth = 20, size = 0.25, aes(x = x, y = y, z = mean.rich), data = sdms) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(SDM), nrow = 3, ncol = 1)
    
    
### Additional interesting plots:
# Faith's index ~ SR + lat
ggplot(data = ann, aes(x = mean.rich, y = mean.faith, colour = abs(y))) + 
    geom_point() + scale_colour_viridis(name = "Latitude", direction = -1) + 
    #geom_smooth(method = "lm", formula = y ~ poly(x,2), se = T, colour = "black") +
    xlab("Mean annual richness") + ylab("Mean annual Faith's index") + theme_classic()

# Split per 'biomes' and compute differences in slopes 
# Test lm between FR and SR
summary(lm(mean.faith ~ mean.rich+abs(y), data = ann)) # 0.732 
summary(lm(mean.faith ~ mean.rich*abs(y), data = ann)) # 0.801 ! --> interac

### Split data into tropical bands ('biomes')
ann$biome <- NA
ann[which(abs(ann$y) < 30),"biome"] <- "Tropical (0-30°)"
ann[which(abs(ann$y) > 30),"biome"] <- "Temperate (30-45°)"
ann[which(abs(ann$y) > 45),"biome"] <- "Polar & Subpolar (>45°)"
summary(factor(ann$biome))

plot <- ggplot(data = ann, aes(x = mean.rich, y = mean.faith, colour = factor(biome))) + geom_point() + 
  scale_colour_manual(name = "Biome", values = c("#085a73","#08acd5","#ff6262")) +
  geom_smooth(aes(x = mean.rich, y = mean.faith), method = "lm", colour = "black", se = T) + 
  ylab("Mean annual Faith's index") + xlab("Mean annual richness") +
  theme_classic() + facet_wrap(.~ factor(biome))
  
ggsave(plot = plot, filename = "plot_FRxSRxBiomes.jpg", dpi = 300, height = 4, width = 8)           
  
### ANCOVA to tets for differences in slope between biomes
summary(lm(mean.faith ~ mean.rich+factor(biome), data = ann)) # 0.79
summary(lm(mean.faith ~ mean.rich*factor(biome), data = ann)) # 0.87!
summary( aov(mean.faith ~ mean.rich*factor(biome), data = ann) )
#                            Df Sum Sq Mean Sq F value Pr(>F)    
# mean.rich                   1  894.7   894.7  146340 <2e-16 ***
# factor(biome)               2  416.8   208.4   34090 <2e-16 ***
# mean.rich:factor(biome)     2  139.6    69.8   11417 <2e-16 ***
# Residuals               35674  218.1     0.0 
### --> SIGNIF INTERACTION! 

# Compare without the interaction to see if it affects the fit of the linear model
mod1 <- aov(mean.faith ~ mean.rich*factor(biome), data = ann) 
mod2 <- aov(mean.faith ~ mean.rich+factor(biome), data = ann) 
anova(mod1,mod2)
#   Res.Df    RSS Df Sum of Sq     F    Pr(>F)    
#1  35674 218.10                                 
#2  35676 357.71 -2    -139.6 11417 < 2.2e-16 ***
### The anova clearly shows that removing the interaction DOES significantly affect the fit of the model 

# Check linear relationships
# summary(lm(mean.faith ~ mean.rich, data = ann[ann$biome == "Tropical (0-30°)",])) # R-squared: 0.8534
# summary(lm(mean.faith ~ mean.rich, data = ann[ann$biome == "Temperate (30-45°)",])) # R-squared: 0.7054
# summary(lm(mean.faith ~ mean.rich, data = ann[ann$biome == "Polar & Subpolar (>45°)",])) # R-squared: 0.9521


# Faith ~ STDEV
# ggplot(data = ann, aes(x = mean.faith, y = sd.faith, colour = abs(y))) +
#     geom_point() + scale_colour_viridis(name = "Latitude", direction = -1) +
#     xlab("Mean annual Faith's index") + ylab("Uncertainty (stdev)") + theme_classic()
#
# # SR ~ STDEV
# ggplot(data = ann, aes(x = mean.rich, y = sd.rich, colour = abs(y))) +
#     geom_point() + scale_colour_viridis(name = "Latitude", direction = -1) +
#     xlab("Mean annual richness") + ylab("Uncertainty (stdev)") + theme_classic()


### ------------------------------------------------------------------

### 2) Same as above but for SES Faith
# ?ses.pd
# ntaxa: Number of taxa in community
# pd.obs: Observed PD in community
# pd.rand.mean: Mean PD in null communities
# pd.rand.sd: Standard deviation of PD in null communities
# pd.obs.rank: Rank of observed PD vs. null communities
# pd.obs.z: Standardized effect size of PD vs. null communities (= (pd.obs - pd.rand.mean) / pd.rand.sd)
# pd.obs.p: P-value (quantile) of observed PD vs. null communities (= mpd.obs.rank / runs + 1)

setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/Faith/SES") # dir() # should be of length 36
files <- dir()
# f <- files[2]
res <- mclapply(files, function(f) {
            d <- get(load(f))
            # extract month and sdm from filename 
            filename <- str_replace_all(f,".Rdata","")
            # unlist(strsplit(x = filename, split = "_", fixed = T))
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
            return(d)
        }, mc.cores = length(files)
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
### Compute frequency of pd.obs.p < .05 (more meaningful than mean pval)
ann <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.rich = mean(SR, na.rm = T), sd.rich = sd(SR, na.rm = T),
            mean.faith = mean(FR, na.rm = T), sd.faith = sd(FR, na.rm = T),
            mean.ses = mean(pd.obs.z, na.rm = T), sd.ses = sd(pd.obs.z, na.rm = T),
            freq.pval = sum(pd.obs.p[pd.obs.p < .05])*100
    )
)
# dim(ann); summary(ann)

# First, check if FR matches maps from above
# ggplot() + geom_raster(aes(x = x, y = y, fill = mean.faith), data = ann) + scale_fill_viridis(name = "Mean annual\nFaith's index") +
#      geom_contour(colour = "grey30", binwidth = 0.15, size = 0.25, aes(x = x, y = y, z = mean.faith), data = ann) +
#      geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#      coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#      panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#      scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#      scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
### --> matches previous maps above, focus on SES and pvalues
        
### Draw maps
p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.ses), data = ann) +
    scale_fill_gradient2(name = "Mean annual\nSES", low = "#3288bd", mid = "white", high = "#b2182b") +
    geom_contour(colour = "grey30", binwidth = 1, size = 0.25, aes(x = x, y = y, z = mean.ses), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.ses), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.5, size = 0.25, aes(x = x, y = y, z = sd.ses), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = freq.pval), data = ann) +
    scale_fill_distiller(name = "Freq. of\np < .05", palette = "Blues", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 25, size = 0.25, aes(x = x, y = y, z = freq.pval), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.faith), data = ann) +
     scale_fill_viridis(name = "Mean annual\nFaith's index") +
     geom_contour(colour = "grey30", binwidth = 0.25, size = 0.25, aes(x = x, y = y, z = mean.faith), data = ann) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

panel.maps.ann <- ggarrange(p4,p1,p3, align = 'hv', ncol = 1, nrow = 3)    


### Monthly 
mon <- data.frame(tab %>% group_by(month,cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.ses = mean(pd.obs.z, na.rm = T), sd.ses = sd(pd.obs.z, na.rm = T)
        ) # eo summarize
)
# dim(mon); summary(mon)
# Re-order months correctly
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
mon$month <- factor(mon$month, months)

p6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.ses), data = mon) +
    scale_fill_gradient2(name = "Monthly SES", low = "#3288bd", mid = "white", high = "#b2182b") +
    geom_contour(colour = "grey30", binwidth = 2, size = 0.25, aes(x = x, y = y, z = mean.ses), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)


### SDM specific
sdms <- data.frame(tab %>% group_by(SDM,cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.ses = mean(pd.obs.z, na.rm = T), sd.ses = sd(pd.obs.z, na.rm = T)
        ) # eo summarize
)
# Re-order months correctly
sdms$SDM <- factor(sdms$SDM, c("GLM","GAM","ANN"))

p7 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.ses), data = sdms) +
    scale_fill_gradient2(name = "Mean annual\nSES", low = "#3288bd", mid = "white", high = "#b2182b") +
    geom_contour(colour = "grey30", binwidth = 1, size = 0.25, aes(x = x, y = y, z = mean.ses), data = sdms) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(SDM), ncol = 1)


### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------

##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 04/07/23: R script to map 'funrar' indices (functional distinctiveness and species scarcity) computed based on a Gower distance matrix and the community tables made during the MSc thesis of Jonas Wydler (data from Benedetti et al., 2023 - JBIO) © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - load the tables containing the funrar indices values 
# - compute and map the mean (ensemble) annual indices, the associated stdev
# - compute and map the monthly averages, and the SDM-specific predictions
# - same but with P/A data

### !!! BEWARE THIS PART OF THE ANALYSES WAS NOT KEPT FOR THE MANUSCRIPT !!!

### Latest update: 04/07/23

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

setwd("/net/kryo/work/fabioben/BIOCeans5D/data") # working dir

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 1) Load the tables of Faith's indices and compute ensemble averages and uncertainties around it
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/funrar/") # dir() # should be of length 36
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
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) # length(unique(tab$cell_id))

### Compute mean annual indices
ann <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.dist = mean(distinct, na.rm = T), sd.dist = sd(distinct, na.rm = T),
            mean.scar = mean(scarcity, na.rm = T), sd.scar = sd(scarcity, na.rm = T))
)
# dim(ann); summary(ann)

### Draw maps
p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.dist), data = ann) + scale_fill_viridis(name = "Mean annual\ndistinctiveness") +
    geom_contour(colour = "grey30", binwidth = .005, size = 0.25, aes(x = x, y = y, z = mean.dist), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.dist), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.001, size = 0.25, aes(x = x, y = y, z = sd.dist), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.scar), data = ann) + scale_fill_viridis(name = "Mean annual\nscarcity") +
    geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = mean.scar), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.scar), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = sd.scar), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

panel.maps.ann <- ggarrange(p3,p4,p1,p2, align = 'hv', ncol = 2, nrow = 2)    


### Same as above but monthly
mon <- data.frame(tab %>% group_by(month,cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.dist = mean(distinct, na.rm = T), sd.dist = sd(distinct, na.rm = T),
            mean.scar = mean(scarcity, na.rm = T), sd.scar = sd(scarcity, na.rm = T))
)
# dim(mon); summary(mon)
# Re-order months correctly
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
mon$month <- factor(mon$month, months)

p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.dist), data = mon) + scale_fill_viridis(name = "Monthly distinctiveness") +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = mean.dist), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.scar), data = mon) + scale_fill_viridis(name = "Monthly scarcity") +
    geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = mean.scar), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)


### Same as above but per SDM
sdms <- data.frame(tab %>% group_by(SDM,cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.dist = mean(distinct, na.rm = T), sd.dist = sd(distinct, na.rm = T),
            mean.scar = mean(scarcity, na.rm = T), sd.scar = sd(scarcity, na.rm = T))
)
# dim(sdms); summary(sdms)
sdms$SDM <- factor(sdms$SDM, c("GLM","GAM","ANN"))

p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.dist), data = sdms) + scale_fill_viridis(name = "Mean annual\ndistinctiveness") +
    geom_contour(colour = "grey30", binwidth = 0.005, size = 0.25, aes(x = x, y = y, z = mean.dist), data = sdms) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(SDM), nrow = 3, ncol = 1)

p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = ), data = sdms) + scale_fill_viridis(name = "Mean annual\nscarcity") +
    geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = ), data = sdms) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(SDM), nrow = 3, ncol = 1)


### ------------------------------------------------------------------

### 2) Examine mean (and stdev) of species' functional uniqueness (not spatial)
files <- dir()[grep("species_funrar",dir())]
# f <- files[1]
res <- mclapply(files, function(f) {
            d <- get(load(f))
            # extract month and sdm from filename 
            filename <- str_replace_all(f,".Rdata","")
            # unlist(strsplit(x = filename, split = "_", fixed = T))
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[7]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[8]
            return(d)
        }, mc.cores = 10
) # eo mclapply
tab <- bind_rows(res)
rm(res); gc()
# dim(tab); head(tab); summary(tab)

months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
tab$month <- factor(tab$month, months)
tab$SDM <- factor(tab$SDM, c("GLM","GAM","ANN"))

### Quickly check distribution of Ui per months and SDM (shouldn't change too much)
#ggplot(data = tab, aes(x = factor(SDM), y = log10(Ui))) + geom_boxplot(fill = "gray") + xlab("SDM") + ylab("Functional uniqueness") + theme_bw() # good
# Per months
#ggplot(data = tab, aes(x = factor(month), y = log10(Ui))) + geom_boxplot(fill = "gray") + xlab("Month") + ylab("Functional uniqueness") + theme_bw() # good
# good as well

# Compute mean Ui and Ri per species and examine ranks and cor
avg <- data.frame( dat %>% group_by(species) %>% summarise(Ui = mean(Ui, na.rm = T), Ri = mean(Ri, na.rm = T)) )
# summary(avg)

### NOTE: not very informative based on HSI. Use the binary outputs instead...

### ------------------------------------------------------------------

### 3) Same as 1°, but based on the Distinctiveness values computed based on PA matrices
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/funrar/PA") # dir() # should be of length 36
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
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) # length(unique(tab$cell_id))

### Compute mean annual indices
ann <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.dist = mean(distinct, na.rm = T), sd.dist = sd(distinct, na.rm = T)
    ) # eo summarise
)
# dim(ann); summary(ann)

### Draw maps
p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.dist), data = ann) + scale_fill_viridis(name = "Mean annual\ndistinctiveness") +
    geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = mean.dist), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.dist), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = sd.dist), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

panel.maps.ann <- ggarrange(p1,p2, align = 'hv', ncol = 2, nrow = 1)    
### Way clearer than with HSI


### Same as above but monthly
mon <- data.frame(tab %>% group_by(month,cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.dist = mean(distinct, na.rm = T), sd.dist = sd(distinct, na.rm = T)
    )
)
# dim(mon); summary(mon)
# Re-order months correctly
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
mon$month <- factor(mon$month, months)

p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.dist), data = mon) + scale_fill_viridis(name = "Monthly distinctiveness") +
    geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = mean.dist), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)


### Same as above but per SDM
sdms <- data.frame(tab %>% group_by(SDM,cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.dist = mean(distinct, na.rm = T), sd.dist = sd(distinct, na.rm = T),
    )
)
# dim(sdms); summary(sdms)
sdms$SDM <- factor(sdms$SDM, c("GLM","GAM","ANN"))

p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean.dist), data = sdms) + scale_fill_viridis(name = "Mean annual\ndistinctiveness") +
    geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = mean.dist), data = sdms) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(SDM), nrow = 3, ncol = 1)



### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
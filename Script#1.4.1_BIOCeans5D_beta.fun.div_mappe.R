
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 12/07/23: R script to map functional beta diversity indices (Jac/Jtu/Jne) computed based on a Gower distance matrix and the community tables made during the MSc thesis of Jonas Wydler (data from Benedetti et al., 2023 - JBIO) © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - load the tables containing the functional beta diversity indices (Jac/Jtu/Jne)
# - compute and map the mean (ensemble) annual indices, the associated stdev
# - compute and map the monthly averages, and the SDM-specific predictions

### Latest update: 13/07/23

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
setwd("/net/kryo/work/fabioben/GODLY/data/fd_indices/beta.div/") # dir() # should be of length 36
files <- dir()[grep("baseline",dir())]
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

### Examine distributions
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
tab$month <- factor(tab$month,months)
ggplot(data = tab, aes(x = factor(SDM), y = beta.jac)) + geom_boxplot() + xlab("SDM") + ylab("Jaccard") + facet_wrap(.~factor(month), ncol = 4)
ggplot(data = tab, aes(x = factor(SDM), y = beta.jtu)) + geom_boxplot() + xlab("SDM") + ylab("Turnover") + facet_wrap(.~factor(month), ncol = 4)
ggplot(data = tab, aes(x = factor(SDM), y = beta.jne)) + geom_boxplot() + xlab("SDM") + ylab("Nestedness") + facet_wrap(.~factor(month), ncol = 4)
# beta ratio
ggplot(data = tab, aes(x = factor(SDM), y = beta.jtu/beta.jac)) + geom_boxplot() + xlab("SDM") + ylab("ßratio") + facet_wrap(.~factor(month), ncol = 4)

tab$ratio <- tab$beta.jtu/tab$beta.jac

### Compute mean annual indices
ann <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
        jac = mean(beta.jac, na.rm = T), jac.std = sd(beta.jac, na.rm = T),
        jtu = mean(beta.jtu, na.rm = T), jtu.std = sd(beta.jtu, na.rm = T),
        jac = mean(beta.jne, na.rm = T), jne.std = sd(beta.jne, na.rm = T),
        beta.ratio = mean(ratio, na.rm = T), beta.ratio.std = sd(ratio, na.rm = T)
    ) 
)
# dim(ann); summary(ann)


### Draw maps
p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac), data = ann) + scale_fill_viridis(name = "Mean annual\nJaccard") +
    geom_contour(colour = "grey30", binwidth = .05, size = 0.25, aes(x = x, y = y, z = jac), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = jac.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.02, size = 0.25, aes(x = x, y = y, z = jac.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu), data = ann) + scale_fill_viridis(name = "Mean annual\nTurnover") +
    geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = jtu), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = jtu.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.02, size = 0.25, aes(x = x, y = y, z = jtu.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
p5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = beta.ratio), data = ann) + scale_fill_viridis(name = "Mean annual\nßratio") +
    geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = beta.ratio), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = beta.ratio.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = beta.ratio.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

panel.maps.ann <- ggarrange(p1,p2,p3,p4,p5,p6, align = 'hv', ncol = 2, nrow = 3)    

# To save individual maps 
ggsave(plot = p5, filename = "map_mean_ann_Bratio.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = p3, filename = "map_mean_ann_Bjtu.jpg", dpi = 300, height = 4, width = 7)


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
        jac = mean(beta.jac, na.rm = T), jac.std = sd(beta.jac, na.rm = T),
        jtu = mean(beta.jtu, na.rm = T), jtu.std = sd(beta.jtu, na.rm = T),
        jac = mean(beta.jne, na.rm = T), jne.std = sd(beta.jne, na.rm = T),
        beta.ratio = mean(ratio, na.rm = T), beta.ratio.std = sd(ratio, na.rm = T)
    )    
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


### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
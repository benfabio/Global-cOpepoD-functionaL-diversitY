
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 07/08/23: R script to map 'dbFD()' indices (FEve,FDiv,FDis,RaoQ,FRic) computed based on a Gower distance matrix and the community tables made during the MSc thesis of Jonas Wydler (data from Benedetti et al., 2023 - JBIO) © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - load the tables containing the dbFD() indices values (based on HSI+Gower dist matrix+PCoA)
# - compute and map the mean (ensemble) annual indices, the associated stdev
# - compute and map the monthly averages, and the SDM-specific predictions
# - same but with P/A data

### Latest update: 24/08/23

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
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid/") # dir() # should be of length 36
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

# Add cell id
tab$cell_id <- factor(paste(tab$x,tab$y,sep = "_")) # length(unique(tab$cell_id))

### To check ranges of values per SDM/months
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
tab$month <- factor(tab$month,months)
ggplot(data = tab, aes(x = factor(SDM), y = FRic)) + geom_boxplot() + xlab("SDM") + ylab("FRic") + facet_wrap(.~factor(month), ncol = 4)
ggplot(data = tab, aes(x = factor(SDM), y = FEve)) + geom_boxplot() + xlab("SDM") + ylab("FEve") + facet_wrap(.~factor(month), ncol = 4)
ggplot(data = tab, aes(x = factor(SDM), y = FDis)) + geom_boxplot() + xlab("SDM") + ylab("FDis") + facet_wrap(.~factor(month), ncol = 4)
ggplot(data = tab, aes(x = factor(SDM), y = RaoQ)) + geom_boxplot() + xlab("SDM") + ylab("RaoQ") + facet_wrap(.~factor(month), ncol = 4)
ggplot(data = tab, aes(x = factor(SDM), y = FDiv)) + geom_boxplot() + xlab("SDM") + ylab("FDiv") + facet_wrap(.~factor(month), ncol = 4)

### STRONGEST SDM-DRIVEN VARIABILITY FOR: FEve (higher for ANN > GLM > GAM) and FRic (higher for ANN > GAM > GLM)
#summary(tab[tab$SDM == "ANN",])
# ggplot() + geom_raster(aes(x = x, y = y, fill = FRic), data = tab[tab$SDM == "ANN",]) + scale_fill_viridis(name = "FRic") +
#     geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = FRic), data = tab[tab$SDM == "ANN",]) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) + facet_wrap(.~ month)

### Compute mean annual indices
ann <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            FRic.avg = mean(FRic, na.rm = T), FRic.std = sd(FRic, na.rm = T),
            FEve.avg = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
            FDis.avg = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
            RaoQ.avg = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
            FDiv.avg = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T) )
)
# dim(ann); summary(ann)

### Draw maps: FRic (should not be interesting), FEve, FDiv, FDis, Rao'Q
p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nFRic") +
    geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = FRic.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = FRic.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FEve.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nFEve") +
    geom_contour(colour = "grey30", binwidth = .05, size = 0.25, aes(x = x, y = y, z = FEve.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FEve.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = FEve.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDis.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nFDis") +
    geom_contour(colour = "grey30", binwidth = 0.005, size = 0.25, aes(x = x, y = y, z = FDis.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDis.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.001, size = 0.25, aes(x = x, y = y, z = FDis.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p7 <- ggplot() + geom_raster(aes(x = x, y = y, fill = RaoQ.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nRao's Q") +
    geom_contour(colour = "grey30", binwidth = 0.001, size = 0.25, aes(x = x, y = y, z = RaoQ.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p8 <- ggplot() + geom_raster(aes(x = x, y = y, fill = RaoQ.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = .001, size = 0.25, aes(x = x, y = y, z = RaoQ.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p9 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDiv.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nFDiv") +
    geom_contour(colour = "grey30", binwidth = 0.005, size = 0.25, aes(x = x, y = y, z = FDiv.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p10 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDiv.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = .005, size = 0.25, aes(x = x, y = y, z = FDiv.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

panel.maps.ann <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, align = 'hv', ncol = 2, nrow = 5)    

# To save individual maps 
ggsave(plot = p3, filename = "map_mean_ann_FEve.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = p5, filename = "map_mean_ann_FDis.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = p7, filename = "map_mean_ann_RaoQ.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = p9, filename = "map_mean_ann_FDiv.jpg", dpi = 300, height = 4, width = 7)


### Compute mean monthly indices
mon <- data.frame(tab %>% group_by(cell_id,month) %>% 
        summarize(x = unique(x), y = unique(y),
            FRic.avg = mean(FRic, na.rm = T), FRic.std = sd(FRic, na.rm = T),
            FEve.avg = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
            FDis.avg = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
            RaoQ.avg = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
            FDiv.avg = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T) )
)
# dim(mon); summary(mon)
# Re-order months correctly
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
mon$month <- factor(mon$month,months)

p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.avg), data = mon) + scale_fill_viridis(name = "Monthly FRic") +
    geom_contour(colour = "grey30", binwidth = .05, size = 0.25, aes(x = x, y = y, z = FRic.avg), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FEve.avg), data = mon) + scale_fill_viridis(name = "Monthly FEve") +
    geom_contour(colour = "grey30", binwidth = .1, size = 0.25, aes(x = x, y = y, z = FEve.avg), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)

p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDis.avg), data = mon) + scale_fill_viridis(name = "Monthly FDis") +
    geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = FDis.avg), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)

p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = RaoQ.avg), data = mon) + scale_fill_viridis(name = "Monthly RaoQ") +
    geom_contour(colour = "grey30", binwidth = 0.005, size = 0.25, aes(x = x, y = y, z = RaoQ.avg), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)
    
p5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDiv.avg), data = mon) + scale_fill_viridis(name = "Monthly FDiv") +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = FDiv.avg), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)


### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 2) Same as above, but for indices based on PA data
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/PA_Gawdis_PCoA_Euclid/") # dir() # should be of length 36
files <- dir()[grep("FDindices_baseline",dir())]
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
tab$cell_id <- factor(paste(tab$x,tab$y,sep = "_")) # length(unique(tab$cell_id))

### Compute mean annual indices
ann <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            FRic.avg = mean(FRic, na.rm = T), FRic.std = sd(FRic, na.rm = T),
            FEve.avg = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
            FDis.avg = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
            RaoQ.avg = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
            FDiv.avg = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T) )
)
# dim(ann); summary(ann)

p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nFRic") +
    geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = FRic.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = FRic.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FEve.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nFEve") +
    geom_contour(colour = "grey30", binwidth = .05, size = 0.25, aes(x = x, y = y, z = FEve.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FEve.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = FEve.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDis.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nFDis") +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = FDis.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDis.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.005, size = 0.25, aes(x = x, y = y, z = FDis.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p7 <- ggplot() + geom_raster(aes(x = x, y = y, fill = RaoQ.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nRao's Q") +
    geom_contour(colour = "grey30", binwidth = 0.02, size = 0.25, aes(x = x, y = y, z = RaoQ.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p8 <- ggplot() + geom_raster(aes(x = x, y = y, fill = RaoQ.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = .02, size = 0.25, aes(x = x, y = y, z = RaoQ.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p9 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDiv.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nFDiv") +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = FDiv.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p10 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDiv.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = FDiv.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

panel.maps.ann <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, align = 'hv', ncol = 2, nrow = 5)  


### Monthly variability/ SDM variability

# Compute mean monthly indices
mon <- data.frame(tab %>% group_by(cell_id,month) %>% 
        summarize(x = unique(x), y = unique(y),
            FRic.avg = mean(FRic, na.rm = T), FRic.std = sd(FRic, na.rm = T),
            FEve.avg = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
            FDis.avg = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
            RaoQ.avg = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
            FDiv.avg = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T) )
)
# Re-order months correctly
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
mon$month <- factor(mon$month,months)

ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.avg), data = mon) + scale_fill_viridis(name = "Monthly FRic") +
    geom_contour(colour = "grey30", binwidth = .05, size = 0.25, aes(x = x, y = y, z = FRic.avg), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)
### --> ISSUES WIH JAN/SEPT/DEC; Need to check input

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FEve.avg), data = mon) + scale_fill_viridis(name = "Monthly FEve") +
    geom_contour(colour = "grey30", binwidth = .1, size = 0.25, aes(x = x, y = y, z = FEve.avg), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)
### OK, checks out

p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDis.avg), data = mon) + scale_fill_viridis(name = "Monthly FDis") +
    geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = FDis.avg), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)
### OK, checks out

p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = RaoQ.avg), data = mon) + scale_fill_viridis(name = "Monthly RaoQ") +
    geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = RaoQ.avg), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)
### OK, checks out
    
p5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDiv.avg), data = mon) + scale_fill_viridis(name = "Monthly FDiv") +
    geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = FDiv.avg), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)
### OK, checks out

# Same as above but per SDM
sdms <- data.frame(tab %>% group_by(SDM,cell_id) %>% 
            summarize(x = unique(x), y = unique(y),
                FRic.avg = mean(FRic, na.rm = T), FRic.std = sd(FRic, na.rm = T),
                FEve.avg = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
                FDis.avg = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
                RaoQ.avg = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
                FDiv.avg = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T)
        )
)
# dim(sdms); summary(sdms)
sdms$SDM <- factor(sdms$SDM, c("GLM","GAM","ANN"))

# Plot iteratively for now
ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.avg), data = sdms) + scale_fill_viridis(name = "Mean annual\nFRic") +
    geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = FRic.avg), data = sdms) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(SDM), nrow = 3, ncol = 1)
### ISSUE WITH GAM AND ANN FOR FRic

### FEve/FDis/RaoQ: check out and smaller SDM variability
ggplot() + geom_raster(aes(x = x, y = y, fill = FDiv.avg), data = sdms) + scale_fill_viridis(name = "Mean annual\nFDiv") +
    geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = FDiv.avg), data = sdms) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(SDM), nrow = 3, ncol = 1)
    
### Large SDM uncertainty with FDiv!


### Identify the problem with FRic: is all months for GAM and ANN?
tab$month <- factor(tab$month,months)

# Ranges per SDMs
#summary(tab[tab$SDM == "GLM","FRic"])
#summary(tab[tab$SDM == "GAM","FRic"])
#summary(tab[tab$SDM == "ANN","FRic"])

ggplot() + geom_raster(aes(x = x, y = y, fill = FRic), data = tab[tab$SDM == "GLM",]) + scale_fill_viridis(name = "FRic") +
    geom_contour(colour = "grey30", binwidth = .1, size = 0.25, aes(x = x, y = y, z = FRic), data = tab[tab$SDM == "GLM",]) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)
    
ggplot() + geom_raster(aes(x = x, y = y, fill = FRic), data = tab[tab$SDM == "GAM",]) + scale_fill_viridis(name = "FRic") +
    geom_contour(colour = "grey30", binwidth = .1, size = 0.25, aes(x = x, y = y, z = FRic), data = tab[tab$SDM == "GAM",]) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)

ggplot() + geom_raster(aes(x = x, y = y, fill = FRic), data = tab[tab$SDM == "ANN",]) + scale_fill_viridis(name = "FRic") +
    geom_contour(colour = "grey30", binwidth = .1, size = 0.25, aes(x = x, y = y, z = FRic), data = tab[tab$SDM == "ANN",]) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)

# Issue is that GLM-based FRic ranges from 0 to 0.43 while GAM-based and ANN-based range from 0 to 0.39...like if a factor 100 was missing

#summary(tab[tab$SDM == "GLM" & tab$month == "jan","FRic"])
#summary(tab[tab$SDM == "GLM" & tab$month == "apr","FRic"])
#summary(tab[tab$SDM == "GLM" & tab$month == "dec","FRic"]) # issue with december

# Examine distrib of FRic per SDMx months with boxplots
ggplot(data = tab, aes(x = factor(SDM), y = FRic)) + geom_boxplot() + xlab("SDM") + ylab("FRic") + facet_wrap(.~factor(month), ncol = 4)

# We know that FRic is biased by SR: higher SR --> higher FRic, so much higher FRic from GLM could be due to higher predicted SR! 
ggplot(data = tab, aes(x = factor(SDM), y = SR)) + geom_boxplot() + xlab("SDM") + ylab("SR") + facet_wrap(.~factor(month), ncol = 4)
### --> GLM does show stronges ranges and higher SR. Could explain.

### This one makes sense.
ggplot() + geom_raster(aes(x = x, y = y, fill = FRic), data = tab[tab$month == "dec",]) + scale_fill_viridis(name = "FRic") +
    geom_contour(colour = "grey30", binwidth = .1, size = 0.25, aes(x = x, y = y, z = FRic), data = tab[tab$month == "dec",]) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ SDM, ncol = 1)
    
    
### Pick a month where strong differences in FRic are apparent: 'apr' and compare comm tables. 
## Then do the same with 'dec' since it showed comparable FRic levels
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/community_tables_05_01_21/contemp")
glm.comm <- read.table("table_mon_composition_baseline_GLM_dec.txt")
gam.comm <- read.table("table_mon_composition_baseline_GAM_dec.txt")
ann.comm <- read.table("table_mon_composition_baseline_ANN_dec.txt")
# dim(glm.comm); dim(gam.comm); dim(ann.comm) # Same dimensions, check
# setdiff(colnames(glm.comm), colnames(gam.comm))
# setdiff(colnames(gam.comm), colnames(ann.comm))

# Ranges of HSI? 
# summary(glm.comm[,c(5:10)])
# summary(gam.comm[,c(5:10)])
# summary(ann.comm[,c(5:10)])
# Seems OK

# Load the HSI cutoffs and examine the 1/0 data (SR and comp.)
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/sdm/evaluation_scores_05_01_21")
score.files <- dir()#; files
scores <- lapply(score.files, function(f) {d<-get(load(f)); return(d)})
tab.scores <- bind_rows(scores)
rm(score.files); gc()
tab.scores$SDM <- NA
tab.scores[grepl("GLM",rownames(tab.scores)),"SDM"] <- "GLM"
tab.scores[grepl("ANN",rownames(tab.scores)),"SDM"] <- "ANN"
tab.scores[grepl("GAM",rownames(tab.scores)),"SDM"] <- "GAM"
cutoffs <- data.frame( tab.scores %>% group_by(species,SDM) %>% summarise(cutoff = mean(Cutoff_TSS, na.rm = T)/1000) )
#summary(cutoffs) # seems ok

library("biomod2")
# Convert to 1/0 for each comm
# For computing Faith's index, need to use 1/0 community data (ensemble of thresholds)
glm.cutoff <- cutoffs[cutoffs$SDM == "GLM",]
gam.cutoff <- cutoffs[cutoffs$SDM == "GAM",]
ann.cutoff <- cutoffs[cutoffs$SDM == "ANN",]
    
# Use those cutoffs to derive a PA (1/0) comm table
glm.comm_PA <- glm.comm
gam.comm_PA <- gam.comm
ann.comm_PA <- ann.comm

spp <- colnames(glm.comm_PA)[c(4:length(glm.comm_PA))] #; spp
# s <- "Euchaeta_spinosa"
spp <- spp[!spp == "Euchaeta_spinosa"]
spp <- spp[!spp == "Sapphirina_ovatolanceolata"]
spp <- spp[!spp == "Scaphocalanus_echinatus"]
spp <- spp[!spp == "Scolecithricella_abyssalis"]
spp <- spp[!spp == "Spinocalanus_magnus"]

for(s in spp) {
    
    message(s)
    
    t <- glm.cutoff[glm.cutoff$species == s,"cutoff"]
    glm.comm_PA[,s] <- bm_BinaryTransformation(data = glm.comm_PA[,s], threshold = t)
    
    t <- gam.cutoff[gam.cutoff$species == s,"cutoff"]
    gam.comm_PA[,s] <- bm_BinaryTransformation(data = gam.comm_PA[,s], threshold = t)
    
    t <- ann.cutoff[ann.cutoff$species == s,"cutoff"]
    ann.comm_PA[,s] <- bm_BinaryTransformation(data = ann.comm_PA[,s], threshold = t)
    
} # eo for loop - c in commons

# Check if there are species (i.e., columns) with only 0; retrive name of given species and remove from 'comm_fdiv_PA'
### Compare SR per model
sr.glm <- rowSums(glm.comm_PA[,c(4:length(glm.comm_PA))], na.rm = T) 
# summary(sr.glm)
sr.gam <- rowSums(gam.comm_PA[,c(4:length(gam.comm_PA))], na.rm = T) 
# summary(sr.gam)
sr.ann <- rowSums(ann.comm_PA[,c(4:length(ann.comm_PA))], na.rm = T) 
# summary(sr.ann)
### comparable ranges...map

ann.comm_PA$SR <- sr.ann
gam.comm_PA$SR <- sr.gam
glm.comm_PA$SR <- sr.glm

p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = SR), data = glm.comm_PA) + scale_fill_viridis(name = "Monthly SR\n(GLM)", limits = c(0,175)) +
    geom_contour(colour = "grey30", binwidth = 20, size = 0.25, aes(x = x, y = y, z = SR), data = glm.comm_PA) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) 

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = SR), data = gam.comm_PA) + scale_fill_viridis(name = "Monthly SR\n(GAM)", limits = c(0,175)) +
    geom_contour(colour = "grey30", binwidth = 20, size = 0.25, aes(x = x, y = y, z = SR), data = gam.comm_PA) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) 

p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = SR), data = ann.comm_PA) + scale_fill_viridis(name = "Monthly SR\n(ANN)", limits = c(0,175)) +
    geom_contour(colour = "grey30", binwidth = 20, size = 0.25, aes(x = x, y = y, z = SR), data = ann.comm_PA) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) 
    
ggarrange(p1,p2,p3, align = 'hv', ncol = 1, nrow = 3)  

### Examine corr coefficients
#cor(ann.comm_PA$SR, gam.comm_PA$SR, method = "spearman") # basically the same variable
#cor(glm.comm_PA$SR, gam.comm_PA$SR, method = "spearman") # cor is 0.89 in apr so very high too; 0.88 in december
# Plot
# plot(x = glm.comm_PA$SR, y = gam.comm_PA$SR)
# abline(0,1,col='red')

### SR differences alone do not explain wild differences in FRic...Examine values of standardized FRic

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 10/08/23: Examine distribution of standardized FRic data
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/PA_Gawdis_PCoA_Euclid/Stand.Fric") # dir() # should be of length 36
files <- dir()[grep("FDindices_baseline",dir())]
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
tab$cell_id <- factor(paste(tab$x,tab$y,sep = "_")) # length(unique(tab$cell_id))

### Compute mean annual indices
ann <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            FRic.avg = mean(FRic, na.rm = T), FRic.std = sd(FRic, na.rm = T),
            FEve.avg = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
            FDis.avg = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
            RaoQ.avg = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
            FDiv.avg = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T) )
)
# dim(ann); summary(ann)

p1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nFRic") +
    geom_contour(colour = "grey30", binwidth = .05, size = 0.25, aes(x = x, y = y, z = FRic.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.05, size = 0.25, aes(x = x, y = y, z = FRic.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FEve.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nFEve") +
    geom_contour(colour = "grey30", binwidth = .05, size = 0.25, aes(x = x, y = y, z = FEve.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FEve.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = FEve.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDis.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nFDis") +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = FDis.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDis.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = 0.005, size = 0.25, aes(x = x, y = y, z = FDis.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p7 <- ggplot() + geom_raster(aes(x = x, y = y, fill = RaoQ.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nRao's Q") +
    geom_contour(colour = "grey30", binwidth = 0.02, size = 0.25, aes(x = x, y = y, z = RaoQ.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p8 <- ggplot() + geom_raster(aes(x = x, y = y, fill = RaoQ.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = .02, size = 0.25, aes(x = x, y = y, z = RaoQ.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p9 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDiv.avg), data = ann) + scale_fill_viridis(name = "Mean annual\nFDiv") +
    geom_contour(colour = "grey30", binwidth = 0.01, size = 0.25, aes(x = x, y = y, z = FDiv.avg), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

p10 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDiv.std), data = ann) +
    scale_fill_distiller(name = "Uncertainty", palette = "YlOrRd", direction = 1) +
    geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = FDiv.std), data = ann) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

panel.maps.ann <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, align = 'hv', ncol = 2, nrow = 5)  

### Monthly variability/ SDM variability
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
tab$month <- factor(tab$month,months)

# Draw boxplots
#ggplot(data = tab, aes(x = factor(SDM), y = SR)) + geom_boxplot() + xlab("SDM") + ylab("SR") + facet_wrap(.~factor(month), ncol = 4)
# And the other way around
# ggplot(data = tab, aes(x = factor(month), y = SR)) + geom_boxplot() + xlab("SDM") + ylab("SR") + facet_wrap(.~factor(SDM), ncol = 3)

#ggplot(data = tab, aes(x = factor(SDM), y = FRic)) + geom_boxplot() + xlab("SDM") + ylab("FRic (standardized)") + facet_wrap(.~factor(month), ncol = 4)
#ggplot(data = tab, aes(x = factor(month), y = FRic)) + geom_boxplot() + xlab("SDM") + ylab("FRic (standardized)") + facet_wrap(.~factor(SDM), ncol = 3)
#ggplot(data = tab, aes(x = factor(SDM), y = FEve)) + geom_boxplot() + xlab("SDM") + ylab("FEve") + facet_wrap(.~factor(month), ncol = 4)
#ggplot(data = tab, aes(x = factor(SDM), y = FDis)) + geom_boxplot() + xlab("SDM") + ylab("FDis") + facet_wrap(.~factor(month), ncol = 4)
#ggplot(data = tab, aes(x = factor(SDM), y = FDiv)) + geom_boxplot() + xlab("SDM") + ylab("FDiv") + facet_wrap(.~factor(month), ncol = 4)


# Compute mean monthly indices
mon <- data.frame(tab %>% group_by(cell_id,month) %>% 
        summarize(x = unique(x), y = unique(y),
            SR.avg = mean(SR, na.rm = T), SR.std = sd(SR, na.rm = T), 
            FRic.avg = mean(FRic, na.rm = T), FRic.std = sd(FRic, na.rm = T),
            FEve.avg = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
            FDis.avg = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
            RaoQ.avg = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
            FDiv.avg = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T)
    )
)
# Re-order months correctly
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
mon$month <- factor(mon$month,months)

ggplot() + geom_raster(aes(x = x, y = y, fill = SR.avg), data = mon) + scale_fill_viridis(name = "Monthly SR") +
    geom_contour(colour = "grey30", binwidth = 20, size = 0.25, aes(x = x, y = y, z = SR.avg), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)
    
ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.avg), data = mon) + scale_fill_viridis(name = "Monthly FRic") +
    geom_contour(colour = "grey30", binwidth = .1, size = 0.25, aes(x = x, y = y, z = FRic.avg), data = mon) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ month)


# SDM specific estimates
sdms <- data.frame(tab %>% group_by(SDM,cell_id) %>% 
            summarize(x = unique(x), y = unique(y),
                SR.avg = mean(SR, na.rm = T), SR.std = sd(SR, na.rm = T), 
                FRic.avg = mean(FRic, na.rm = T), FRic.std = sd(FRic, na.rm = T),
                FEve.avg = mean(FEve, na.rm = T), FEve.std = sd(FEve, na.rm = T),
                FDis.avg = mean(FDis, na.rm = T), FDis.std = sd(FDis, na.rm = T),
                RaoQ.avg = mean(RaoQ, na.rm = T), RaoQ.std = sd(RaoQ, na.rm = T),
                FDiv.avg = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T)
        )
)
# dim(sdms); summary(sdms)
sdms$SDM <- factor(sdms$SDM, c("GLM","GAM","ANN"))

# Draw maps
ggplot() + geom_raster(aes(x = x, y = y, fill = SR.avg), data = sdms) + scale_fill_viridis(name = "SR") +
    geom_contour(colour = "grey30", binwidth = 20, size = 0.25, aes(x = x, y = y, z = SR.avg), data = sdms) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(SDM), nrow = 3, ncol = 1)
    
ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.avg), data = sdms) + scale_fill_viridis(name = "FRic (standardized)") +
    geom_contour(colour = "grey30", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = FRic.avg), data = sdms) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(SDM), nrow = 3, ncol = 1)

# And iteratively for the other indices
ggplot() + geom_raster(aes(x = x, y = y, fill = FEve.avg), data = sdms) + scale_fill_viridis(name = "FEve") +
    geom_contour(colour = "grey30", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = FEve.avg), data = sdms) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(SDM), nrow = 3, ncol = 1)
    
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------

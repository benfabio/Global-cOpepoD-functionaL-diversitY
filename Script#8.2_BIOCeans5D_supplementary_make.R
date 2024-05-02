
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 24/11/23: R script to make the supplementary material for the manuscript in preparation for Global Change Biology © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Latest update: 18/04/24 (calculating relative standard deviation for intra-annual variability (aka seasonality) and SDM-based variability)

### ------------------------------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("viridis")
library("parallel")
library("ggpubr")
library("ggthemes")
library("pals") 
library("marmap")
library("FactoMineR")
library("ggpmisc")
library("raster")

world <- map_data("world") 
world2 <- map_data("world2")

setwd("/net/kryo/work/fabioben/BIOCeans5D/data")

### ------------------------------------------------------------------------------------------------------------------------------------------------------

fd <- get(load("table_mean_ann_FD_indices_baseline+BCP+biom_22.11.23.RData"))

### 1) Faith ~ FRic; FDis ~ Rao's Q
formula1 <- y~x # summary(lm)
# p1 <- ggplot(fd[fd$FRic > 0.05,], aes(x = Faith, y = FRic, colour = abs(y))) +
#   geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
#   scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
#   stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "bottom", label.x = "right", size = 3) +
#   xlab("Faith index") + ylab("FRic") + theme_bw()

p2 <- ggplot(fd, aes(x = FDis, y = RaoQ.scaled, colour = abs(y))) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "left", size = 3) + 
  xlab("FDis") + ylab("Rao's Q (scaled)") + theme_bw()

setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
ggsave(plot = p2, filename = "Fig.X_FdisxRaoQ_01.12.23.jpg", dpi = 300, width = 5, height = 4)

#SM1 <- ggarrange(p1,p2, align = 'hv', ncol = 2, nrow = 1, labels = letters[1:2], common.legend = T)
#ggsave(plot = SM1, filename = "Fig.X_panel_FDisxRaoQ+FaithxFRic_27.11.23.jpg", dpi = 300, width = 10, height = 6)

npp <- ggplot(fd, aes(x = log10(NPP), y = log10(NPPv2), colour = abs(y))) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "left", size = 3) + 
  xlab("NPP (standard VGPM) - log10(mgC.m-2.d-1)") + ylab("NPP (deVries & Weber, 2017) - log10(mgC.m-2.d-1)") + theme_bw()

ggsave(plot = npp, filename = "Fig.X_NPPv1xNPPv2_27.11.23.jpg", dpi = 300, width = 5, height = 4)

### ------------------------------------------------------------------

### 2) Maps of EF proxies (POC, MESOZOOPL, E RATIO etc.)
vars2map <- colnames(fd)[c(21:25,29,30,32:35,39:40)]; vars2map # 13 vars --> make two 7-map manels (2 colsx4 rows)
# v <- "CHL"

m1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = MESOZOO), data = fd) +
    scale_fill_gradientn(name = "Mesozooplankton\nbiomass\n(mgC.m-3)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = 1, size = .4, aes(x = x, y = y, z = MESOZOO), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(CHL)), data = fd) +
    scale_fill_gradientn(name = "Chlorophyll-a\nconcentration\nlog10(mg.m-3)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .25, size = .4, aes(x = x, y = y, z = log10(CHL)), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(DIATO)), data = fd) +
    scale_fill_gradientn(name = "Diatom biomass\nlog10(mg.m-3)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .25, size = .4, aes(x = x, y = y, z = log10(DIATO)), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(DINO)), data = fd) +
    scale_fill_gradientn(name = "Dinoflagellate\nbiomass\nlog10(mg.m-3)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .25, size = .4, aes(x = x, y = y, z = log10(DINO)), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(GREEN)), data = fd) +
    scale_fill_gradientn(name = "Green algae\nbiomass\nlog10(mg.m-3)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .25, size = .4, aes(x = x, y = y, z = log10(GREEN)), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(HAPTO)), data = fd) +
    scale_fill_gradientn(name = "Haptophyte\nbiomass\nlog10(mg.m-3)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .25, size = .4, aes(x = x, y = y, z = log10(HAPTO)), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m7 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(PROCHL)), data = fd) +
    scale_fill_gradientn(name = "Prochlorococcus\nbiomass\nlog10(mg.m-3)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .2, size = .4, aes(x = x, y = y, z = log10(PROCHL)), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m8 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(PROKAR)), data = fd) +
    scale_fill_gradientn(name = "Prokaryote\nbiomass\nlog10(mg.m-3)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .2, size = .4, aes(x = x, y = y, z = log10(PROKAR)), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m9 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(NPPv2)), data = fd) +
    scale_fill_gradientn(name = "NPP\nlog10(mgC.m-2.d-1)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .25, size = .4, aes(x = x, y = y, z = log10(NPPv2)), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m10 <- ggplot() + geom_raster(aes(x = x, y = y, fill = Slope2), data = fd) +
    scale_fill_gradientn(name = "PSD slope", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .5, size = .4, aes(x = x, y = y, z = Slope2), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m11 <- ggplot() + geom_raster(aes(x = x, y = y, fill = eratio), data = fd[fd$eratio < .4,]) +
    scale_fill_gradientn(name = "E ratio", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .05, size = .4, aes(x = x, y = y, z = eratio), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m12 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(POCflux)), data = fd) +
    scale_fill_gradientn(name = "POC flux\nlog10(mgC.m-2.d-1)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .3, size = .4, aes(x = x, y = y, z = log10(POCflux)), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m13 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(FPOCex)), data = fd) +
    scale_fill_gradientn(name = "FPOCex\nlog10(mgC.m-2.d-1)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .3, size = .4, aes(x = x, y = y, z = log10(FPOCex)), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
panel1 <- ggarrange(m1,m2,m3,m4,m5,m6,m7,m8, ncol = 2, nrow = 4, align = "h", labels = letters)
panel2 <- ggarrange(m10,m9,m12,m13,m11, ncol = 2, nrow = 3, align = "v", labels = letters[9:20])

ggsave(plot = panel1, filename = "Fig.X_panel_BCP.1_27.11.23.jpg", dpi = 300, width = 12, height = 10)
ggsave(plot = panel1, filename = "Fig.X_panel_BCP.1_27.11.23.pdf", dpi = 300, width = 12, height = 10)
ggsave(plot = panel2, filename = "Fig.X_panel_BCP.2_27.11.23.jpg", dpi = 300, width = 12, height = 7)
ggsave(plot = panel2, filename = "Fig.X_panel_BCP.2_27.11.23.pdf", dpi = 300, width = 12, height = 7)


### BONUS: the fields from Clements et al., 2023
m14 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(Biovolume)), data = fd) +
    scale_fill_gradientn(name = "PSD biovolume\nlog10(ppm)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .2, size = .4, aes(x = x, y = y, z = log10(Biovolume)), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m15 <- ggplot() + geom_raster(aes(x = x, y = y, fill = Slope), data = fd) +
    scale_fill_gradientn(name = "PSD slope", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .25, size = .4, aes(x = x, y = y, z = Slope), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

m16 <- ggplot() + geom_raster(aes(x = x, y = y, fill = log10(Export)), data = fd) +
    scale_fill_gradientn(name = "Export flux\nlog10(mgC.m-2.d-1)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .2, size = .4, aes(x = x, y = y, z = log10(Export)), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

panel3 <- ggarrange(m14,m15,m16, ncol = 1, nrow = 3, align = "hv", labels = letters)
ggsave(plot = panel3, filename = "Fig.X_panel_BCP.3_Clements&al._27.11.23.jpg", dpi = 300, width = 7, height = 9)


### Compare export fluxes estimates beteen Clements et al., 2023 and DeVries & Weber 2017
exp <- ggplot(fd[which(log10(fd$POCflux) > 1),], aes(x = log10(POCflux), y = log10(Export), colour = abs(y))) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "left", size = 3) + 
  geom_abline(slope = 1, linetype = "dashed") + 
  ylab("POC export flux (Clements et al., 2023)\nlog10(mgC.m-2.d-1)") +
  xlab("POC export flux (deVries & Weber, 2017)\nlog10(mgC.m-2.d-1)") + theme_bw()

ggsave(plot = exp, filename = "Fig.X_exports_27.11.23.jpg", dpi = 300, width = 6.5, height = 5)


### ------------------------------------------------------------------

### 3) PCA plot and maps highlighting the covariance between div indices and EF proxies

names2keep <- colnames(fd)[c(4:9,11,12,14,16,21:25,29,30,32:35,39,40)] ; names2keep
data4pca <- fd[,c("x","y",names2keep)]
# Log transform
data4pca[,c(13:25)] <- log10(data4pca[,c(13:25)])
colnames(data4pca)[3:25] <- c("Richness","Faith","SES Faith","Trait dissim.","Trait turnover","Trait nestedness",
        "FEve","FDis","FDiv","Rao's Q","CHL-A","DIATO","DINO","GREEN","HAPTO","PROCHL",
        "PROKAR","NPP","FPOC","POC FLUX","E RATIO","SLOPE","MESOZOO")
data4pca <- na.omit(data4pca) 
# Perform PCA
pca <- PCA(X = data4pca[,c(3:length(data4pca))], scale.unit = T, graph = F, ncp = 5, quanti.sup = c(11:23))
# summary(pca)
#                        Dim.1   Dim.2   Dim.3   Dim.4
# Variance               6.857   1.640   0.779   0.398
# % of var.             68.572  16.398   7.790   3.985
# Cumulative % of var.  68.572  84.971  92.760  96.745
# plot.PCA(pca, axes = c(1,2), choix = "var") 
# plot.PCA(pca, axes = c(2,3), choix = "var") 
### OK. Make better plots and maps
eig <- data.frame(perc = pca$eig[,"percentage of variance"], nb = c(1:nrow(pca$eig))) # eig
pca1 <- paste0("PC1 (",floor(eig$perc[1]*100)/100,"%)")
pca2 <- paste0("PC2 (",floor(eig$perc[2]*100)/100,"%)")
pca3 <- paste0("PC3 (",floor(eig$perc[3]*100)/100,"%)")

# Function for nicer pca plots
library("ggrepel")

augment.PCA <- function(x, dims = c(1:4), which="col") {
  .get <- function(x, element, dims) {
    y <- as.data.frame(x[[element]]$coord[,dims])
    if (nrow(y) == 0) {
      y <- NULL
    } else {
      y$type <- element
    }
    return(y)
  }
  if (which == "col") {
    y <- rbind(.get(x, "var", dims), .get(x, "quanti.sup", dims))
  } else {
    y <- rbind(.get(x, "ind", dims), .get(x, "quali.sup", dims))
  }
  y$var <- row.names(y)
  row.names(y) <- NULL
  return(y)
}

pcad <- augment.PCA(pca)
# str(pcad)

p1 <- ggplot(pcad, aes(colour = factor(type))) + coord_fixed() + 
  annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
  geom_segment(aes(x=0, xend = Dim.1, y=0, yend = Dim.2), arrow=arrow(angle=20, length=unit(0.01,"npc"))) +
  scale_colour_manual(name = "", values = c("#c51b7d","#276419"), guide = "none") + 
  geom_text_repel(aes(x=Dim.1, y=Dim.2, label=var)) + xlab(pca1) + ylab(pca2) + theme_bw()

p2 <- ggplot(pcad, aes(colour = factor(type))) + coord_fixed() + 
  annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
  geom_segment(aes(x=0, xend = Dim.2, y=0, yend = Dim.3), arrow=arrow(angle=20, length=unit(0.01,"npc"))) +
  scale_colour_manual(name = "", values = c("#c51b7d","#276419"), guide = "none") + 
  geom_text_repel(aes(x=Dim.2, y=Dim.3, label=var)) + xlab(pca2) + ylab(pca3) + theme_bw()

panel.pca <- ggarrange(p1,p2, ncol = 2, nrow = 1, align = "hv", labels = letters)
ggsave(plot = panel.pca, filename = "Fig.X_panel_PCs_FD+BCP._27.11.23.jpg", dpi = 300, width = 14, height = 7)


### And for mapping PCs
data4pca[,c("PC1","PC2","PC3")] <- pca$ind$coord[,c(1:3)]

pc1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC1), data = data4pca) +
    scale_fill_gradient2(name = "PC1", low = "#3288bd", high = "#d53e4f", mid = "white") +
    geom_contour(colour = "black", binwidth = 1.5, size = .4, aes(x = x, y = y, z = PC1), data = data4pca) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

pc2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC2), data = data4pca) +
    scale_fill_gradient2(name = "PC2", low = "#3288bd", high = "#d53e4f", mid = "white") +
    geom_contour(colour = "black", binwidth = 1.5, size = .4, aes(x = x, y = y, z = PC2), data = data4pca) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

pc3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC3), data = data4pca) +
    scale_fill_gradient2(name = "PC3", low = "#3288bd", high = "#d53e4f", mid = "white") +
    geom_contour(colour = "black", binwidth = 1, size = .4, aes(x = x, y = y, z = PC3), data = data4pca) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

panel.PCs <- ggarrange(pc1,pc2,pc3, ncol = 1, nrow = 3, align = "hv", labels = letters) # panel.PCs
ggsave(plot = panel.PCs, filename = "Fig.X_panel_maps_PCs_27.11.23.jpg", dpi = 300, width = 6.5, height = 9)

### ------------------------------------------------------------------

### 4.1) Bivariate BEF plots with facet per clusters/regions
rownames(fd) <- fd$cell_id
names2keep <- colnames(fd)[c(4:8,10:12,14,16,21:25,29,30,32:35,39,40)] ; names2keep
data4pca <- fd[,c("x","y",names2keep)]
maxou <- max(data4pca$SR)
data4pca$SR <- data4pca$SR/maxou
data4pca[,c(13:22,25)] <- log10(data4pca[,c(13:22,25)])
data4pca <- na.omit(data4pca) 
# summary(data4pca)

# Center and scale values prioro to clustering
scaled.data4pca <- data4pca
colnames(scaled.data4pca)[3:25] <- c("Species Richn.","Faith","SES Faith","Trait dissim.","Trait turnover",
        "Beta ratio","FEve","FDis","FDiv","Rao's Q","CHL-A","DIATO","DINO","GREEN","HAPTO","PROCHLO",
        "PROKAR","NPP","FPOC","POC FLUX","E RATIO","PSD SLOPE","MESOZOO")
scaled.data4pca[,c(3:length(scaled.data4pca))] <- base::scale(scaled.data4pca[,c(3:length(scaled.data4pca))], center = T, scale = T)

# Clusters 
vars_clust_ward <- hclust(dist(t(scaled.data4pca[,c(3:length(scaled.data4pca))])), method = "ward.D2") 
space_clust <- hclust(dist(scaled.data4pca[,c(3:length(scaled.data4pca))]), method = "ward.D2")
klusters <- cutree(space_clust, k = 6)
scaled.data4pca$k <- as.numeric(klusters) 

# For practical reasons, make a raster out of the clusters (alos with k = 5:8)
library("raster")
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/clusters")
maps <- lapply(c(2:8), function(k) {
    
        klusters <- cutree(space_clust, k = k)
        scaled.data4pca$k <- as.numeric(klusters) 
        spg <- scaled.data4pca[,c("x","y","k")]
        colnames(spg)[1] <- "x"
        coordinates(spg) <- ~ x + y
        gridded(spg) <- TRUE
        ras <- raster(spg)
        crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    
        save(x = ras, file = paste("raster_clusters_ward_k",k,"_28.11.23.RData", sep = ""))
    
        scaled.data4pca$x2 <- scaled.data4pca$x
        scaled.data4pca[scaled.data4pca$x < 0,"x2"] <- (scaled.data4pca[scaled.data4pca$x < 0,"x"])+360

        mollweide <- ggplot() + geom_tile(aes(x = x2, y = y, fill = factor(k)), data = scaled.data4pca) +
            geom_contour(colour = "grey25", binwidth = 1, size = .4, aes(x = x2, y = y, z = k), data = scaled.data4pca) +
            scale_fill_brewer(name = "Regions", palette = "Paired") +
            geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
            theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                 panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
            theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
            scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
            scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
            coord_map("mollweide", orientation = c(90,-180,0))
    
    return(mollweide)
    
    } # eo FUN

) # eo lapply
# length(maps) # OK    
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots/")
library("ggpubr")
panel <- ggarrange(maps[[1]],maps[[2]],maps[[3]],maps[[4]],maps[[5]],maps[[6]],maps[[7]], ncol = 2, nrow = 4, align = "hv", labels = letters) # panel.PCs
ggsave(plot = panel, filename = "Fig.X_panel_maps_klusersk2-8_19.12.23.jpg", dpi = 300, width = 8, height = 15)


# Final choice is k = 6
klusters <- cutree(space_clust, k = 6)
scaled.data4pca$k <- as.numeric(klusters) 
spg <- scaled.data4pca[,c("x","y","k")]
colnames(spg)[1] <- "x"
coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
ras <- raster(spg)
crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Plot their latitudinal distribution
fd$region <- extract(ras, fd[,c('x','y')]) # summary(factor(fd$region))
ggplot(data = fd[!is.na(fd$region),], aes(x = factor(region), y = abs(y))) + 
     scale_y_continuous(breaks = seq(0,90,5)) + geom_violin(fill = "gray") +
     geom_boxplot(fill = "white", width = .2) + theme_bw()
#
# round((summary(factor(fd[!is.na(fd$region),"region"]))/nrow(fd[!is.na(fd$region),]))*100,1)

# Define the color palette for regions
pal.clusters <- rev(parula(6)) # 1st colour is too bright --> adjust
pal.clusters[1] <- "#fbe30e"

# Plot MESO x div indices, facet per 'region'
formula1 <- y~x # summary(lm)

fd2 <- fd %>% drop_na(region)

p1 <- ggplot(fd2, aes(x = SR, y = MESOZOO, colour = factor(region))) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_manual(name = "Regions", values = pal.clusters, guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "right", size = 3, colour = "black") + 
  xlab("Species richness") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw() + facet_wrap(.~factor(region), scales = "free")

p2 <- ggplot(fd2, aes(x = Faith, y = MESOZOO, colour = factor(region) )) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_manual(name = "Regions", values = pal.clusters, guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "right", size = 3, colour = "black") + 
  xlab("Faith index") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw() + facet_wrap(.~factor(region), scales = "free")

p3 <- ggplot(fd2, aes(x = FEve, y = MESOZOO, colour = factor(region) )) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_manual(name = "Regions", values = pal.clusters, guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "right", size = 3, colour = "black") + 
  xlab("Functional eveness (FEve)") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw() + facet_wrap(.~factor(region), scales = "free")

p4 <- ggplot(fd2, aes(x = FDis, y = MESOZOO, colour = factor(region) )) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_manual(name = "Regions", values = pal.clusters, guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "right", size = 3, colour = "black") + 
  xlab("Functional dispersion (FDis)") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw() + facet_wrap(.~factor(region), scales = "free")

p5 <- ggplot(fd2, aes(x = FDiv, y = MESOZOO, colour = factor(region) )) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_manual(name = "Regions", values = pal.clusters, guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "right", size = 3, colour = "black") + 
  xlab("Functional divergence (FDiv)") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw() + facet_wrap(.~factor(region), scales = "free")

p6 <- ggplot(fd2, aes(x = Jac, y = MESOZOO, colour = factor(region) )) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_manual(name = "Regions", values = pal.clusters, guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "right", size = 3, colour = "black") + 
  xlab("Trait dissimilarity (Jaccard)") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw() + facet_wrap(.~factor(region), scales = "free")

p7 <- ggplot(fd2, aes(x = Jtu, y = MESOZOO, colour = factor(region) )) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_manual(name = "Regions", values = pal.clusters, guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "right", size = 3, colour = "black") + 
  xlab("Trait turnover") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw() + facet_wrap(.~factor(region), scales = "free")

# Save
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
ggsave(plot = p1, filename = "Fig.X_facet_clusters_sr_biomass_28.11.23.jpg", dpi = 300, height = 8, width = 12)
ggsave(plot = p2, filename = "Fig.X_facet_clusters_faith_biomass_28.11.23.jpg", dpi = 300, height = 8, width = 12)
ggsave(plot = p3, filename = "Fig.X_facet_clusters_feve_biomass_28.11.23.jpg", dpi = 300, height = 8, width = 12)
ggsave(plot = p4, filename = "Fig.X_facet_clusters_fdis_biomass_28.11.23.jpg", dpi = 300, height = 8, width = 12)
ggsave(plot = p5, filename = "Fig.X_facet_clusters_fdiv_biomass_28.11.23.jpg", dpi = 300, height = 8, width = 12)
ggsave(plot = p6, filename = "Fig.X_facet_clusters_jac_biomass_28.11.23.jpg", dpi = 300, height = 8, width = 12)
ggsave(plot = p7, filename = "Fig.X_facet_clusters_jtu_biomass_28.11.23.jpg", dpi = 300, height = 8, width = 12)


### ------------------------------------------------------------------

### 4.2) Models from Fig. 3 but when adding absolute latitude as a factor --> ANCOVA analysis

### 4.2.a) FaithxSR
### We've seen how FR scales with SR linearly:
#summary(lm(Faith ~ SR, data = fd)) # R-squared:  0.536
#summary(lm(Faith ~ poly(SR,2), data = fd)) # R-squared: 0.57 

# Now, astudy the impact of latitude on this linear relationship
#summary(lm(Faith ~ SR*abs(y), data = fd)) # R2 = 0.802 # strong interactions between SR and latitude in driving Faith
#summary(lm(Faith ~ poly(SR,2)*abs(y), data = fd)) # R-squared: 0.875 instead of 0.57 
# Increase of ((0.875-0.57)/0.57)*100 = 53.51%

# Perform ANCOVA
#summary(aov(Faith ~ SR*abs(y), data = fd))
#                Df Sum Sq Mean Sq F value Pr(>F)    
# SR              1  894.7   894.7   96430 <2e-16 ***
# abs(y)          1  327.7   327.7   35316 <2e-16 ***
# SR:abs(y)       1  115.9   115.9   12488 <2e-16 ***
#mod1 <- aov(Faith ~ poly(SR,2), data = fd) 
#mod2 <- aov(Faith ~ poly(SR,2)*abs(y), data = fd) 
#anova(mod1,mod2)
#    Res.Df  RSS Df Sum of Sq     F    Pr(>F)    
#1   35677 719.25                                 
#2   35674 207.36  3    511.89 29354 < 2.2e-16 ***


### 4.2.b) SESxSR
#summary(lm(SES.Faith ~ SR, data = fd)) # R2 = 0.59 
#summary(lm(SES.Faith ~ SR*abs(y), data = fd)) # R2 = 0.86 (instead of 0.59)
# Increase of ((0.86-0.59)/0.59)*100 = 45.76%
#summary(aov(SES.Faith ~ SR*abs(y), data = fd))
#                Df Sum Sq Mean Sq F value Pr(>F)    
# SR              1  58552   58552  150798 <2e-16 ***
# abs(y)          1  22872   22872   58906 <2e-16 ***
# SR:abs(y)       1   3268    3268    8416 <2e-16 ***
# Residuals   35676  13852       0 
# mod1 <- aov(SES.Faith ~ SR, data = fd) 
# mod2 <- aov(SES.Faith ~ SR*abs(y), data = fd) 
# anova(mod1,mod2)


### 4.2.c) FEve x SR
#summary(lm(FEve ~ SR, data = fd)) # 0.5547 
#summary(lm(FEve ~ SR*abs(y), data = fd)) # R2 = 0.76 (instead of 0.55)
# Increase of ((0.76-0.5547)/0.5547)*100 = 37.011%
#summary(aov(FEve ~ SR*abs(y), data = fd))
#                 Df Sum Sq Mean Sq F value Pr(>F)    
# SR              1  96.03   96.03   81395 <2e-16 ***
# abs(y)          1   4.41    4.41    3737 <2e-16 ***
# SR:abs(y)       1  30.58   30.58   25919 <2e-16 ***
# Residuals   35676  42.09    0.00          
# mod1 <- aov(FEve ~ SR, data = fd) 
# mod2 <- aov(FEve ~ SR*abs(y), data = fd) 
# anova(mod1,mod2)


### 4.2.d) FDis x SR
#summary(lm(FDis ~ poly(SR,2), data = fd)) # R2 = 0.7438
#summary(lm(FDis ~ poly(SR,2)*abs(y), data = fd)) # R2 = 0.7959 
# Increase of ((0.7959-0.7438)/0.7438)*100 = 7.00%
#summary(aov(FDis ~ poly(SR,2)*abs(y), data = fd))
#                 Df Sum Sq Mean Sq F value Pr(>F)    
# SR              1 1.2903  1.2903  104154 <2e-16 ***
# abs(y)          1 0.0404  0.0404    3264 <2e-16 ***
# SR:abs(y)       1 0.0828  0.0828    6681 <2e-16 ***
# Residuals   35676 0.4420  0.0000 
# mod1 <- aov(FDis ~ poly(SR,2), data = fd) 
# mod2 <- aov(FDis ~ poly(SR,2)*abs(y), data = fd) 
# anova(mod1,mod2)


### 4.2.e) FDiv x SR
#summary(lm(FDiv ~ poly(SR,2), data = fd)) # R2 = 0.67
#summary(lm(FDiv ~ poly(SR,2)*abs(y), data = fd)) # R2 = 0.85
# Increase of ((0.85-0.67)/0.67)*100 = 26.86%
#summary(aov(FDiv ~ SR*abs(y), data = fd))
#                 Df Sum Sq Mean Sq F value Pr(>F)    
# SR              1  3.283   3.283  129856 <2e-16 ***
# abs(y)          1  1.051   1.051   41587 <2e-16 ***
# SR:abs(y)       1  0.057   0.057    2267 <2e-16 ***
# Residuals   35676  0.902   0.000  
# mod1 <- aov(FDiv ~ poly(SR,2), data = fd) 
# mod2 <- aov(FDiv ~ poly(SR,2)*abs(y), data = fd) 
# anova(mod1,mod2)


### 4.2.f) Jac x SR
#summary(lm(Jac ~ poly(SR,2), data = fd)) # R2 = 0.94
#summary(lm(Jac ~ poly(SR,2)*abs(y), data = fd)) # R2 = 0.95
# Increase of ((0.95-0.94)/0.94)*100 = 1.06%
#summary(aov(Jac ~ SR*abs(y), data = fd))
#                 Df Sum Sq Mean Sq F value Pr(>F)    
# SR              1  98.13   98.13  485910 <2e-16 ***
# abs(y)          1   2.61    2.61   12938 <2e-16 ***
# SR:abs(y)       1   3.72    3.72   18405 <2e-16 ***
# Residuals   35676   7.20    0.00   
#mod1 <- aov(Jac ~ poly(SR,2), data = fd) 
#mod2 <- aov(Jac ~ poly(SR,2)*abs(y), data = fd) 
#anova(mod1,mod2)


### 4.2.g) Jtu x SR
#summary(lm(Jtu ~ poly(SR,2), data = fd)) # R2 = 0.85
#summary(lm(Jtu ~ poly(SR,2)*abs(y), data = fd)) # R2 = 0.93
# Increase of ((0.93-0.85)/0.85)*100 = 9.41%
#summary(aov(Jtu ~ SR*abs(y), data = fd))
#                Df Sum Sq Mean Sq F value Pr(>F)    
# SR              1  70.12   70.12  435596 <2e-16 ***
# abs(y)          1   6.97    6.97   43294 <2e-16 ***
# SR:abs(y)       1   4.62    4.62   28695 <2e-16 ***
# Residuals   35676   5.74    0.00 
#mod1 <- aov(Jtu ~ poly(SR,2), data = fd) 
#mod2 <- aov(Jtu ~ poly(SR,2)*abs(y), data = fd) 
#anova(mod1,mod2)


### 4.2.h) Jne x SR
#summary(lm(Jne ~ poly(SR,2), data = fd)) # R2 = 0.2632
#summary(lm(Jne ~ poly(SR,2)*abs(y), data = fd)) # R2 = 0.5101
# Increase of ((0.5101-0.2632)/0.2632)*100 = 93.81%
#summary(aov(Jne ~ SR*abs(y), data = fd))
#                 Df Sum Sq Mean Sq F value Pr(>F)    
# SR              1  2.446  2.4461 11808.1 <2e-16 ***
# abs(y)          1  0.986  0.9856  4757.7 <2e-16 ***
# SR:abs(y)       1  0.027  0.0271   130.9 <2e-16 ***
# Residuals   35676  7.390  0.0002  
#mod1 <- aov(Jne ~ poly(SR,2), data = fd) 
#mod2 <- aov(Jne ~ poly(SR,2)*abs(y), data = fd) 
#anova(mod1,mod2)


### With mesozooplankton biomass
# SR
#summary(lm(MESOZOO ~ SR, data = fd)) # R2 = 0.46 
#summary(lm(MESOZOO ~ SR*abs(y), data = fd)) # R2 = 0.50 
# ((0.50-0.46)/0.46)*100 = 8.7% 

# Faith
#summary(lm(MESOZOO ~ Faith, data = fd)) # R2 = 0.43
#summary(lm(MESOZOO ~ Faith*abs(y), data = fd)) # R2 = 0.58
# ((0.58-0.43)/0.43)*100 = 34.9%
#summary(aov(MESOZOO ~ Faith*abs(y), data = fd))
#                 Df Sum Sq Mean Sq F value Pr(>F)    
# Faith            1  21725   21725   34902 <2e-16 ***
# abs(y)           1   2121    2121    3408 <2e-16 ***
# Faith:abs(y)     1   5802    5802    9321 <2e-16 ***
# Residuals    34544  21502       1   

# FEve
#summary(lm(MESOZOO ~ FEve, data = fd)) # R2 = 0.1839 
#summary(lm(MESOZOO ~ FEve*abs(y), data = fd)) # R2 = 0.1866 
# ((0.1866-0.1839)/0.1839)*100 = 1.5%
# mod1 <- aov(MESOZOO ~ FEve, data = fd) 
# mod2 <- aov(MESOZOO ~ FEve*abs(y), data = fd) 
# anova(mod1,mod2)

# FDis
#summary(lm(MESOZOO ~ FDis, data = fd)) # R2 = 0.62
#summary(lm(MESOZOO ~ FDis*abs(y), data = fd)) # R2 = 0.66
# ((0.66-0.62)/0.62)*100 = 6.45%
#mod1 <- aov(MESOZOO ~ FDis, data = fd) 
#mod2 <- aov(MESOZOO ~ FDis*abs(y), data = fd) 
#anova(mod1,mod2)

# FDiv
#summary(lm(MESOZOO ~ FDiv, data = fd)) # R2 = 0.1634 
#summary(lm(MESOZOO ~ FDiv*abs(y), data = fd)) # R2 = 0.2737
# ((0.2737-0.1634)/0.1634)*100 = 67.5%


# Jac
#summary(lm(MESOZOO ~ Jac, data = fd)) # R2 = 0.4792
#summary(lm(MESOZOO ~ Jac*abs(y), data = fd)) # R2 = 0.5939 
# ((0.5939-0.4792)/0.4792)*100 = 23.9%

#summary(aov(MESOZOO ~ Jac*abs(y), data = fd))
#                Df Sum Sq Mean Sq F value Pr(>F)    
# Jac             1  24511   24511   40770 <2e-16 ***
# abs(y)          1   4791    4791    7968 <2e-16 ***
# Jac:abs(y)      1   1081    1081    1798 <2e-16 ***
# Residuals   34544  20768       1  

# Jtu
#summary(lm(MESOZOO ~ Jtu, data = fd)) # R2 = 0.41
#summary(lm(MESOZOO ~ Jtu*abs(y), data = fd)) # R2 = 0.56
#summary(aov(MESOZOO ~ Jtu*abs(y), data = fd))
#                Df Sum Sq Mean Sq F value Pr(>F)    
# Jtu             1  20884   20884   31963 <2e-16 ***
# abs(y)          1   6449    6449    9870 <2e-16 ***
# Jtu:abs(y)      1   1248    1248    1910 <2e-16 ***
# Residuals   34544  22570       1    

# Jne
#summary(lm(MESOZOO ~ Jne, data = fd)) # R2 = 0.16
#summary(lm(MESOZOO ~ Jne*abs(y), data = fd)) # R2 = 0.31
#summary(aov(MESOZOO ~ Jne*abs(y), data = fd))
#                Df Sum Sq Mean Sq F value Pr(>F)    
# Jne             1   8208    8208    8082 <2e-16 ***
# abs(y)          1   3530    3530    3476 <2e-16 ***
# Jne:abs(y)      1   4329    4329    4262 <2e-16 ***
# Residuals   34544  35083       1 


### ------------------------------------------------------------------

### 5) Maps of uncertainty/sd() of contemporary div indices: intra-annual variability (month to month)
### 18/04/24: Same but adding relative standard deviation in the summarize()

# 5.1) Faith's and SR
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
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) # length(unique(tab$cell_id))
# Compute monthly mean indices first
mon.sr.faith <- data.frame(tab %>% group_by(cell_id,month) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.faith = mean(Faith, na.rm = T),
            mean.rich = mean(SR, na.rm = T))
)
# Then annual mean and seasonality (std across month)
var.sr.faith <- data.frame(mon.sr.faith %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            ann.faith = mean(mean.faith, na.rm = T), sd.faith = sd(mean.faith, na.rm = T), rel.sd.faith = (sd.faith*100)/ann.faith,
            ann.rich = mean(mean.rich, na.rm = T), sd.rich = sd(mean.rich, na.rm = T), rel.sd.rich = (sd.rich*100)/ann.rich)
)

# 5.2) SES Faith
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
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) # length(unique(tab$cell_id))
colnames(tab)[1] <- "SR"
colnames(tab)[2] <- "FR"
mon.ses <- data.frame(tab %>% group_by(cell_id,month) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.ses = mean(pd.obs.z, na.rm = T)
    )
) # eo ddf
var.ses <- data.frame(mon.ses %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            ann.ses = mean(mean.ses, na.rm = T), sd.ses = sd(mean.ses, na.rm = T), rel.sd.ses = (sd.ses*100)/ann.ses)
)


# 5.3) Beta div indices (Jac/Jtu/Jne) 
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
# Compute monthly means
mon.beta <- data.frame(tab %>% group_by(cell_id,month) %>% 
        summarize(x = unique(x), y = unique(y),
        jac = mean(beta.jac, na.rm = T),
        jtu = mean(beta.jtu, na.rm = T),
        jne = mean(beta.jne, na.rm = T)
    ) 
) # eo ddf
var.beta <- data.frame(mon.beta %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            ann.jac = mean(jac, na.rm = T), sd.jac = sd(jac, na.rm = T), rel.sd.jac = (sd.jac*100)/ann.jac,
            ann.jtu = mean(jtu, na.rm = T), sd.jtu = sd(jtu, na.rm = T), rel.sd.jtu = (sd.jtu*100)/ann.jtu,
            ann.jne = mean(jne, na.rm = T), sd.jne = sd(jne, na.rm = T), rel.sd.jne = (sd.jne*100)/ann.jne
        )
)
# summary(var.beta)


# 5.4) dbFD indices
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid/")
files <- dir()[grep("_baseline",dir())]
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
# Compute monthly means
mon.dbFD <- data.frame(tab %>% group_by(cell_id,month) %>% 
        summarize(x = unique(x), y = unique(y),
            FEve.avg = mean(FEve, na.rm = T),
            FDis.avg = mean(FDis, na.rm = T),
            FDiv.avg = mean(FDiv, na.rm = T)
        )
)
var.dbFD <- data.frame(mon.dbFD %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            ann.feve = mean(FEve.avg, na.rm = T), sd.feve = sd(FEve.avg, na.rm = T), rel.sd.feve = (sd.feve*100)/ann.feve,
            ann.fdis = mean(FDis.avg, na.rm = T), sd.fdis = sd(FDis.avg, na.rm = T), rel.sd.fdis = (sd.fdis*100)/ann.fdis,
            ann.fdiv = mean(FDiv.avg, na.rm = T), sd.fdiv = sd(FDiv.avg, na.rm = T), rel.sd.fdiv = (sd.fdiv*100)/ann.fdiv
    )
)
# summary(var.dbFD)

### Maps to make: SR, Faith, SES Faith, FEve, FDis, FDiv, Jac, Jtu, Jne: 9x9 maps panel
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rel.sd.rich), data = var.sr.faith) +
    scale_fill_gradientn(name = "Intra-annual\nvariability\n(%Richness)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = 10, size = .4, aes(x = x, y = y, z = rel.sd.rich), data = var.sr.faith) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rel.sd.faith), data = var.sr.faith) +
    scale_fill_gradientn(name = "Intra-annual\nvariability\n(%Faith)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = 5, size = .4, aes(x = x, y = y, z = rel.sd.faith), data = var.sr.faith) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

# map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.ses), data = var.ses) +
#     scale_fill_gradientn(name = "Intra-annual\nvariability\n(SES Faith)", colours = parula(100), guide = "colourbar") +
#     geom_contour(colour = "black", binwidth = .5, size = .4, aes(x = x, y = y, z = sd.ses), data = var.ses) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rel.sd.feve), data = var.dbFD) +
    scale_fill_gradientn(name = "Intra-annual\nvariability\n(%FEve)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = 10, size = .4, aes(x = x, y = y, z = rel.sd.feve), data = var.dbFD) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rel.sd.fdis), data = var.dbFD) +
    scale_fill_gradientn(name = "Intra-annual\nvariability\n(%FDis)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = 1, size = .4, aes(x = x, y = y, z = rel.sd.fdis), data = var.dbFD) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rel.sd.fdiv), data = var.dbFD) +
    scale_fill_gradientn(name = "Intra-annual\nvariability\n(%FDiv)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = 1, size = .4, aes(x = x, y = y, z = rel.sd.fdiv), data = var.dbFD) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map7 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rel.sd.jac), data = var.beta) +
    scale_fill_gradientn(name = "Intra-annual\nvariability\n(%Jaccard)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = 5, size = .4, aes(x = x, y = y, z = rel.sd.jac), data = var.beta) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map8 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rel.sd.jtu), data = var.beta) +
    scale_fill_gradientn(name = "Intra-annual\nvariability\n(%Turnover)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = 5, size = .4, aes(x = x, y = y, z = rel.sd.jtu), data = var.beta) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map9 <- ggplot() + geom_raster(aes(x = x, y = y, fill = rel.sd.jne), data = var.beta) +
    scale_fill_gradientn(name = "Intra-annual\nvariability\n(%Nestedness)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = 20, size = .4, aes(x = x, y = y, z = rel.sd.jne), data = var.beta) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
figX <- ggarrange(map1,map2,map4,map5,map6,map7,map8,map9, align = 'hv', ncol = 2, nrow = 4, labels = letters)
#ggsave(plot = figX, filename = "Fig.X_map_sd_monthly_fd_indices_19.12.23.jpg", dpi = 300, width = 12, height = 12)
#ggsave(plot = figX, filename = "Fig.X_map_sd_monthly_fd_indices_27.11.23.pdf", dpi = 300, width = 16.5, height = 7)
ggsave(plot = figX, filename = "Fig.X_map_relative_sd_monthly_fd_indices_18.04.24.jpg", dpi = 300, width = 12, height = 12)


### ------------------------------------------------------------------

### 6) Maps of uncertainty/sd() of contemporary div indices: SDM-driven variability (month to month)
### Re-use the same code generated above

# 6.1) Faith's and SR
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/Faith")
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
# Add cell id
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) # length(unique(tab$cell_id))
### Compute mean annual indices
mon.sr.faith <- data.frame(tab %>% group_by(cell_id,SDM) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.faith = mean(Faith, na.rm = T),
            mean.rich = mean(SR, na.rm = T))
)
var.sr.faith <- data.frame(mon.sr.faith %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            ann.faith = mean(mean.faith, na.rm = T), sd.faith = sd(mean.faith, na.rm = T),
            ann.rich = mean(mean.rich, na.rm = T), sd.rich = sd(mean.rich, na.rm = T))
)

# 6.2) SES Faith
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/SES")
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
# Add cell id
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) # length(unique(tab$cell_id))
colnames(tab)[1] <- "SR"
colnames(tab)[2] <- "FR"
mon.ses <- data.frame(tab %>% group_by(cell_id,SDM) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.ses = mean(pd.obs.z, na.rm = T)
    )
) # eo ddf
var.ses <- data.frame(mon.ses %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            ann.ses = mean(mean.ses, na.rm = T), sd.ses = sd(mean.ses, na.rm = T)
        )
)

# 6.3) Beta div indices (Jac/Jtu/Jne) 
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
# Compute mean annual indices
mon.beta <- data.frame(tab %>% group_by(cell_id,SDM) %>% 
        summarize(x = unique(x), y = unique(y),
        jac = mean(beta.jac, na.rm = T),
        jtu = mean(beta.jtu, na.rm = T),
        jne = mean(beta.jne, na.rm = T)
    ) 
) # eo ddf
var.beta <- data.frame(mon.beta %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            ann.jac = mean(jac, na.rm = T), sd.jac = sd(jac, na.rm = T),
            ann.jtu = mean(jtu, na.rm = T), sd.jtu = sd(jtu, na.rm = T),
            ann.jne = mean(jne, na.rm = T), sd.jne = sd(jac, na.rm = T)
    )
)

# 6.4) dbFD indices
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid/")
files <- dir()[grep("_baseline",dir())]
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
mon.dbFD <- data.frame(tab %>% group_by(cell_id,SDM) %>% 
        summarize(x = unique(x), y = unique(y),
            FEve.avg = mean(FEve, na.rm = T),
            FDis.avg = mean(FDis, na.rm = T),
            FDiv.avg = mean(FDiv, na.rm = T)
        )
)
var.dbFD <- data.frame(mon.dbFD %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            ann.feve = mean(FEve.avg, na.rm = T), sd.feve = sd(FEve.avg, na.rm = T),
            ann.fdis = mean(FDis.avg, na.rm = T), sd.fdis = sd(FDis.avg, na.rm = T),
            ann.fdiv = mean(FDiv.avg, na.rm = T), sd.fdiv = sd(FDiv.avg, na.rm = T)
    )
)


### Maps to make: SR, Faith, SES Faith, FEve, FDis, FDiv, Jac, Jtu, Jne: 9x9 maps panel
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.rich), data = var.sr.faith) +
    scale_fill_gradientn(name = "SDM-driven\nvariability\n(Richness)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = 10, size = .4, aes(x = x, y = y, z = sd.rich), data = var.sr.faith) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.faith), data = var.sr.faith) +
    scale_fill_gradientn(name = "SDM-driven\nvariability\n(Faith)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .1, size = .4, aes(x = x, y = y, z = sd.faith), data = var.sr.faith) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.ses), data = var.ses) +
    scale_fill_gradientn(name = "SDM-driven\nvariability\n(SES Faith)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .5, size = .4, aes(x = x, y = y, z = sd.ses), data = var.ses) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.feve), data = var.dbFD) +
    scale_fill_gradientn(name = "SDM-driven\nvariability\n(FEve)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .05, size = .4, aes(x = x, y = y, z = sd.feve), data = var.dbFD) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map5 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.fdis), data = var.dbFD) +
    scale_fill_gradientn(name = "SDM-driven\nvariability\n(FDis)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .003, size = .4, aes(x = x, y = y, z = sd.fdis), data = var.dbFD) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map6 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.fdiv), data = var.dbFD) +
    scale_fill_gradientn(name = "SDM-driven\nvariability\n(FDiv)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .005, size = .4, aes(x = x, y = y, z = sd.fdiv), data = var.dbFD) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map7 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.jac), data = var.beta) +
    scale_fill_gradientn(name = "SDM-driven\nvariability\n(Jaccard)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .025, size = .4, aes(x = x, y = y, z = sd.jac), data = var.beta) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map8 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.jtu), data = var.beta) +
    scale_fill_gradientn(name = "SDM-driven\nvariability\n(Turnover)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .02, size = .4, aes(x = x, y = y, z = sd.jtu), data = var.beta) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map9 <- ggplot() + geom_raster(aes(x = x, y = y, fill = sd.jne), data = var.beta) +
    scale_fill_gradientn(name = "SDM-driven\nvariability\n(Nestedness)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .025, size = .4, aes(x = x, y = y, z = sd.jne), data = var.beta) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
figX <- ggarrange(map1,map2,map3,map4,map5,map6,map7,map8,map9, align = 'hv', ncol = 2, nrow = 5, labels = letters)
ggsave(plot = figX, filename = "Fig.X_map_sd_SDM_fd_indices_19.12.23.jpg", dpi = 300, width = 12, height = 12)
#ggsave(plot = figX, filename = "Fig.X_map_sd_SDM_fd_indices_27.11.23.pdf", dpi = 300, width = 16.5, height = 7)


### ------------------------------------------------------------------

### 7.1) Maps of uncertainty/sd() for changes in future div indices  
### Need to adjust the code here because of the way the predictions were stored 
indices <- c("SR","FR","FEve","FDis","FDiv","Jac","Jne","Jtu")
# i <- "SR" # for testing

res <- lapply(indices, function(i) {
    
    message(paste("\n","Calculating ∆ values for ",i,"\n", sep = ""))
    
    if(i %in% c("SR","FR")) {
        setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/Faith")
    } else if (i == "SES.FR") {
        setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/SES")
    } else if (i %in% c("Jac","Jne","Jtu")) {
        setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/beta.div")
    } else {
        setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid")    
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
    ens <- data.frame(tab %>% group_by(cell_id,ESM,SDM) %>% 
            summarize(x = unique(x), y = unique(y),
                mean.perc = mean(perc, na.rm = T)
        )
    ) # eo ddf
    # Derive variability out of SDM+ESM
    var.sdm <- data.frame(ens %>% group_by(cell_id) %>% 
            summarize(x = unique(x), y = unique(y),
                mean = mean(mean.perc, na.rm = T),
                std = sd(mean.perc, na.rm = T)
        )
    ) # eo ddf
    # dim(var.sdm); summary(var.sdm)
    
    # Map for stdev
    if(i == "SR") {
        title <- paste("Uncertainty\n(∆","Richness",")", sep = "")
    } else if(i == "FR") {
        title <- paste("Uncertainty\n(∆","Faith",")", sep = "")
    } else if(i == "FEve") {
        title <- paste("Uncertainty\n(∆","FEve",")", sep = "")
    } else if(i == "FDis") {
        title <- paste("Uncertainty\n(∆","FDis",")", sep = "")
    } else if(i == "FDiv") {
        title <- paste("Uncertainty\n(∆","FDiv",")", sep = "")
    } else if(i == "Jac") {
        title <- paste("Uncertainty\n(∆","Jaccard",")", sep = "")
    } else if(i == "Jtu") {
        title <- paste("Uncertainty\n(∆","Turnover",")", sep = "")
    } else if(i == "Jne") {
        title <- paste("Uncertainty\n(∆","Nestedness",")", sep = "")
    } 
    
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = std), data = var.sdm) +
         scale_fill_gradientn(name = title, guide = "colourbar", colours = parula(50)) +
         geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
         coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
         scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
         scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)    
    
    return(map)
    
    } # eo FUN
    
) #  eo lapply

# Arrange on grid
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
figX <- ggarrange(res[[1]],res[[2]],res[[3]],res[[4]],res[[5]],res[[6]],res[[7]],res[[8]], align = "hv", ncol = 2, nrow = 4, labels = letters)
ggsave(plot = figX, filename = "Fig.X_map_sd_future_indices_19.12.23.jpg", dpi = 300, width = 12, height = 10)
#ggsave(plot = figX, filename = "Fig.X_map_sd_future_indices_27.11.23.pdf", dpi = 300, width = 12, height = 10) # does not return the '∆'


### ------------------------------------------------------------------

### 12/01/24: 7.2) Maps of agreement of changes in future div indices  
indices <- c("SR","FR","FEve","FDis","FDiv","Jac","Jne","Jtu")
i <- "FDiv" # for testing
 
res <- lapply(indices, function(i) {
    
    message(paste("\n","Calculating ∆ values for ",i,"\n", sep = ""))
    
    if(i %in% c("SR","FR")) {
        setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/Faith")
    } else if (i == "SES.FR") {
        setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/SES")
    } else if (i %in% c("Jac","Jne","Jtu")) {
        setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/beta.div")
    } else {
        setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid")    
    } # eo 1st if else loop

    files <- dir()[grep(i,dir())]
    # f <- files[1]
    
    res <- mclapply(files, function(f) {
                d <- get(load(f))
                filename <- str_replace_all(f,".Rdata","")
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
        
    # Compute ensemble mean for 'perc'
    ens <- data.frame(tab %>% group_by(cell_id,ESM,SDM) %>% 
            summarize(x = unique(x), y = unique(y),
                mean.perc = mean(perc, na.rm = T)
        )
    ) # eo ddf
    
    ens <- ens %>% drop_na(mean.perc)
    
    # Count positve and negative values and compute agreement (% of )
    agree <- data.frame(ens %>% group_by(cell_id) %>% 
            summarize(x = unique(x), y = unique(y),
                n = n(),
                n.posi = sum(mean.perc > 0),
                n.nega = sum(mean.perc < 0),
                p.posi = (n.posi/n)*100,
                p.nega = (n.nega/n)*100
        )
    ) # eo ddf
    # dim(agree); summary(agree)
    
    agree$agree <- 0
    agree[agree$p.posi > 75 | agree$p.nega > 75,"agree"] <- 1
    # summary(agree$agree)    
        
    # Map for stdev
    if(i == "SR") {
        title <- paste("75% Agreement in\n∆","Richness", sep = "")
    } else if(i == "FR") {
        title <- paste("75% Agreement in\n∆","Faith", sep = "")
    } else if(i == "FEve") {
        title <- paste("75% Agreement in\n∆","FEve", sep = "")
    } else if(i == "FDis") {
        title <- paste("75% Agreement in\n∆","FDis", sep = "")
    } else if(i == "FDiv") {
        title <- paste("75% Agreement in\n∆","FDiv", sep = "")
    } else if(i == "Jac") {
        title <- paste("75% Agreement in\n∆","Jaccard", sep = "")
    } else if(i == "Jtu") {
        title <- paste("75% Agreement in\n∆","Turnover", sep = "")
    } else if(i == "Jne") {
        title <- paste("75% Agreement in\n∆","Nestedness", sep = "")
    } 
    
    map <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(agree)), data = agree[agree$n == 15,]) +
         scale_fill_manual(name = title, values = c('#d53e4f','#3288bd')) +
         geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
         coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
         scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
         scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)    
    
    return(map)
    
    } # eo FUN
    
) #  eo lapply

# Arrange on grid
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
figX <- ggarrange(res[[1]],res[[2]],res[[3]],res[[4]],res[[5]],res[[6]],res[[7]],res[[8]], align = "hv", ncol = 2, nrow = 4, labels = letters)
ggsave(plot = figX, filename = "Fig.X_map_∆_agreement_12.01.24.jpg", dpi = 300, width = 12, height = 10)


### ------------------------------------------------------------------

### 8) Maps of p-value distributions for SES Faith; map of future SES Faith pattern to show consistency
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/SES")
files <- dir()[grep("baseline",dir())] ; files
# f <- files[2]
res <- mclapply(files, function(f) {
            d <- get(load(f))
            filename <- str_replace_all(f,".Rdata","")
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
            return(d)
        }, mc.cores = 5
) # eo mclapply
tab <- bind_rows(res)
rm(res); gc()
tab$cell_id <- factor(paste(tab$x, tab$y, sep = "_")) # length(unique(tab$cell_id))
colnames(tab)[1] <- "SR"
colnames(tab)[2] <- "FR"

# Report basic statistics of 'pd.obs.p'
#median(tab$pd.obs.p); IQR(tab$pd.obs.p)

distrib1 <- ggplot(tab, aes(pd.obs.p)) + geom_density(position = "stack", fill = 'grey') + 
     geom_vline(xintercept = .05, colour = "red", linetype = "dashed") +
     xlab("P value") + ylab("Density") + theme_bw()

# Separate tropics from extra tropics and add on density plot with fill as a new aes()
tab$region <- NA
tab[which(abs(tab$y) <= 30),"region"] <- "Tropics (0-30°)"
tab[which(abs(tab$y) > 30),"region"] <- "Extra-tropics (>30°)"

distrib2 <- ggplot(tab, aes(x = pd.obs.p, fill = factor(region))) + geom_density(position = "fill", bw = .001) + 
     scale_fill_manual(name = "", values = c("#3288bd","#d53e4f")) + 
     geom_vline(xintercept = .05, colour = "white", linetype = "dashed") +
     xlab("P value") + ylab("Density") + theme_bw() + theme(legend.position = "none")

### Per SDM?
distrib3 <- ggplot(tab, aes(x = pd.obs.p, fill = factor(region))) + geom_density(position = "fill", bw = .001) + 
     scale_fill_manual(name = "", values = c("#3288bd","#d53e4f")) + 
     geom_vline(xintercept = .05, colour = "white", linetype = "dashed") +
     xlab("P value") + ylab("Density") + theme_bw() + theme(legend.position = "bottom") +
     facet_wrap(.~factor(SDM))

panel.distrib <- ggarrange(distrib1,distrib2,distrib3, align = "hv", ncol = 1, nrow = 3, labels = letters)

### Compute annual means
ann.ses.faith <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
        mean.ses = mean(pd.obs.z, na.rm = T),
        sd.ses = sd(pd.obs.z, na.rm = T),
        N = n(), N.pval = length(pd.obs.p[pd.obs.p < 0.05]),
        freq.pval = (N.pval/N)*100
    )
)
# summary(ann.ses.faith)

map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = freq.pval), data = ann.ses.faith) +
     scale_fill_gradientn(name = "Frequency of SES Faith\nwith P < 0.05", guide = "colourbar", colours = parula(100)) +
     geom_contour(colour = "black", binwidth = 20, size = .4, aes(x = x, y = y, z = freq.pval), data = ann.ses.faith) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)    


### Check how plots change for future conditions?
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/SES")
files <- dir()[grep("2100-2000",dir())]
# f <- files[1]
res <- mclapply(files, function(f) {
            d <- get(load(f))
            filename <- str_replace_all(f,".Rdata","")
            d$ESM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
            d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
            d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[7]
            return(d)
        }, mc.cores = 5
) # eo mclapply
fut.tab <- bind_rows(res)
rm(res); gc()
# Re-rotate longitudes back to -179.5/+179.5
fut.tab$x2 <- fut.tab$x
fut.tab[fut.tab$x > 179.5,"x2"] <- (fut.tab[fut.tab$x > 179.5,"x"]) - 360
fut.tab$cell_id <- factor(paste(fut.tab$x2, fut.tab$y, sep = "_"))
fut.tab <- fut.tab[order(fut.tab$cell_id),]

fut.distrib1 <- ggplot(fut.tab, aes(pd.obs.p)) + geom_density(position = "stack", fill = 'grey') + 
     geom_vline(xintercept = .05, colour = "red", linetype = "dashed") +
     xlab("P value") + ylab("Density") + theme_bw()

# Separate tropics from extra tropics and add on density plot with fill as a new aes()
fut.tab$region <- NA
fut.tab[which(abs(fut.tab$y) <= 30),"region"] <- "Tropics (0-30°)"
fut.tab[which(abs(fut.tab$y) > 30),"region"] <- "Extra-tropics (>30°)"

fut.distrib2 <- ggplot(fut.tab, aes(x = pd.obs.p, fill = factor(region))) + geom_density(position = "fill", bw = .001) + 
     scale_fill_manual(name = "", values = c("#3288bd","#d53e4f")) + 
     geom_vline(xintercept = .05, colour = "white", linetype = "dashed") +
     xlab("P value") + ylab("Density") + theme_bw() + theme(legend.position = "none")

### Per SDM?
fut.distrib3 <- ggplot(fut.tab, aes(x = pd.obs.p, fill = factor(region))) + geom_density(position = "fill", bw = .001) + 
     scale_fill_manual(name = "", values = c("#3288bd","#d53e4f")) + 
     geom_vline(xintercept = .05, colour = "white", linetype = "dashed") +
     xlab("P value") + ylab("Density") + theme_bw() + theme(legend.position = "bottom") +
     facet_wrap(.~factor(SDM))

panel.fut.distrib <- ggarrange(fut.distrib1,fut.distrib2,fut.distrib3, align = "hv", ncol = 1, nrow = 3, labels = letters)

fut.ses.faith <- data.frame(fut.tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x2), y = unique(y),
            mean.ses = mean(pd.obs.z, na.rm = T),
            sd.ses = sd(pd.obs.z, na.rm = T),
            N = n(),
            N.pval = length(pd.obs.p[pd.obs.p < 0.05]),
            freq.pval = (N.pval/N)*100
    )
)
# summary(fut.ses.faith)

map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = freq.pval), data = fut.ses.faith) +
     scale_fill_gradientn(name = "Frequency of SES Faith\nwith P < 0.05", guide = "colourbar", colours = parula(100)) +
     geom_contour(colour = "black", binwidth = 20, size = .4, aes(x = x, y = y, z = freq.pval), data = fut.ses.faith) +
     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)    

panel.maps <- ggarrange(map1,map2, align = "hv", ncol = 1, nrow = 2, labels = letters)

### Save plots 
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots/")
ggsave(plot = panel.maps, filename = "Fig.X_maps_freq.pval_ses.faith_28_11_23.jpg", dpi = 300, width = 8, height = 6)
ggsave(plot = panel.distrib, filename = "Fig.X_panel_freq.pval_ses.faith_28_11_23.jpg", dpi = 300, width = 7, height = 10)
ggsave(plot = panel.fut.distrib, filename = "Fig.X_panel_fut_freq.pval_ses.faith_28_11_23.jpg", dpi = 300, width = 7, height = 10)

  
### ------------------------------------------------------------------

### 9) Boxplots illustrating future changes in diversity metrics per regions
indices <- c("SR","FR","FEve","FDis","FDiv","Jac","Jtu","Jne")

res <- lapply(indices, function(i) {
    
        message(paste("\n","Mapping ∆",i,"\n", sep = ""))
    
        if(i %in% c("SR","FR")) {
            setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/Faith")
        } else if (i == "SES.FR") {
            setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/SES")
        } else if (i %in% c("Jac","Jne","Jtu")) {
            setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/beta.div")
        } else {
            setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid")    
        } # eo 1st if else loop
    
        # Vector files of interest
        files <- dir()[grep(i,dir())] 
    
        res <- mclapply(files, function(f) {
                d <- get(load(f))
                filename <- str_replace_all(f,".Rdata","")
                d$ESM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
                d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
                d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[7]
                return(d)
            }, mc.cores = 10
        ) # eo mclapply
        tab <- bind_rows(res)
        rownames(tab) <- NULL
        rm(res); gc()
        
        # Compute ensemble mean for 'perc'
        ens <- data.frame(tab %>% group_by(cell_id) %>% 
                summarize(x = unique(x), y = unique(y),
                    mean = mean(perc, na.rm = T),
                    std = sd(perc, na.rm = T)
                )
        ) # eo ddf
        # summary(ens)
        ens$index <- i
        
        return(ens)
    
    } # eo FUN
    
) # eo lapply 
ddf <- bind_rows(res)
rm(res); gc()

### Load the raster containing the 6 regions and extract at ddf position
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/clusters")
ras <- get(load("raster_clusters_ward_k6_28.11.23.RData"))
ddf$region <- extract(ras, ddf[,c("x","y")])
ddf2 <- ddf %>% drop_na(region)

# Defien the colour palette for the regions
pal.clusters <- rev(parula(6)) # 1st colour is too bright --> adjust
pal.clusters[1] <- "#fbe30e"

# Re-name and re-order the indices
ddf2[ddf2$index == "SR","index"] <- "Richness"
ddf2[ddf2$index == "FR","index"] <- "Faith"
ddf2[ddf2$index == "Jac","index"] <- "Trait dissimilarity"
ddf2[ddf2$index == "Jtu","index"] <- "Trait turnover"
ddf2[ddf2$index == "Jne","index"] <- "Trait nestedness"
# unique(ddf2$index)

ddf2$index <- factor(ddf2$index, levels = c("Richness","Faith","FEve","FDis","FDiv","Trait dissimilarity","Trait turnover","Trait nestedness"))

plot <- ggplot(ddf2[which(ddf2$mean < 1),], aes(x = factor(region), y = mean*100, fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters, guide = "colourbar") +
    geom_hline(yintercept = 0) + geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("Average % difference (2100-2012)") + 
    facet_wrap(.~factor(index), scales = "free", nrow = 4)

setwd("/net/kryo/work/fabioben/BIOCeans5D/plots/")
ggsave(plot = plot, filename = "Fig.X_facet_deltas_fd_regions_k6_28.11.23.jpg", width = 7, height = 10, dpi = 300)
ggsave(plot = plot, filename = "Fig.X_facet_deltas_fd_regions_k6_28.11.23.pdf", width = 7, height = 10, dpi = 300)

### Pairwise Wilcoxon tests to test wif ∆indices differ significantly between regions
dat <- ddf2[which(ddf2$mean < 1 & ddf2$index == "Richness"),]
pairwise.wilcox.test(dat$mean, factor(dat$region), p.adjust.method = "bonferroni", paired = F)
### 1-4-5 do not differ

dat <- ddf2[which(ddf2$mean < 1 & ddf2$index == "Faith"),]
pairwise.wilcox.test(dat$mean, factor(dat$region), p.adjust.method = "bonferroni", paired = F)
### all pairs differ

dat <- ddf2[which(ddf2$mean < 1 & ddf2$index == "FEve"),]
pairwise.wilcox.test(dat$mean, factor(dat$region), p.adjust.method = "bonferroni", paired = F)
### all pairs differ

dat <- ddf2[which(ddf2$mean < 1 & ddf2$index == "FDis"),]
pairwise.wilcox.test(dat$mean, factor(dat$region), p.adjust.method = "bonferroni", paired = F)
### all pairs differ

dat <- ddf2[which(ddf2$mean < 1 & ddf2$index == "FDiv"),]
pairwise.wilcox.test(dat$mean, factor(dat$region), p.adjust.method = "bonferroni", paired = F)
### all pairs differ

dat <- ddf2[which(ddf2$mean < 1 & ddf2$index == "Trait dissimilarity"),]
pairwise.wilcox.test(dat$mean, factor(dat$region), p.adjust.method = "bonferroni", paired = F)
### all pairs differ

dat <- ddf2[which(ddf2$mean < 1 & ddf2$index == "Trait turnover"),]
pairwise.wilcox.test(dat$mean, factor(dat$region), p.adjust.method = "bonferroni", paired = F)
### 5 & 6 do not differ

dat <- ddf2[which(ddf2$mean < 1 & ddf2$index == "Trait nestedness"),]
pairwise.wilcox.test(dat$mean, factor(dat$region), p.adjust.method = "bonferroni", paired = F)
### 1 & 6 do not differ



### ------------------------------------------------------------------

### 10) Boxplots to support Fig 5: main differences across regions
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/clusters")
ras <- get(load("raster_clusters_ward_k6_28.11.23.RData"))
fd$region <- extract(ras, fd[,c("x","y")])
fd2 <- fd %>% drop_na(region)

# Define the colour palette for the regions
pal.clusters <- rev(parula(6)) 
pal.clusters[1] <- "#fbe30e"

# To make sure the clusters on Fig. 5 match the actual clusters in 'region'
# ggplot() + geom_raster(aes(x = x, y = y, fill = factor(region)), data = fd2) +
#     scale_fill_manual(name = "Regions", values = pal.clusters) +
#     geom_contour(colour = "black", binwidth = 1, size = .4, aes(x = x, y = y, z = region), data = fd2) +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

# ggplot(fd2, aes(x = factor(region), y = y, fill = factor(region))) +
#     scale_fill_manual(name = "Regions", values = pal.clusters, guide = "colourbar") +
#     geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) +
#     geom_hline(yintercept = 0, colour = "red") + xlab("Regions") + ylab("Latitude")
    

### Differences between regions 1-3 (tropicals) vs. 4-6 (higher latitude) are pretty obvious on the first order

### What does separate region #1 vs. #2 and #3 ?
dat.trop <- fd2[fd2$region %in% c(1:3),]
dat.extra <- fd2[fd2$region %in% c(4:6),]

### 10.A) Differences between 1-2-3 (tropical clusters)
# 4 shows higher FPOCex/POCflux and NPP
b1 <- ggplot(dat.extra, aes(x = factor(region), y = log10(NPP), fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(4:6)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("NPP\nlog10(mgC.m-2.d-1)") + theme_bw()

b2 <- ggplot(dat.extra, aes(x = factor(region), y = log10(FPOCex), fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(4:6)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("FPOCex\nlog10(mgC.m-2.d-1)") + theme_bw()

# 5 shows higher Faith, FDis and lower Jac/Jtu and productivity (CHLA)
b3 <- ggplot(dat.extra, aes(x = factor(region), y = Faith, fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(4:6)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("Faith index") + theme_bw()

b4 <- ggplot(dat.extra, aes(x = factor(region), y = FDis, fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(4:6)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("FDis") + theme_bw()
# 
b5 <- ggplot(dat.extra, aes(x = factor(region), y = Jac, fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(4:6)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("Trait dissimilarity\n(Jaccard)") + theme_bw()

b6 <- ggplot(dat.extra, aes(x = factor(region), y = log10(CHL), fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(4:6)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("Chlorophyll a\nlog10(mgC.m-3)") + theme_bw()

b7 <- ggplot(dat.extra, aes(x = factor(region), y = log10(MESOZOO), fill = factor(region))) +
        scale_fill_manual(name = "Regions", values = pal.clusters[c(4:6)], guide = "colourbar") +
        geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
        xlab("Regions") + ylab("Mesozooplankton biomass\nlog10(mgC.m-3)") + theme_bw()

# Assemble in panel and save
panel <- ggarrange(b3,b4,b5,b7,b2,b6,b1, ncol = 2, nrow = 4, align = "hv", labels = letters)

setwd("/net/kryo/work/fabioben/BIOCeans5D/plots/")
ggsave(plot = panel, filename = "Fig.X_panel_boxplots_tropical_regions_29.11.23.jpg", height = 12, width = 7, dpi = 300)


### 10.B) Differences between 4-5-6
# 2 stands out from 1 and 3 because higher diversity (SR, Faith) and lower NPP/ FPOCex/POCflux
b1 <- ggplot(dat.trop, aes(x = factor(region), y = SR, fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(1:3)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("Species richness") + theme_bw()

b2 <- ggplot(dat.trop, aes(x = factor(region), y = Faith, fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(1:3)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("Faith index") + theme_bw()   

b3 <- ggplot(dat.trop, aes(x = factor(region), y = log10(MESOZOO), fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(1:3)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("Mesozooplankton biomass\nlog10(mgC.m-3)") + theme_bw()

# 3 stands out from 1 because of FEve and POC fluxes and CHLA
b4 <- ggplot(dat.trop, aes(x = factor(region), y = log10(FPOCex), fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(1:3)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("FPOCex\nlog10(mgC.m-2.d-1)") + theme_bw()

b5 <- ggplot(dat.trop, aes(x = factor(region), y = FEve, fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(1:3)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("FEve") + theme_bw()

b6 <- ggplot(dat.trop, aes(x = factor(region), y = log10(CHL), fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(1:3)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("Chlorophyll a\nlog10(mgC.m-3)") + theme_bw()

b7 <- ggplot(dat.trop, aes(x = factor(region), y = log10(PROCHL), fill = factor(region))) +
    scale_fill_manual(name = "Regions", values = pal.clusters[c(1:3)], guide = "colourbar") +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", width = .2) + 
    xlab("Regions") + ylab("Prochlorococcus biomass\nlog10(mgC.m-3)") + theme_bw()

# Assemble in panel and save
panel <- ggarrange(b1,b2,b5,b3,b6,b7,b4, ncol = 2, nrow = 4, align = "hv", labels = letters[8:15])
ggsave(plot = panel, filename = "Fig.X_panel_boxplots_extratropical_regions_29.11.23.jpg", height = 12, width = 7, dpi = 300)


### ------------------------------------------------------------------

### 11) IntDi distribution & variations across FGs from Benedetti et al. (2023); + relationship with shift metrics
setwd("/net/kryo/work/fabioben/BIOCeans5D/data")

### 11.1) Examine Intdi distribution per copepod species and per FG
traits <- read.csv("traits_table_Benedetti2023.csv", h = T, sep = ";", dec = ",")
colnames(traits)[8] <- "Body.length" 
require("naniar")
traits <- traits %>% replace_with_na_all(condition = ~.x == "")
traits$Spawning.mode <- as.factor(traits$Spawning.mode)
traits$Trophic.group <- as.factor(traits$Trophic.group)
traits$Feeding.mode <- as.factor(traits$Feeding.mode)
names <- colnames(traits)[c(8:11,16)] ; names
traits$na_count <- apply(traits[,names], 1, function(x) sum(is.na(x)))
traits_red <- traits[!is.na(traits$Body.length),]
traits_red2 <- traits_red[traits_red$na_count < 2,]

traits_red2$Myelination <- as.integer(as.logical(traits_red2$Myelination))
traits_red2$Omnivore <- as.integer(as.logical(traits_red2$Omnivore))
traits_red2$Carnivore <- as.integer(as.logical(traits_red2$Carnivore))
traits_red2$Herbivore <- as.integer(as.logical(traits_red2$Herbivore))
traits_red2$Detritivore <- as.integer(as.logical(traits_red2$Detritivore))
traits_red2$Current <- as.integer(as.logical(traits_red2$Current))
traits_red2$Cruise <- as.integer(as.logical(traits_red2$Cruise))
traits_red2$Ambush <- as.integer(as.logical(traits_red2$Ambush))

traits_red2 <- data.frame(traits_red2)
# dim(traits_red2)

# Compute Gower distance matrix
gow <- gawdis(x = traits_red2[,c(8,9,10,12:15,17:19)], groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))

# Load list of IntDi ranks per species
list_mat <- get(load("list_gawdis_mat_IntDi_24_07_23.Rdata"))

# Compute average matrix
mean_matrix_gaw <- Reduce("+",list_mat)/length(list_mat) 
# Compute average IntDi values
IntDi <- apply(mean_matrix_gaw, 1, sum, na.rm = T) / (nrow(mean_matrix_gaw) + 1)
IntDi <- cbind(as.data.frame(IntDi), traits_red2$Species)
colnames(IntDi)[2] <- "Species"
IntDi[order(IntDi$IntDi, decreasing = T),]
# IntDi[IntDi$IntDi < .2,] # Why 0? 
# Fix values == 0 ?
# IntDi[IntDi$IntDi == 0,"IntDi"] <- .1

### Plot
IntDi$Species <- factor(IntDi$Species, levels = IntDi[order(IntDi$IntDi),"Species"]) # 
a <- ggplot(IntDi[IntDi$IntDi > 0,], aes(x = Species, y = IntDi)) + geom_point() + ylab("Integrated distinctiveness (IntDi)") +
    scale_x_discrete(name = "Species", labels = NULL) + theme_classic() + theme(axis.ticks = element_blank()) 

#setwd("/net/kryo/work/fabioben/GODLY/plots") # working dir
#ggsave(plot = a, filename = "Fig.X_IntDi_distrib_29.11.23.jpg", height = 3, width = 4)

### Plot distribution of IntDi per FGs of Benedetti et al. 2022
traits_red2$IntDi <- IntDi$IntDi

b <- ggplot(aes(y = IntDi, x = factor(FG)), data = traits_red2[traits_red2$IntDi > 0 & !is.na(traits_red2$FG),]) +
    geom_boxplot(colour = "black", fill = "gray") +  ylab("Integrated distinctiveness (IntDi)") +
    xlab("Copepod FG (Benedetti et al., 2022)") + theme_classic()

panel <- ggarrange(a,b, ncol = 1, nrow = 2, align = "hv", labels = letters)
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots") # working dir
ggsave(plot = panel, filename = "Fig.X_panel_IntDi_distrib_FGs_29.11.23.jpg", height = 7, width = 5)

### Test 
#pairwise.wilcox.test(traits_red2[traits_red2$IntDi > 0 & !is.na(traits_red2$FG),"IntDi"], traits_red2[traits_red2$IntDi > 0 & !is.na(traits_red2$FG),"FG"], p.adjust.method = "bonferroni", paired = F)
#    1       2       3       4       5       6       7       8       9      
# 2  0.00024 -       -       -       -       -       -       -       -
# 3  0.09376 2.4e-07 -       -       -       -       -       -       -
# 4  0.00098 7.5e-08 1.00000 -       -       -       -       -       -
# 5  0.01097 1.00000 7.2e-07 2.2e-06 -       -       -       -       -
# 6  1.00000 0.54331 1.00000 1.00000 0.30305 -       -       -       -
# 7  6.1e-07 6.6e-10 1.5e-08 5.2e-12 8.6e-12 2.5e-05 -       -       -
# 8  1.6e-06 1.3e-09 1.7e-05 3.4e-10 6.2e-11 0.00072 1.00000 -       -
# 9  1.9e-06 5.5e-09 9.4e-14 1.0e-10 1.0e-10 4.4e-06 1.1e-14 1.7e-13 -
# 10 0.00041 3.4e-05 1.6e-06 8.4e-06 8.4e-06 0.00061 1.4e-06 0.18371 3.3e-06
# 11 5.2e-06 3.3e-08 6.9e-12 1.2e-09 1.2e-09 1.1e-05 3.4e-11 0.00015 1.00000
#    10
# 2  -
# 3  -
# 4  -
# 5  -
# 6  -
# 7  -
# 8  -
# 9  -
# 10 -
# 11 1.00000

### Plotting PCoA space
library("vegan")
pcoa2 <- wcmdscale(d = gow, eig = T)
pcoa.scores <- data.frame(pcoa2$points)
scores1 <- paste0("PCoA 1 (",floor(pcoa2$eig[1]*100)/100,"%)")
scores2 <- paste0("PCoA 2 (",floor(pcoa2$eig[2]*100)/100,"%)")
scores3 <- paste0("PCoA 3 (",floor(pcoa2$eig[3]*100)/100,"%)")
scores4 <- paste0("PCoA 4 (",floor(pcoa2$eig[4]*100)/100,"%)")

# dim(pcoa.scores)
IntDi$PCoA1 <- pcoa.scores$Dim1
IntDi$PCoA2 <- pcoa.scores$Dim2
IntDi$PCoA3 <- pcoa.scores$Dim3
IntDi$PCoA4 <- pcoa.scores$Dim4

# round(cor(IntDi[IntDi$IntDi > 0,"IntDi"], IntDi[IntDi$IntDi > 0,"PCoA1"], method = "spearman"),3) # -0.279
# round(cor(IntDi[IntDi$IntDi > 0,"IntDi"], IntDi[IntDi$IntDi > 0,"PCoA2"], method = "spearman"),3) # 0.669
# round(cor(IntDi[IntDi$IntDi > 0,"IntDi"], IntDi[IntDi$IntDi > 0,"PCoA3"], method = "spearman"),3) # -0.053 N.S.
# round(cor(IntDi[IntDi$IntDi > 0,"IntDi"], IntDi[IntDi$IntDi > 0,"PCoA4"], method = "spearman"),3) # 0.261
# PCoA 2 seems to be the main axis of variation

plot <- ggplot(data=IntDi[IntDi$IntDi > 0,], aes(x = PCoA2, y = IntDi)) +
  geom_point(colour = "grey30") + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1,
      label.y = "top", label.x = "left", size = 3, colour = "black") + 
  xlab("PCoA 2") + ylab("IntDi") + theme_bw() 


# Plot in PCoA space
p1 <- ggplot(aes(x = PCoA1, y = PCoA2), data = IntDi[IntDi$IntDi > 0,]) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..)), fill = "grey") + 
    geom_point(aes(fill = IntDi), pch = 21, colour = "black") +
    scale_fill_gradientn(name = "IntDi", colours = parula(100), guide = "colourbar") +
    scale_alpha_continuous(range = c(0, 1)) + 
    theme_bw() + xlab(scores1) + ylab(scores2)+ guides(alpha = "none")

p2 <- ggplot(aes(x = PCoA3, y = PCoA4), data = IntDi[IntDi$IntDi > 0,]) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..)), fill = "grey") + 
    geom_point(aes(fill = IntDi), pch = 21, colour = "black") + 
    scale_fill_gradientn(name = "IntDi", colours = parula(100), guide = "colourbar") +
    scale_alpha_continuous(range = c(0, 1)) + 
    theme_bw() + xlab(scores3) + ylab(scores4)+ guides(alpha = "none")

pcoa.plot <- ggarrange(p1,p2, labels = letters, align = "hv", ncol = 2, nrow = 1)
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots/")
ggsave(plot = pcoa.plot, filename = "Fig.X_PCoA+IntDi_01.12.23.jpg", dpi = 300, height = 4, width = 10)
ggsave(plot = pcoa.plot, filename = "Fig.X_PCoA+IntDi_01.12.23.pdf", dpi = 300, height = 4, width = 10)
ggsave(plot = plot, filename = "Fig.X_IntDixPCoA.2_01.12.23.jpg", dpi = 300, height = 3, width = 4)


### 11.2) Examine covariance between IntDi index and intensity of shift metrics
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/shift_metrics")
res <- lapply(dir()[grep("contemp",dir())], function(f) { d <- get(load(f)) ; return(d) })
base <- dplyr::bind_rows(res); rm(res); gc()

# For the 'global' centroids, compute each future shift relative to the SDM-specific contemporary positions
files <- dir()[grepl("fut",dir())]
res <- lapply(files, function(f) {
        
        message(paste("\n","Compute shifts in mean annual centroids for ",f,"\n", sep = ""))
        d <- get(load(f)) 
        # Choose reference of interest wihtin 'base'
        centroids.base <- base[which(base$SDM == unique(d$SDM) & base$region == "Global"),]
        
        # For each species, compute spatial shifts
        spp <- unique(d$species) 
        # s <- "Acartia_Acartiura_clausi" # for testing
        res2 <- mclapply(spp, function(s) {
                    message(paste("Computing shifts for ",s, sep = ""))
                    require("geosphere")
                    ref <- centroids.base[centroids.base$species == s,]
                    fut <- d[d$species == s & d$region == "Global",]
                    dist <- distm(ref[,c("lon","lat")], fut[,c("lon","lat")], fun = distHaversine)
                    # Convert to km
                    dist.km <- dist/1000
                    return(data.frame(species = s, shift.km = dist.km))
            }, mc.cores = 20
                
        ) # eo mclapply
        shifts <- dplyr::bind_rows(res2)
        rm(res2); gc()
        
        shifts$SDM <- unique(d$SDM)
        shifts$ESM <- unique(d$ESM)
        
        # Return
        return(shifts)
    }
         
) # eo lapply
shifts <- dplyr::bind_rows(res)
rm(res); gc()

# Summarize ensemble prediction
ens <- data.frame( shifts %>% group_by(species) %>% summarise(mean.shift = mean(shift.km), sd = sd(shift.km)) )
ens$species <- as.character(ens$species)
# Supply to 'traits_red2' and analyze covariance with IntDi
traits_red2$mean.shift <- NA
traits_red2$sd.shift <- NA

for(s in unique(ens$species)) {
    message(s)
    traits_red2[traits_red2$Species == s,"mean.shift"] <- ens[ens$species == s,"mean.shift"]
    traits_red2[traits_red2$Species == s,"sd.shift"] <- ens[ens$species == s,"sd"]
} # eo for loop - s in species

# Test covariance between IntDi and mean shifts amplitude
ggplot(data = traits_red2[traits_red2$IntDi > 0,], aes(x = IntDi, y = mean.shift)) + geom_point() + theme_bw()
df <- traits_red2 %>% filter(!(is.na(mean.shift)))
# cor(df$IntDi, df$mean.shift, method = "spearman") # -0.074
# summary(lm(IntDi ~ mean.shift, data = df)) # Adjusted R-squared: 0.03 
# No relationship

# No signif. trend; Plot distrib of shifts per FG (boxplots)
ggplot(aes(y = mean.shift, x = factor(FG)), data = traits_red2[traits_red2$IntDi > 0 & !is.na(traits_red2$FG),]) +
    geom_boxplot(colour = "black", fill = "gray") + ylab("Centroid shift (km)") + xlab("FG (Benedetti et al., 2023)") +
    theme_classic()
### --> all medians kind of overlap with IQR; no strong variance in shifts per known FG

# pairwise.wilcox.test(traits_red2[traits_red2$IntDi > 0 & !is.na(traits_red2$FG),"mean.shift"],
#     traits_red2[traits_red2$IntDi > 0 & !is.na(traits_red2$FG),"FG"],
#     p.adjust.method = "bonferroni", paired = F)
# # n.s. 


### 11.3) Same as above, but based on changes in HSI. Examine Global/NH/SH separately

setwd("/net/kryo/work/fabioben/BIOCeans5D/data/shift_metrics")
res <- lapply(dir()[grep("contemp_mean_ann_hsi",dir())], function(f) { d <- get(load(f)) ; return(d) })
base <- dplyr::bind_rows(res); rm(res); gc()

res <- lapply(dir()[grepl("fut_mean_ann_hsi",dir())], function(f) {
        
        message(paste("\n","Compute shifts in mean HSI for ",f,"\n", sep = ""))
        d <- get(load(f)) 
        # Choose reference of interest wihtin 'base'
        hsi.base <- base[which(base$SDM == unique(d$SDM) & base$region == "SH"),]
        hsi.fut <- d[which(d$region == "NH"),]
        # dim(hsi.fut); dim(hsi.base)
        
        changes <- data.frame(species = hsi.base$species, HSI.base = hsi.base$mean, HSI.fut = hsi.fut$mean)
        changes$perc <- ((changes$HSI.fut-changes$HSI.base)/changes$HSI.base)*100
        
        changes$SDM <- unique(d$SDM)
        changes$ESM <- unique(d$ESM)
        
        # Return
        return(changes)
    }
         
) # eo lapply
changes <- dplyr::bind_rows(res)
rm(res); gc()
# summary(changes); dim(changes)

# Summarize ensemble prediction
ens.changes <- data.frame( changes %>% group_by(species) %>% summarise(mean.change = mean(perc), sd = sd(perc)) )
ens.changes$species <- as.character(ens.changes$species)
#  summary(ens.changes)

# Supply to 'traits_red2' 
traits_red2$mean.ch <- NA
traits_red2$sd.ch <- NA
for(s in unique(ens.changes$species)) {
    message(s)
    traits_red2[traits_red2$Species == s,"mean.ch"] <- ens.changes[ens.changes$species == s,"mean.change"]
    traits_red2[traits_red2$Species == s,"sd.ch"] <- ens.changes[ens.changes$species == s,"sd"]
} # eo for loop - s in species

# Examine covariance
ggplot(data = traits_red2[traits_red2$IntDi > 0,], aes(x = IntDi, y = mean.ch)) + geom_point() + theme_bw()
df <- traits_red2 %>% filter(!(is.na(mean.ch)))
#cor(df$IntDi, df$mean.ch, method = "spearman") # -0.03 Globally; -0.04 for the NH only; -0.06 for the SH only
#summary(lm(IntDi ~ mean.ch, data = df)) # Adjusted R-squared: 0.008; N.S.
# No relationship

# No signif. trend; Plot distrib of shifts per FG (boxplots)
ggplot(aes(y = mean.ch, x = factor(FG)), data = traits_red2[traits_red2$IntDi > 0 & !is.na(traits_red2$FG),]) +
    geom_boxplot(colour = "black", fill = "gray") + ylab("Change in relative HSI (%)") + xlab("FG (Benedetti et al., 2023)") +
    theme_classic()
### --> all medians kind of overlap with IQR; no strong variance in shifts per known FG

# pairwise.wilcox.test(traits_red2[traits_red2$IntDi > 0 & !is.na(traits_red2$FG),"mean.ch"],
#     traits_red2[traits_red2$IntDi > 0 & !is.na(traits_red2$FG),"FG"],
#     p.adjust.method = "bonferroni", paired = F)
# n.s. 


### ------------------------------------------------------------------

### 12) Co-ranking AUC criterion ~ N PCoA dimensions to estumate FD indices

library("dimRed")
library("coRanking")
library("vegan")
library("gawdis")
library("FD")
library("naniar")

setwd("/net/kryo/work/fabioben/BIOCeans5D/data") 

traits <- read.csv("traits_table_Benedetti2023.csv", h = T, sep = ";", dec = ",")
colnames(traits)[8] <- "Body.length" 
traits <- traits %>% replace_with_na_all(condition = ~.x == "")
traits$Spawning.mode <- as.factor(traits$Spawning.mode)
traits$Trophic.group <- as.factor(traits$Trophic.group)
traits$Feeding.mode <- as.factor(traits$Feeding.mode)
names <- colnames(traits)[c(8:11,16)] ; names
traits$na_count <- apply(traits[,names], 1, function(x) sum(is.na(x)))
traits_red <- traits[!is.na(traits$Body.length),]
traits_red2 <- traits_red[traits_red$na_count < 2,]
traits_red2$Myelination <- as.integer(as.logical(traits_red2$Myelination))
traits_red2$Omnivore <- as.integer(as.logical(traits_red2$Omnivore))
traits_red2$Carnivore <- as.integer(as.logical(traits_red2$Carnivore))
traits_red2$Herbivore <- as.integer(as.logical(traits_red2$Herbivore))
traits_red2$Detritivore <- as.integer(as.logical(traits_red2$Detritivore))
traits_red2$Current <- as.integer(as.logical(traits_red2$Current))
traits_red2$Cruise <- as.integer(as.logical(traits_red2$Cruise))
traits_red2$Ambush <- as.integer(as.logical(traits_red2$Ambush))
traits_red2 <- data.frame(traits_red2)
# dim(traits_red2)

# Compute Gower's distance matrix with fuzzy coding
gow <- gawdis(x = traits_red2[,c(8,9,10,12:15,17:19)], groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))

# Without Cailliez correction
pcoa1 <- wcmdscale(d = gow, eig = T)
pcoa1.scores <- data.frame(pcoa1$points) # dim(pcoa.scores) # 124 dimensions
# With Cailliez correction
pcoa2 <- ape::pcoa(gow, correction = "cailliez")

# dim(pcoa.scores); dim(pcoa2$vectors)

res <- mclapply(c(1:20), function(n) {
    
        # Message
        message(paste("Testing ",n," PCoA dimensions", sep = ""))
        euclid_dist1 <- dist(pcoa1.scores[,1:n], method = "euclidean")
        euclid_dist2 <- dist(pcoa2$vectors[,1:n], method = "euclidean")
        
        # Compute co-ranking scores for non Cailliez-corrected version
        Q1 <- coranking(Xi = gow, X = euclid_dist1, input_Xi = c("dist"), input_X = c("dist"), use = "R")
        NX1 <- coRanking::R_NX(Q1)
        AUC1 <- coRanking::AUC_ln_K(NX1)
        
        # Same but for Cailliez-corrected version
        Q2 <- coranking(Xi = gow, X = euclid_dist2, input_Xi = c("dist"), input_X = c("dist"), use = "R")
        NX2 <- coRanking::R_NX(Q2)
        AUC2 <- coRanking::AUC_ln_K(NX2)
        
        # Retrun
        return(data.frame(n = n, AUC1 = round(AUC1,3), AUC2 = round(AUC2,3)))
        
    }, mc.cores = 15
    
) # eo mclapply
# Rbind
table <- bind_rows(res)

# Plot
plot <- ggplot(data = table, aes(x = n, y = AUC1)) + geom_point(colour = "black") + 
  geom_vline(xintercept = 4, colour = "red") + xlab("Number of PCoA dimensions") + 
  ylab("Quality of traits space (AUC)") + theme_bw()

setwd("/net/kryo/work/fabioben/BIOCeans5D/plots") # working dir
ggsave(plot = plot, filename = "Fig.X_AUC_PCoA_dimensions_29.11.23.jpg", height = 3, width = 4)



### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 01/12/23: 13) Comparing FD outputs between HSI-based communities and threshold-based PA communities
### Re-using Script#4.1

### 13.1) Comparing mean annual Faith x mean annual FRic (standardized this time) - both indices base don PA data 
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/Faith/Faith")
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
tab$cell_id <- factor(paste(tab$x,tab$y,sep = "_")) # length(unique(tab$cell_id))
ann.faith <- data.frame(tab %>% group_by(cell_id) %>% 
        summarize(x = unique(x), y = unique(y),
            mean.faith = mean(Faith, na.rm = T), sd.faith = sd(Faith, na.rm = T),
            mean.rich = mean(SR, na.rm = T), sd.rich = sd(SR, na.rm = T))
)
# dim(ann.faith) ; summary(ann.faith)

setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/PA_Gawdis_PCoA_Euclid/Stand.Fric") 
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
# dim(ann.ind.pa); summary(ann.ind.pa) 

# ggplot() + geom_raster(aes(x = x, y = y, fill = FRic.avg), data = ann.ind.pa) +
#      scale_fill_gradientn(name = "FRic", colours = parula(100), guide = "colourbar") +
#      geom_contour(colour = "grey30", binwidth = .05, size = 0.25, aes(x = x, y = y, z = FRic.avg), data = ann.ind.pa) +
#      geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#      coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#      panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#      scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#      scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

ann.faith$FRic <- ann.ind.pa$FRic.avg
#cor(ann.faith$mean.faith, ann.faith$FRic, method = "spearman") # 0.772

# Standardize them by their max values; won't change the correlation coeff of course
ann.faith$Faith.std <- (ann.faith$mean.faith)/max(ann.faith$mean.faith, na.rm = T)
ann.faith$FRic.std <- (ann.faith$FRic)/max(ann.faith$FRic, na.rm = T) 
# Bi-plot
library("ggpmisc")
formula1 <- y ~ x 
formula2 <- y ~ poly(x,2)

p <- ggplot(ann.faith, aes(x = mean.faith, y = FRic, colour = abs(y))) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "bottom", label.x = "right", size = 3, colour = "black") + 
  xlab("Faith index") + ylab("FRcic (standardize)") + theme_bw() 

setwd("/net/kryo/work/fabioben/BIOCeans5D/plots/")
ggsave(plot = p, filename = "Fig.X_FRicxFaith_PA_01.12.23.pdf", dpi = 300, width = 5, height = 4)
ggsave(plot = p, filename = "Fig.X_FRicxFaith_PA_01.12.23.jpg", dpi = 300, width = 5, height = 4)

    
### Per biome (just a curiousity)
# Tropics
# cor(ann.faith[which(abs(ann.faith$y) < 30),"Faith.std"], ann.faith[which(abs(ann.faith$y) < 30),"FRic.std"], method = "spearman") # 0.94
# summary(lm(FRic.std ~ Faith.std, data = ann.faith[which(abs(ann.faith$y) < 30),])) # R-squared: 0.90
# # High latitudes
# cor(ann.faith[which(abs(ann.faith$y) > 60),"Faith.std"], ann.faith[which(abs(ann.faith$y) > 60),"FRic.std"], method = "spearman") # 0.875
# summary(lm(FRic.std ~ Faith.std, data = ann.faith[which(abs(ann.faith$y) > 60),])) # R-squared: 0.75
# # Temperate latitudes
# cor(ann.faith[which(abs(ann.faith$y) < 60 & abs(ann.faith$y) > 30),"Faith.std"], ann.faith[which(abs(ann.faith$y) < 60 & abs(ann.faith$y) > 30),"FRic.std"], method = "spearman")
# # 0.76
# summary(lm(FRic.std ~ Faith.std, data = ann.faith[which(abs(ann.faith$y) < 60 & abs(ann.faith$y) > 30),])) # R-squared: 0.628
### --> Departure from correlation outside the tropics because the GLM_bias in the high latitude FRic


### 13.2) PA-based FD indices (FEve/RaoQ/FDis/FDiv) vs. HSI-based indices

setwd("/net/kryo/work/fabioben/BIOCeans5D/data/fd_indices/db_FD/HSI_Gawdis_PCoA_Euclid")
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
            FDiv.hsi = mean(FDiv, na.rm = T), FDiv.std = sd(FDiv, na.rm = T)
        )
) # eo summarise 
# dim(ann.ind.hsi)
ann.ind.hsi$FRic.pa <- ann.ind.pa$FRic.avg
ann.ind.hsi$FEve.pa <- ann.ind.pa$FEve.avg
ann.ind.hsi$FDis.pa <- ann.ind.pa$FDis.avg
ann.ind.hsi$RaoQ.pa <- ann.ind.pa$RaoQ.avg
ann.ind.hsi$FDiv.pa <- ann.ind.pa$FDiv.avg

### Compute correlations coeff
### NB: FRic based on HSI doe snot make sense at all because FRic only sees PA data and any HSI > 0 is automatically converted to 1 by dbFD() 

#cor(ann.ind.hsi$FEve.pa, ann.ind.hsi$FEve.hsi, method = "spearman") # 0.61
#summary(lm(FEve.pa ~ FEve.hsi, data = ann.ind.hsi)) # R-squared:  0.298

#cor(ann.ind.hsi$FDis.pa, ann.ind.hsi$FDis.hsi, method = "spearman") # 0.87
#summary(lm(FDis.pa ~ FDis.hsi, data = ann.ind.hsi)) # R-squared:  0.655

#cor(ann.ind.hsi$RaoQ.pa, ann.ind.hsi$RaoQ.hsi, method = "spearman") # 0.77
#summary(lm(RaoQ.pa ~ RaoQ.hsi, data = ann.ind.hsi)) # R-squared:  0.579

#cor(ann.ind.hsi$FDiv.pa, ann.ind.hsi$FDiv.hsi, method = "spearman") # 0.28
#summary(lm(FDiv.pa ~ FDiv.hsi, data = ann.ind.hsi)) # R-squared:  0.104

### Draw biplots
p1 <- ggplot(ann.ind.hsi, aes(x = FEve.hsi, y = FEve.pa, colour = abs(y))) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  #geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "bottom", label.x = "right", size = 3, colour = "black") + 
  xlab("FEve (based on HSI)") + ylab("FEve (based on 1/0)") + theme_bw() 

p2 <- ggplot(ann.ind.hsi, aes(x = FDis.hsi, y = FDis.pa, colour = abs(y))) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  #geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "bottom", label.x = "right", size = 3, colour = "black") + 
  xlab("FDis (based on HSI)") + ylab("FDis (based on 1/0)") + theme_bw() 

p3 <- ggplot(ann.ind.hsi, aes(x = RaoQ.hsi, y = RaoQ.pa, colour = abs(y))) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  #geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "bottom", label.x = "right", size = 3, colour = "black") + 
  xlab("Rao's Q (based on HSI)") + ylab("Rao's Q (based on 1/0)") + theme_bw() 

p4 <- ggplot(ann.ind.hsi, aes(x = FDiv.hsi, y = FDiv.pa, colour = abs(y))) +
  geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  #geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "bottom", label.x = "right", size = 3, colour = "black") + 
  xlab("FDiv (based on HSI)") + ylab("FDiv (based on 1/0)") + theme_bw() 

panel <- ggarrange(p1,p2,p3,p4, align = "hv", ncol = 2, nrow = 2, labels = letters)
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
ggsave(plot = panel, filename = "Fig.X_panel_fd_HSIvsPA_01.12.23.jpg", dpi = 200, width = 11, height = 9)


### Maps to check
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FEve.pa), data = ann.ind.hsi) +
    scale_fill_gradientn(name = "FEve (1/0)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "grey30", binwidth = .025, size = 0.25, aes(x = x, y = y, z = FEve.pa), data = ann.ind.hsi) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDis.pa), data = ann.ind.hsi) +
    scale_fill_gradientn(name = "FDis (1/0)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = FDis.pa), data = ann.ind.hsi) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map3 <- ggplot() + geom_raster(aes(x = x, y = y, fill = RaoQ.pa), data = ann.ind.hsi) +
    scale_fill_gradientn(name = "Rao's Q (1/0)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "grey30", binwidth = .02, size = 0.25, aes(x = x, y = y, z = RaoQ.pa), data = ann.ind.hsi) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

map4 <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDiv.pa), data = ann.ind.hsi) +
    scale_fill_gradientn(name = "FDiv (1/0)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "grey30", binwidth = .01, size = 0.25, aes(x = x, y = y, z = FDiv.pa), data = ann.ind.hsi) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

panel <- ggarrange(map1,map2,map3,map4, align = "hv", ncol = 2, nrow = 2, labels = letters)
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
ggsave(plot = panel, filename = "Fig.X_maps_fd_PA_01.12.23.jpg", dpi = 200, width = 10, height = 4.5)


# Check FEve without and with tropics
#cor(ann.ind.hsi[which(abs(ann.ind.hsi$y) < 20),"FEve.hsi"], ann.ind.hsi[which(abs(ann.ind.hsi$y) < 20),"FEve.pa"], method = "spearman") # 0.46
#cor(ann.ind.hsi[which(abs(ann.ind.hsi$y) < 30),"FEve.hsi"], ann.ind.hsi[which(abs(ann.ind.hsi$y) < 30),"FEve.pa"], method = "spearman") # 0.47
#cor(ann.ind.hsi[which(abs(ann.ind.hsi$y) < 50),"FEve.hsi"], ann.ind.hsi[which(abs(ann.ind.hsi$y) < 50),"FEve.pa"], method = "spearman") # 0.65
#cor(ann.ind.hsi[which(abs(ann.ind.hsi$y) < 60),"FEve.hsi"], ann.ind.hsi[which(abs(ann.ind.hsi$y) < 60),"FEve.pa"], method = "spearman") # 0.67
#cor(ann.ind.hsi[which(abs(ann.ind.hsi$y) < 70),"FEve.hsi"], ann.ind.hsi[which(abs(ann.ind.hsi$y) < 70),"FEve.pa"], method = "spearman") # 0.64
# So, better if you remove latitudes lower than 60° but not that much

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
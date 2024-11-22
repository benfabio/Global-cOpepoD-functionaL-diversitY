
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 23/11/23: R script to make the final main figures for the manuscript in preparation for Global Change Biology © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Figures and panels to make:
# - Fig.1: Panel of mean annual 'alpha' diversity indices (contemporary conditions) - SR/Faith/SES Faith/FEve/FDis or Rao's Q/FDiv
# - Fig.2: Panel of mean annual 'beta' diversity indices (contemporary conditions) - Jac/Jtu/Jne (+ beta ratio?) - make 2 versions
# - Fig.3: Panel of bivariate plots between FD indices and SR 
# - Fig.4: Emergent global B-EF relationships (most interesting ones) - chose MESOZOOPLANK BIOMASS 
# - Fig.5: Clustergram analysis + map of emergent regions
# - Fig.6: Panel of changes in mean annual diversity indices - to chose - SR/Faith/FEve/FDis or Rao's Q/FDiv if changes are strong/ beta div?

### Latest update: 19/11/24 (revising some figures for Minor Revisions)

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
library("ecodist")
library("raster")

world <- map_data("world") 
world2 <- map_data("world2")

setwd("/net/kryo/work/fabioben/BIOCeans5D/data")

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 1°) Making the main figures 

### A) Fig.1: Panel of mean annual diverity indices
fd <- get(load("table_mean_ann_FD_indices_baseline+BCP+biom_22.11.23.RData"))
# dim(fd); colnames(fd)

sr <- ggplot() + geom_raster(aes(x = x, y = y, fill = SR), data = fd) +
    scale_fill_gradientn(name = "Species\nrichness",colours = parula(100), guide = "colourbar") +
    #scale_fill_distiller(name = "Species richness", values = parula(100)) +
    geom_contour(colour = "black", binwidth = 20, size = .4, aes(x = x, y = y, z = SR), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

faith <- ggplot() + geom_raster(aes(x = x, y = y, fill = Faith), data = fd) +
    scale_fill_gradientn(name = "Faith\nindex",colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .25, size = .4, aes(x = x, y = y, z = Faith), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

ses.faith <- ggplot() + geom_raster(aes(x = x, y = y, fill = SES.Faith), data = fd) +
    scale_fill_gradient2(name = "SES Faith\nindex", guide = "colourbar", low = "#3D50C3", mid = "white", high = "#f50735") +
    geom_contour(colour = "black", binwidth = 1, size = .4, aes(x = x, y = y, z = SES.Faith), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

feve <- ggplot() + geom_raster(aes(x = x, y = y, fill = FEve), data = fd) +
    scale_fill_gradientn(name = "FEve",colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .08, size = .4, aes(x = x, y = y, z = FEve), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

fdis <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDis), data = fd) +
    scale_fill_gradientn(name = "FDis",colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .01, size = .4, aes(x = x, y = y, z = FDis), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

rao <- ggplot() + geom_raster(aes(x = x, y = y, fill = RaoQ.scaled), data = fd) +
    scale_fill_gradientn(name = "Rao's Q",colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .03, size = .4, aes(x = x, y = y, z = RaoQ.scaled), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
#
fdiv <- ggplot() + geom_raster(aes(x = x, y = y, fill = FDiv), data = fd) +
    scale_fill_gradientn(name = "FDiv",colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .01, size = .4, aes(x = x, y = y, z = FDiv), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

### Use FDis instead of Rao's Q in panel
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
fig1 <- ggarrange(sr,faith,ses.faith,feve,fdis,fdiv, align = 'hv', ncol = 2, nrow = 3, labels = letters)
ggsave(plot = fig1, filename = "Fig.1_mean_ann_fd_indices_12.12.23.jpg", dpi = 300, width = 10, height = 6.5)
ggsave(plot = fig1, filename = "Fig.1_mean_ann_fd_indices_12.12.23.pdf", dpi = 300, width = 10, height = 6.5)

# Also save individual maps 
ggsave(plot = sr, filename = "map_mean_ann_sr_23.11.23.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = faith, filename = "map_mean_ann_faith_23.11.23.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = ses.faith, filename = "map_mean_ann_ses.faith_23.11.23.jpg", dpi = 300, height = 4, width = 7) # Need to plot the p-values too here --> SM
ggsave(plot = feve, filename = "map_mean_ann_feve_23.11.23.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = fdis, filename = "map_mean_ann_fdis_23.11.23.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = rao, filename = "map_mean_ann_rao_23.11.23.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = fdiv, filename = "map_mean_ann_fdiv_23.11.23.jpg", dpi = 300, height = 4, width = 7)


### 13/12/23: Writing up the results section
# cor(fd[which(abs(fd$y) < 30),"SR"], fd[which(abs(fd$y) < 30),"Faith"], method = "spearman") # 0.94
# cor(fd[which(abs(fd$y) > 30),"SR"], fd[which(abs(fd$y) > 30),"Faith"], method = "spearman") # 0.87
# cor(fd[which(abs(fd$y) > 45),"SR"], fd[which(abs(fd$y) > 45),"Faith"], method = "spearman") # 0.87
# summary(fd[which(abs(fd$y) < 30),"SR"])
# summary(fd[which(abs(fd$y) > 30),"SR"])
#
# summary(fd[,"SES.Faith"])
# sd(fd[,"SES.Faith"])
#
# summary(fd[which(fd$SES.Faith > 0),"SES.Faith"])
# summary(fd[which(fd$SES.Faith < 0),"SES.Faith"])
#
# summary(fd[which(abs(fd$y) < 30),"SES.Faith"])
# summary(fd[which(abs(fd$y) > 31),"SES.Faith"])
#
# unique(abs(fd[which(fd$SES.Faith > 0),"y"]))
# nrow(fd[which(fd$SES.Faith > 0),])
# nrow(fd[which(fd$SES.Faith < 0),])
# nrow(fd[which(fd$SES.Faith > 0),])/(nrow(fd))
#
# # summary(fd[,"FEve"]); sd(fd[,"FEve"])
# mean(fd[which(abs(fd$y) < 30),"FEve"]); sd(fd[which(abs(fd$y) < 30),"FEve"])
# mean(fd[which(abs(fd$y) > 30),"FEve"]); sd(fd[which(abs(fd$y) > 30),"FEve"])
#
# summary(fd[,"FDiv"]); sd(fd[,"FDiv"])
# mean(fd[which(abs(fd$y) < 30),"FDiv"]); sd(fd[which(abs(fd$y) < 30),"FDiv"])
# mean(fd[which(abs(fd$y) > 30),"FDiv"]); sd(fd[which(abs(fd$y) > 30),"FDiv"])
#
# cor(fd$SR, fd$FDis, method = "spearman")
#
# mean(fd$Jac)
# sd(fd$Jac)
# mean(fd[which(abs(fd$y) < 30),"Jac"]); sd(fd[which(abs(fd$y) < 30),"Jac"])
# mean(fd[which(abs(fd$y) > 60),"Jac"]); sd(fd[which(abs(fd$y) > 60),"Jac"])
#
# cor(fd$Jac, fd$Jtu, method = "spearman")
#
# mean(fd$Bratio)
# sd(fd$Bratio)

### Find inflection points for the 2nd degree polynomials
# fit <- lm(formula = Faith ~ poly(SR,2), data = fd)
# new <- data.frame(SR = seq(from = min(fd$SR), to = max(fd$SR), by = 1))
# out <- predict(fit, newdata = new)
# diff <- diff(out)
# (diff/diff[1])*100
# new[52,] # 83
# new[78,] # 110
#
# fit <- lm(formula = Jne ~ poly(SR,2), data = fd)
# new <- data.frame(SR = seq(from = min(fd$SR), to = max(fd$SR), by = 1))
# out <- predict(fit, newdata = new)
# diff <- diff(out)
# (diff/diff[1])*100
# new[43,] # 82


### --------------------------------------------------------------------------

### B) Fig.2: Panel of mean annual beta functional diversity indices

jac <- ggplot() + geom_raster(aes(x = x, y = y, fill = Jac), data = fd) +
    scale_fill_gradientn(name = "Trait dissimilarity\n(Jaccard index)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .05, size = .4, aes(x = x, y = y, z = Jac), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

jtu <- ggplot() + geom_raster(aes(x = x, y = y, fill = Jtu), data = fd) +
    scale_fill_gradientn(name = "Trait turnover\n(Jtu)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .05, size = .4, aes(x = x, y = y, z = Jtu), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

jne <- ggplot() + geom_raster(aes(x = x, y = y, fill = Jne), data = fd) +
    scale_fill_gradientn(name = "Trait nestedness\n(Jne)", colours = parula(100), guide = "colourbar") +
    geom_contour(colour = "black", binwidth = .02, size = .4, aes(x = x, y = y, z = Jne), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

ratio <- ggplot() + geom_raster(aes(x = x, y = y, fill = Bratio), data = fd) +
    scale_fill_gradientn(name = "Ratio (Jtu/Jac)", colours = parula(100), guide = "colourbar", limits = c(0.45,0.88)) +
    geom_contour(colour = "black", binwidth = .05, size = .4, aes(x = x, y = y, z = Bratio), data = fd) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

# Save 
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
fig2 <- ggarrange(jac,jtu,jne,ratio, align = 'hv', ncol = 2, nrow = 2, labels = letters[1:4])
#ggsave(plot = fig2, filename = "Fig.2_mean_ann_beta_fd_indices_23.11.23.jpg", dpi = 300, width = 12, height = 4.5)
#ggsave(plot = fig2, filename = "Fig.2_mean_ann_beta_fd_indices_23.11.23.pdf", dpi = 300, width = 12, height = 4.5)

ggsave(plot = fig2, filename = "Fig.2_mean_ann_beta_fd_indices_revised_18.11.24.jpg", dpi = 300, width = 12, height = 4.5)
ggsave(plot = fig2, filename = "Fig.2_mean_ann_beta_fd_indices_revised_18.11.24.pdf", dpi = 300, width = 12, height = 4.5)


ggsave(plot = jac, filename = "map_mean_ann_jac_23.11.23.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = jtu, filename = "map_mean_ann_jtu_23.11.23.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = jne, filename = "map_mean_ann_jne.faith_23.11.23.jpg", dpi = 300, height = 4, width = 7) 
ggsave(plot = ratio, filename = "map_mean_ann_beta.ratio_23.11.23.jpg", dpi = 300, height = 4, width = 7)


### --------------------------------------------------------------------------

### C) Fig.4: Global emergent links between FD indices and Mesozooplankton biomass

# install.packages("ggpmisc")
library("ggpmisc")
# ?stat_poly_eq
# Define the linear model formula - simple
formula1 <- y~x # summary(lm)

p1 <- ggplot(fd, aes(x = SR, y = MESOZOO, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "right", size = 3.5) + 
  xlab("Species richness") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw()

p2 <- ggplot(fd, aes(x = Faith, y = MESOZOO, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "right", size = 3.5) + 
  xlab("Faith index") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw()

p3 <- ggplot(fd, aes(x = FEve, y = MESOZOO, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "left", size = 3.5) + 
  xlab("Functional eveness (FEve)") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw()

p4 <- ggplot(fd, aes(x = FDis, y = MESOZOO, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "right", size = 3.5) + 
  xlab("Functional dispersion (FDis)") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw()

p5 <- ggplot(fd, aes(x = FDiv, y = MESOZOO, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "left", size = 3.5) + 
  xlab("Functional divergence (FDiv)") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw()

p6 <- ggplot(fd, aes(x = Jac, y = MESOZOO, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "left", size = 3.5) + 
  xlab("Trait dissimilarity (Jaccard)") + ylab("Mesozooplankton biomass (mgC.m-3)") + theme_bw()

# Save
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
fig4 <- ggarrange(p1,p2,p3,p4,p5,p6, align = 'hv', ncol = 2, nrow = 3, labels = letters[2:7], common.legend = T)
fig4v2 <- ggarrange(p1,p2,p3,p4,p5,p6, align = 'hv', ncol = 3, nrow = 2, labels = letters[2:7], common.legend = T)
ggsave(plot = fig4, filename = "Fig.4_B-EF_fd_biomass_08.03.24.jpg", dpi = 300, width = 9.5, height = 13)
ggsave(plot = fig4, filename = "Fig.4_B-EF_fd_biomass_08.03.24.pdf", dpi = 300, width = 9.5, height = 13)
ggsave(plot = fig4v2, filename = "Fig.4v2_B-EF_fd_biomass_08.03.24.jpg", dpi = 300, width = 13, height = 9)
ggsave(plot = fig4v2, filename = "Fig.4v2_B-EF_fd_biomass_08.03.24.pdf", dpi = 300, width = 13, height = 9)


### 15/03/24: Modify order of plots so X axis variables match vertical order of appearance on heatmap
fig4v3 <- ggarrange(p6,p3,p1,p5,p2,p4, align = 'hv', ncol = 3, nrow = 2, labels = letters[2:7], common.legend = T)
ggsave(plot = fig4v3, filename = "Fig.4v3_B-EF_fd_biomass_15.03.24.jpg", dpi = 300, width = 13, height = 9)
ggsave(plot = fig4v3, filename = "Fig.4v3_B-EF_fd_biomass_15.03.24.pdf", dpi = 300, width = 13, height = 9)


### Other variables of interest: "FPOCex","POCflux","eratio" 
extra.vars <- c("FPOCex","POCflux","eratio")
# v <- "FPOCex"
for(v in extra.vars) {
    
    formula1 <- y~x # summary(lm)

    if(v == "eratio") {
        fd2 <- fd[fd$eratio < 0.4,]
    } else {
        fd2 <- fd
    } # eo if else loop
    
    p1 <- ggplot(fd2, aes(x = SR, y = get(v), colour = abs(y))) +
      geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
      scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
      stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1) + 
      xlab("Species richness") + ylab(v) + theme_bw()

    p2 <- ggplot(fd2, aes(x = Faith, y = get(v), colour = abs(y))) +
      geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
      scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
      stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1) + 
      xlab("Faith index") + ylab(v) + theme_bw()

    p3 <- ggplot(fd2, aes(x = FEve, y = get(v), colour = abs(y))) +
      geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
      scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
      stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1) + 
      xlab("Functional eveness (FEve)") + ylab(v) + theme_bw()

    p4 <- ggplot(fd2, aes(x = FDis, y = get(v), colour = abs(y))) +
      geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
      scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
      stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1) + 
      xlab("Functional dispersion (FDis)") + ylab(v) + theme_bw()

    p5 <- ggplot(fd2, aes(x = FDiv, y = get(v), colour = abs(y))) +
      geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
      scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
      stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1) + 
      xlab("Functional divergence (FDiv)") + ylab(v) + theme_bw()

    p6 <- ggplot(fd2, aes(x = Jac, y = get(v), colour = abs(y))) +
      geom_point() + geom_smooth(colour = "black", method = "lm", formula = formula1) +
      scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
      stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1) + 
      xlab("Trait dissimilarity (Jaccard)") + ylab(v) + theme_bw()

    figX <- ggarrange(p1,p2,p3,p4,p5,p6, align = 'hv', ncol = 2, nrow = 3, labels = letters[1:6], common.legend = T)
    ggsave(plot = figX, filename = paste("Fig.XX_B-EF_fd_",v,"_23.11.23.jpg"), dpi = 300, width = 8, height = 11)
    
} # eo for loop


### --------------------------------------------------------------------------

### D) Fig.4: Clustergram (heatmap of z-scores + map of clusters)
rownames(fd) <- fd$cell_id

names2keep <- colnames(fd)[c(4:8,10:12,14,16,21:25,29,30,32:35,39,40)] ; names2keep
data4pca <- fd[,c("x","y",names2keep)]
# Adjust colnames
colnames(data4pca)[12] <- "RaoQ"
colnames(data4pca)[20] <- "NPP"
colnames(data4pca)[24] <- "Slope"

# Normalize SR for z-scores analyses
maxou <- max(data4pca$SR)
data4pca$SR <- data4pca$SR/maxou

# Log transform
data4pca[,c(13:22,25)] <- log10(data4pca[,c(13:22,25)])
data4pca <- na.omit(data4pca) 
#summary(data4pca); dim(data4pca)

# Center and scale values prioro to clustering
scaled.data4pca <- data4pca
colnames(scaled.data4pca)[3:25] <- c("Species Richn.","Faith","SES Faith","Trait dissim.","Trait turnover",
        "Beta ratio","FEve","FDis","FDiv","Rao's Q","CHL-A","DIATO","DINO","GREEN","HAPTO","PROCHLO",
        "PROKAR","NPP","FPOC","POC FLUX","E RATIO","PSD SLOPE","MESOZOO")
scaled.data4pca[,c(3:length(scaled.data4pca))] <- base::scale(scaled.data4pca[,c(3:length(scaled.data4pca))], center = T, scale = T)

# Clusters vars 
vars_clust_ward <- hclust(dist(t(scaled.data4pca[,c(3:length(scaled.data4pca))])), method = "ward.D2")
# Cluster space 
space_clust <- hclust(dist(scaled.data4pca[,c(3:length(scaled.data4pca))]), method = "ward.D2")
# plot(space_clust)
klusters <- cutree(space_clust, k = 6)
scaled.data4pca$k <- as.numeric(klusters) 

# Find appropariate combination of colour from parula
pal.clusters <- rev(parula(6)) # 1st colour is too bright --> adjust
pal.clusters[1] <- "#fbe30e"

fig.4b.mercator <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(k)), data = scaled.data4pca) +
    scale_fill_manual(name = "Regions", values = pal.clusters) +
    geom_contour(colour = "black", binwidth = 1, size = .4, aes(x = x, y = y, z = k), data = scaled.data4pca) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

# Rotate longitudes for Mollweide projection
scaled.data4pca$x2 <- scaled.data4pca$x
scaled.data4pca[scaled.data4pca$x < 0,"x2"] <- (scaled.data4pca[scaled.data4pca$x < 0,"x"])+360

fig.4b.mollweide <- ggplot() + geom_tile(aes(x = x2, y = y, fill = factor(k)), data = scaled.data4pca) +
    geom_contour(colour = "grey25", binwidth = 1, size = .4, aes(x = x2, y = y, z = k), data = scaled.data4pca) +
    scale_fill_manual(name = "Regions", values = pal.clusters) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
         panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
    theme(legend.text=element_text(size=6), legend.title = element_text(size=8)) +
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    coord_map("mollweide", orientation = c(90,-180,0))

setwd("/net/kryo/work/fabioben/GODLY/plots")
ggsave(plot = fig.4b.mercator, filename = "map_clusters_mean_ann_FD+BCP_ward_23.11.23.pdf", dpi = 300, width = 7, height = 4)
ggsave(plot = fig.4b.mollweide, filename = "map_clusters_mean_ann_FD+BCP+BIOM_ward_23.11.23_v3.pdf", dpi = 300, width = 7, height = 4)
ggsave(plot = fig.4b.mollweide, filename = "map_clusters_mean_ann_FD+BCP+BIOM_ward_23.11.23_v3.jpg", dpi = 300, width = 7, height = 4)

### Clustergram (Fig. 4a)
d_vars_clust_ward <- as.dendrogram(vars_clust_ward)
d_space_clust <- as.dendrogram(space_clust)

# Find best looking divergent color palette
# palette <- rev(brewer.pal(11,"RdBu"))
# palette
# palette <- c("#053061","#2166AC","white","#B2182B","#67001F")
# ramp <- colorRampPalette(palette)(100)
# pal.sineramp(ramp)

### Really nice palette, but not centered around 0 for now
# https://stackoverflow.com/questions/10985224/r-heatmap-with-diverging-colour-palette/10986203#10986203
# https://stackoverflow.com/questions/29262824/r-center-color-palette-on-0
nHalf <- 50
Min <- round(min(as.matrix(scaled.data4pca[,c(3:25)])),1)
Max <- round(max(as.matrix(scaled.data4pca[,c(3:25)])),1)
Thresh <- 0

# Make vector of colors for values below threshold
rc1 <- colorRampPalette(colors = c("#053061","#2166AC","white"), space = "Lab")(nHalf)    
# Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = c("white","#B2182B","#67001F"), space = "Lab")(nHalf)
rampcols <- c(rc1,rc2)
## In your example, this line sets the color for values between -1 and 1
rampcols[c(nHalf,nHalf+1)] <- rgb(t(col2rgb("white")), maxColorValue = 256) 
rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1,rb2)
        
### Check if ramp is actually centred on 0?
# ggplot(data = scaled.data4pca, aes(x = y, y = PROCHLO, colour = PROCHLO)) + geom_point() +
#      scale_colour_gradientn(name = "PROCHLO", colours = rampcols, guide = "colourbar") +
#      geom_hline(yintercept = 0) + theme_bw()
### Good color bar, but not centred on 0° for every variable !!

jpeg("Fig.4a_heatmap_ward+ward_24.11.23_v2.jpeg", height = 6, width = 6, units = 'in', res = 600)
heatmap(x = as.matrix(scaled.data4pca[,c(3:25)]),
        Rowv = d_space_clust,
        Colv = d_vars_clust_ward,
        scale = "none",
        col = rampcols, 
        breaks = rampbreaks,
        cexCol = 0.8)
dev.off()

pdf("Fig.4a_heatmap_ward+ward_24.11.23_v2.pdf", height = 6, width = 6)
heatmap(x = as.matrix(scaled.data4pca[,c(3:25)]),
        Rowv = d_space_clust,
        Colv = d_vars_clust_ward,
        scale = "none",
        col = rampcols, 
        breaks = rampbreaks,
        cexCol = 0.8)
dev.off()


### --------------------------------------------------------------------------

### 14/03/24: Alternative versions of the heatmap to highlight latitudinal patterns
### (and forget about the regions, not so new and don't add anything substantial to the MS)

### Version A) Add absolute latitude in the variables of the z-scores matrix (tried with normal latitude: patterns not clear enough)
names2keep <- colnames(fd)[c(4:8,10:12,14,16,21:25,29,30,32:35,39,40)] ; names2keep
data4pca <- fd[,c("x","y",names2keep)]
colnames(data4pca)[2] <- "Latitude"
colnames(data4pca)[12] <- "RaoQ"
colnames(data4pca)[20] <- "NPP"
colnames(data4pca)[24] <- "Slope"
maxou <- max(data4pca$SR)
data4pca$SR <- data4pca$SR/maxou
data4pca[,c(13:22,25)] <- log10(data4pca[,c(13:22,25)])
data4pca$Latitude <- abs(data4pca$Latitude)
data4pca <- na.omit(data4pca) 
# Center and scale values prioro to clustering
scaled.data4pca <- data4pca
colnames(scaled.data4pca)[3:25] <- c("Species Richn.","Faith","SES Faith","Trait dissim.","Trait turnover",
        "Beta ratio","FEve","FDis","FDiv","Rao's Q","CHL-A","DIATO","DINO","GREEN","HAPTO","PROCHLO",
        "PROKAR","NPP","FPOC","POC FLUX","E RATIO","PSD SLOPE","MESOZOO")
scaled.data4pca[,c(2:length(scaled.data4pca))] <- base::scale(scaled.data4pca[,c(2:length(scaled.data4pca))], center = T, scale = T)

# Clusters vars & cells, without latitude's z-scores
vars_clust_ward <- hclust(dist(t(scaled.data4pca[,c(2:length(scaled.data4pca))])), method = "ward.D2")
space_clust <- hclust(dist(scaled.data4pca[,c(2:length(scaled.data4pca))]), method = "ward.D2")

d_space_clust <- as.dendrogram(space_clust)
d_vars_clust_ward <- as.dendrogram(vars_clust_ward)

setwd("/net/kryo/work/fabioben/GODLY/plots")
pdf("Fig.4a_heatmap_test_latitude_14.03.24.pdf", height = 6, width = 6)
heatmap(x = as.matrix(scaled.data4pca[,c(2:25)]),
        Rowv = d_space_clust,
        Colv = d_vars_clust_ward,
        scale = "none",
        col = rampcols, 
        breaks = rampbreaks,
        cexCol = 0.8)
dev.off()

jpeg("Fig.4a_heatmap_test_latitude_14.03.24.jpg", height = 6, width = 6, units = 'in', res = 600)
heatmap(x = as.matrix(scaled.data4pca[,c(2:25)]),
        Rowv = d_space_clust,
        Colv = d_vars_clust_ward,
        scale = "none",
        col = rampcols, 
        breaks = rampbreaks,
        cexCol = 0.8)
dev.off()


### Version B) Define 5° to 10° bins of latitude (absolute or not) and add a ColSideColors or RowSideColors
### representing latitude on the heatmap
names2keep <- colnames(fd)[c(4:8,10:12,14,16,21:25,29,30,32:35,39,40)] ; names2keep
data4pca <- fd[,c("x","y",names2keep)]
colnames(data4pca)[12] <- "RaoQ"
colnames(data4pca)[20] <- "NPP"
colnames(data4pca)[24] <- "Slope"
maxou <- max(data4pca$SR)
data4pca$SR <- data4pca$SR/maxou
data4pca[,c(13:22,25)] <- log10(data4pca[,c(13:22,25)])
data4pca <- na.omit(data4pca) 

# Center and scale values prior to clustering
scaled.data4pca <- data4pca
colnames(scaled.data4pca)[3:25] <- c("Species Richn.","Faith","SES Faith","Trait dissim.","Trait turnover",
        "Beta ratio","FEve","FDis","FDiv","Rao's Q","CHL-A","DIATO","DINO","GREEN","HAPTO","PROCHLO",
        "PROKAR","NPP","FPOC","POC FLUX","E RATIO","PSD SLOPE","MESOZOO")
scaled.data4pca[,c(3:25)] <- base::scale(scaled.data4pca[,c(3:25)], center = T, scale = T)

# Clusters vars & cells, without latitude's z-scores
vars_clust_ward <- hclust(dist(t(scaled.data4pca[,c(3:25)])), method = "ward.D2")
space_clust <- hclust(dist(scaled.data4pca[,c(3:25)]), method = "ward.D2")
d_space_clust <- as.dendrogram(space_clust)
d_vars_clust_ward <- as.dendrogram(vars_clust_ward)

# # Add bins of absolute latitude
# scaled.data4pca$Y <- abs(scaled.data4pca$y)
# # Cut into bins of 5° or 10°
# mround <- function(x,base){ base*round(x/base ) }
# scaled.data4pca$bins <- mround(scaled.data4pca$Y, 5)
# # Replace 0 by 1
# scaled.data4pca[scaled.data4pca$bins == 0,"bins"] <- 1
# # Define colorpalette for bins and associate bin values to a color value
# # length(scaled.data4pca$bins); length(unique(scaled.data4pca$bins))
# selcol2 <- colorRampPalette(parula(1000))
# colSide <- selcol2(length(unique(scaled.data4pca$bins)))
# lat.clust <- hclust(dist(scaled.data4pca[,"Y"]), method = "average")
# # plot(lat.clust)
# k = 9
# lat.clusters <- cutree(lat.clust, k = k)
# selcol2 <- colorRampPalette(rev(parula(k)))
# clustcol.height <- selcol2(k)
#
# jpeg("Fig.4a_heatmap_test_latitude_RSC_14.03.24.jpg", height = 6, width = 6, units = 'in', res = 600)
# heatmap(x = as.matrix(scaled.data4pca[,c(3:25)]),
#         Rowv = d_space_clust,
#         Colv = d_vars_clust_ward,
#         scale = "none",
#         col = rampcols,
#         breaks = rampbreaks,
#         cexCol = 0.8,
#         RowSideColors = clustcol.height[lat.clusters])
# dev.off()


### Better version: cut into 8 bins of absolute latitude
k <- 9
scaled.data4pca$bins <- cut_number(abs(scaled.data4pca$y), n = k)
data4leftjoin <- data.frame(bins = unique(scaled.data4pca$bins), rank = rank( unique(scaled.data4pca$bins) ))
scaled.data4pca2 <- dplyr::left_join(scaled.data4pca, data4leftjoin, by = "bins")
selcol2 <- colorRampPalette(rev(parula(k)))
clustcol.height <- selcol2(k)

jpeg("Fig.4a_heatmap_test_latitude_RSC_14.03.24.jpg", height = 6, width = 6, units = 'in', res = 600)
heatmap(x = as.matrix(scaled.data4pca[,c(3:25)]),
        Rowv = d_space_clust,
        Colv = d_vars_clust_ward,
        scale = "none",
        col = rampcols, 
        breaks = rampbreaks,
        cexCol = 0.8,
        RowSideColors = clustcol.height[scaled.data4pca2$rank])
dev.off()

### Plot palette
data4leftjoin$color <- clustcol.height
pal <- pal.bands(rev(data4leftjoin$color), n = k)

### --------------------------------------------------------------------------

### D) Fig.5: Panel of ∆FD indices (future - baseline): ∆SR, ∆
indices <- c("SR","FR","FEve","FDis","RaoQ","FDiv","Jac","Jtu","Jne")
# i <- "RaoQ" # for testing

maps <- lapply(indices, function(i) {
    
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
        files <- dir()[grep(i,dir())] # f <- files[1]
    
        res <- mclapply(files, function(f) {
                d <- get(load(f))
                filename <- str_replace_all(f,".Rdata","")
                d$ESM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
                d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
                d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[7]
                return(d)
            }, mc.cores = 20
        ) # eo mclapply
        tab <- bind_rows(res)
        rownames(tab) <- NULL
        rm(res); gc()
        # dim(tab); head(tab); summary(tab)
        
        # Compute ensemble mean for 'perc'
        ens <- data.frame(tab %>% group_by(cell_id) %>% 
                summarize(x = unique(x), y = unique(y),
                    mean = mean(perc, na.rm = T),
                    std = sd(perc, na.rm = T)
                )
        ) # eo ddf
        # dim(ens); summary(ens)
        # ggplot(ens, aes(x=mean)) + geom_histogram(binwidth = 10, colour="black", fill="white")
    
        # Rotate to 0-360° if needed (Niki usually prefers it)
        # ens$x2 <- ens$x
        # ens[ens$x < 0,"x2"] <- (ens[ens$x > 179.5,"x"])+360
                    
        # SR: limit colorbar between -25 and +50; adjust stat_contour 
        # FR: change name to Faith; add one contour level
        # FEve: Bound between -31 and +30, less contour maybe
        # FDis: one more contour
        # FDiv: one more contour
        # Jac & Jtu: one more contour
        # Jne: limit colorbar between -25 and +50; adjust stat_contour 
                    
        if(i == "SR") {
           
            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean*100), data = ens[ens$mean <= 0.3,]) +
                 geom_raster(aes(x = x, y = y), data = ens[ens$mean > 0.3,], fill = "#b10026") +
                 scale_fill_gradient2(name = paste("∆Richness", sep = ""), low = "#3D50C3", mid = "white", high = "#f50735", limits = c(-30,30)) +
                 geom_contour(colour = "grey30", binwidth = 15, size = 0.25, aes(x = x, y = y, z = mean*100), data = ens) +
                 geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
                 coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                 panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
                 scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
                 scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
 
        } else if(i == "FR") {
            
            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean*100), data = ens) +
                 scale_fill_gradient2(name = paste("∆Faith", sep = ""), low = "#3D50C3", mid = "white", high = "#f50735") +
                 geom_contour(colour = "grey30", binwidth = 5, size = 0.25, aes(x = x, y = y, z = mean*100), data = ens) +
                 geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
                 coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                 panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
                 scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
                 scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
            
        } else if(i == "FEve") {
            
            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean*100), data = ens[ens$mean <= 0.3,]) +
                 geom_raster(aes(x = x, y = y), data = ens[ens$mean > 0.3,], fill = "#b10026") +
                 scale_fill_gradient2(name = paste("∆FEve", sep = ""), low = "#3D50C3", mid = "white", high = "#f50735", limits = c(-31,30)) +
                 geom_contour(colour = "grey30", binwidth = 15, size = 0.25, aes(x = x, y = y, z = mean*100), data = ens) +
                 geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
                 coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                 panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
                 scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
                 scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
  
            
        } else if(i == "FDis") {
            
            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean*100), data = ens) +
                 scale_fill_gradient2(name = paste("∆FDis", sep = ""), low = "#3D50C3", mid = "white", high = "#f50735") +
                 geom_contour(colour = "grey30", binwidth = 1, size = 0.25, aes(x = x, y = y, z = mean*100), data = ens) +
                 geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
                 coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                 panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
                 scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
                 scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
            
        } else if(i == "FDiv") {
            
            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean*100), data = ens) +
                scale_fill_gradient2(name = paste("∆FDiv", sep = ""), low = "#3D50C3", mid = "white", high = "#f50735") +
                geom_contour(colour = "grey30", binwidth = .5, size = 0.25, aes(x = x, y = y, z = mean*100), data = ens) +
                geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
                coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
                scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
                scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
            
        } else if(i == "Jac") {
            
            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean*100), data = ens) +
                 scale_fill_gradient2(name = paste("∆Dissimilarity", sep = ""), low = "#3D50C3", mid = "white", high = "#f50735") +
                 geom_contour(colour = "grey30", binwidth = 5, size = 0.25, aes(x = x, y = y, z = mean*100), data = ens) +
                 geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
                 coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                 panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
                 scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
                 scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
            
        } else if(i == "Jtu") {
            
            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean*100), data = ens) +
                 scale_fill_gradient2(name = paste("∆Turnover", sep = ""), low = "#3D50C3", mid = "white", high = "#f50735") +
                 geom_contour(colour = "grey30", binwidth = 5, size = 0.25, aes(x = x, y = y, z = mean*100), data = ens) +
                 geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
                 coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                 panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
                 scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
                 scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
    
        } else if(i == "Jne") {
            
            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean*100), data = ens[ens$mean <= 0.35,]) +
                 geom_raster(aes(x = x, y = y), data = ens[ens$mean > 0.35,], fill = "#b10026") +
                 scale_fill_gradient2(name = paste("∆Nestedness", sep = ""), low = "#3D50C3", mid = "white", high = "#f50735", limits = c(-35,35)) +
                 geom_contour(colour = "grey30", binwidth = 20, size = 0.25, aes(x = x, y = y, z = mean*100), data = ens) +
                 geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
                 coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                 panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
                 scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
                 scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
            
        } else if(i == "RaoQ") {

            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = mean*100), data = ens[abs(ens$mean) <= 0.5,]) +
                 scale_fill_gradient2(name = paste("∆Rao's Q", sep = ""), low = "#3D50C3", mid = "white", high = "#f50735", limits = c(-5,5)) +
                 geom_contour(colour = "grey30", binwidth = 2.5, size = 0.25, aes(x = x, y = y, z = mean*100), data = ens) +
                 geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
                 coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                 panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
                 scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
                 scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
            
        } # eo if else loops
                              
        # Store maps in a list
        return(map)
    
    } # eo FUN
    
) # eo lapply 
# Check maps
# ggarrange(maps[[1]],maps[[2]],maps[[3]],maps[[4]],
#           maps[[6]],maps[[7]],maps[[8]],maps[[9]],
#           align = 'hv', ncol = 2, nrow = 4, labels = letters[1:10]
# )

### Use FDis instead of Rao's Q in panel
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")
fig5 <- ggarrange(maps[[1]],maps[[2]],maps[[3]],maps[[4]],maps[[6]],maps[[7]],maps[[8]],maps[[9]], align = 'hv', ncol = 2, nrow = 4, labels = letters[1:10] )

ggsave(plot = fig5, filename = "Fig.5_mean_ann_∆fd_indices_24.11.23.jpg", dpi = 300, width = 10, height = 8)
ggsave(plot = fig5, filename = "Fig.5_mean_ann_∆fd_indices_24.11.23.pdf", dpi = 300, width = 10, height = 8)

# Also save individual maps 
ggsave(plot = maps[[1]], filename = "map_mean_ann_∆SR_24.11.23.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = maps[[2]], filename = "map_mean_ann_∆FR_24.11.23.jpg.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = maps[[3]], filename = "map_mean_ann_∆FEve_24.11.23.jpg.jpg", dpi = 300, height = 4, width = 7) # Need to plot the p-values too here --> SM
ggsave(plot = maps[[4]], filename = "map_mean_ann_∆FDis_24.11.23.jpg.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = maps[[5]], filename = "map_mean_ann_∆RaoQ_24.11.23.jpg.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = maps[[6]], filename = "map_mean_ann_∆FDiv_24.11.23.jpg.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = maps[[7]], filename = "map_mean_ann_∆Jac_24.11.23.jpg.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = maps[[8]], filename = "map_mean_ann_∆Jtu_24.11.23.jpg.jpg", dpi = 300, height = 4, width = 7)
ggsave(plot = maps[[9]], filename = "map_mean_ann_∆Jne_24.11.23.jpg.jpg", dpi = 300, height = 4, width = 7)


### 17/12/23: Gather all mean (+sd) values in % change in a ddf to describe results in paper
indices <- c("SR","FR","FEve","FDis","RaoQ","FDiv","Jac","Jtu","Jne")
# i <- "SR"
maps <- lapply(indices, function(i) {
    
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
        files <- dir()[grep(i,dir())] # f <- files[1]
    
        res <- mclapply(files, function(f) {
                d <- get(load(f))
                filename <- str_replace_all(f,".Rdata","")
                d$ESM <- unlist(strsplit(x = filename, split = "_", fixed = T))[5]
                d$SDM <- unlist(strsplit(x = filename, split = "_", fixed = T))[6]
                d$month <- unlist(strsplit(x = filename, split = "_", fixed = T))[7]
                return(d)
            }, mc.cores = 20
        ) # eo mclapply
        tab <- bind_rows(res)
        rownames(tab) <- NULL
        rm(res); gc()
        # dim(tab); head(tab)
        
        # Compute ensemble mean for 'perc'
        ens <- data.frame(tab %>% group_by(cell_id) %>% 
                summarize(x = unique(x), y = unique(y),
                    mean = mean(perc, na.rm = T),
                    std = sd(perc, na.rm = T)
                )
        ) # eo ddf
        # dim(ens); summary(ens)

        ens$var <- i

        # Return
        return(ens)
        
    } # eo FUN
    
) # eo lapply
# Rbind
table <- bind_rows(maps)
rm(maps); gc()
# summary(table)
# Convert to 0-100%
table$mean <- (table$mean)*100

# mean(table[table$var == "SR","mean"]) ; sd(table[table$var == "SR","mean"])
# mean(table[table$var == "FR","mean"]) ; sd(table[table$var == "FR","mean"])
# mean(table[table$var == "FEve","mean"]) ; sd(table[table$var == "FEve","mean"])
# mean(table[table$var == "FDis","mean"]) ; sd(table[table$var == "FDis","mean"])
# mean(table[table$var == "FDiv","mean"]) ; sd(table[table$var == "FDiv","mean"])
# mean(table[table$var == "Jac","mean"]) ; sd(table[table$var == "Jac","mean"])
# mean(table[table$var == "Jtu","mean"]) ; sd(table[table$var == "Jtu","mean"])
# mean(table[table$var == "Jne","mean"]) ; sd(table[table$var == "Jne","mean"])

### Add raster of regions
setwd("/net/kryo/work/fabioben/BIOCeans5D/data/clusters")
ras <- get(load("raster_clusters_ward_k6_28.11.23.RData"))
table$region <- extract(ras, table[,c('x','y')]) # summary(factor(table$region))
table2 <- table[!is.na(table$region),]
#summary(factor(table2$region))

# mean(table2[table2$var == "Jac" & abs(table2$y) < 30,"mean"]); sd(table2[table2$var == "Jac" & abs(table2$y) < 30,"mean"])
# mean(table2[table2$var == "Jac" & abs(table2$y) > 30,"mean"]); sd(table2[table2$var == "Jac" & abs(table2$y) > 30,"mean"])

# mean(table2[table2$var == "Jne" & table2$region == 1,"mean"]); sd(table2[table2$var == "Jne" & table2$region == 1,"mean"])
# mean(table2[table2$var == "Jne" & table2$region == 2,"mean"]); sd(table2[table2$var == "Jne" & table2$region == 2,"mean"])
# mean(table2[table2$var == "Jne" & table2$region == 3,"mean"]); sd(table2[table2$var == "Jne" & table2$region == 3,"mean"])
# mean(table2[table2$var == "Jne" & table2$region == 4,"mean"]); sd(table2[table2$var == "Jne" & table2$region == 4,"mean"])
# mean(table2[table2$var == "Jne" & table2$region == 5,"mean"]); sd(table2[table2$var == "Jne" & table2$region == 5,"mean"])
# mean(table2[table2$var == "Jne" & table2$region == 6,"mean"]); sd(table2[table2$var == "Jne" & table2$region == 6,"mean"])

# cor(table2[table2$var == "Jac","mean"], table2[table2$var == "Jtu","mean"])

### --------------------------------------------------------------------------

### E°) Bivariate plots between FD indices and SR
library("ggpmisc")
formula1 <- y ~ x
formula2 <- y ~ poly(x,2)
formula3 <- y ~ poly(x,3)

p1 <- ggplot(fd, aes(x = SR, y = Faith, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula2) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula2, label.y = "bottom", label.x = "right", size = 3) + 
  xlab("Species richness") + ylab("Faith index") + theme_bw()

p2 <- ggplot(fd, aes(x = SR, y = SES.Faith, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "bottom", label.x = "left", size = 3) + 
  xlab("Species richness") + ylab("SES Faith index") + theme_bw()

p3 <- ggplot(fd, aes(x = SR, y = FEve, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "bottom", label.x = "left", size = 3) + 
  xlab("Species richness") + ylab("FEve") + theme_bw()

p4 <- ggplot(fd, aes(x = SR, y = FDis, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula2) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula2, label.y = "bottom", label.x = "right", size = 3) + 
  xlab("Species richness") + ylab("FDis") + theme_bw()
# Alternative: Rao's Q
p4.2 <- ggplot(fd, aes(x = SR, y = RaoQ.scaled, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula2) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula2, label.y = "bottom", label.x = "right", size = 3) + 
  xlab("Species richness") + ylab("Rao's Q") + theme_bw()

p5 <- ggplot(fd, aes(x = SR, y = FDiv, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula2) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula2, label.y = "bottom", label.x = "right", size = 3) + 
  xlab("Species richness") + ylab("FDiv") + theme_bw()

p6 <- ggplot(fd, aes(x = SR, y = Jac, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula2) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula2, label.y = "top", label.x = "right", size = 3) + 
  xlab("Species richness") + ylab("Trait dissimilarity (Jaccard)") + theme_bw()

p7.1 <- ggplot(fd, aes(x = SR, y = Jtu, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula2) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula2, label.y = "top", label.x = "right", size = 3) + 
  xlab("Species richness") + ylab("Trait turnover") + theme_bw()

p7.2 <- ggplot(fd, aes(x = SR, y = Jne, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula2) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula2, label.y = "top", label.x = "right", size = 3) + 
  xlab("Species richness") + ylab("Trait nestedness") + theme_bw()


# And Faith vs. Jac/Jtu/Jne
p8 <- ggplot(fd, aes(x = Faith, y = Jac, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula1) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula1, label.y = "top", label.x = "right", size = 3) + 
  xlab("Faith index") + ylab("Trait dissimilarity (Jaccard)") + theme_bw()

p9 <- ggplot(fd, aes(x = Faith, y = Jtu, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula2) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula2, label.y = "bottom", label.x = "left", size = 3) + 
  xlab("Faith index") + ylab("Trait turnover") + theme_bw()

p10 <- ggplot(fd, aes(x = Faith, y = Jne, colour = abs(y))) +
  geom_point(alpha = .1) + geom_smooth(colour = "black", method = "lm", formula = formula3) +
  scale_colour_gradientn(name = "Latitude", colours = rev(parula(100)), guide = "colourbar") +
  stat_poly_eq(use_label(c("eq","adj.R2","P")), formula = formula3, label.y = "top", label.x = "right", size = 3) + 
  xlab("Faith index") + ylab("Trait nedtedness") + theme_bw()

# Organize in panel and save; split last three plots in a different panel
setwd("/net/kryo/work/fabioben/BIOCeans5D/plots")

figX.1 <- ggarrange(p1,p2,p3,p4,p5,p6,p7.1,p7.2, align = 'hv', ncol = 2, nrow = 4, labels = letters, common.legend = T)
figX.2 <- ggarrange(p8,p9,p10, align = 'hv', ncol = 1, nrow = 3, labels = letters, common.legend = T)
figX.3 <- ggarrange(p9,p10, align = 'hv', ncol = 2, nrow = 1, labels = letters, common.legend = T)

ggsave(plot = figX.1, filename = "Fig.3.1_SRxFD_14.03.24.jpg", dpi = 300, width = 9.5, height = 15.5)
ggsave(plot = figX.1, filename = "Fig.3.1_SRxFD_14.03.24.pdf", dpi = 300, width = 8, height = 11)

ggsave(plot = figX.2, filename = "Fig.3.2_FDxbeta_14.03.24.jpg", dpi = 300, width = 5, height = 12.5)
ggsave(plot = figX.2, filename = "Fig.3.2_FDxbeta_14.03.24.pdf", dpi = 300, width = 5, height = 12.5)

ggsave(plot = figX.3, filename = "Fig.3.3_FDxbeta_19.11.24.jpg", dpi = 300, width = 9, height = 4.75)

### Very good. Should be integrated as Figs.3/4.
### Now fix the position of the stats on panel Fig; 3 (BEF)

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
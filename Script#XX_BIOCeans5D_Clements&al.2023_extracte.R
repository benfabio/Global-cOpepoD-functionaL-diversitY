
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 07/08/23: R script to extarct the Particles/Carbon flux data from Clements et al. (2023) (see: https://darchive.mblwhoilibrary.org/server/api/core/bitstreams/3c0e1f65-337b-5d11-adf5-df111f085b11/content) © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - extract the global reconstructions of particle biovolume, size distribution, and carbon export flux from the seasonal euphotic zone and maximum winter time mixed layer stored in the netCDF of Clements et al. (2023) 
# - make them available in a common raster or data.frame
# - use them to examine the emergent relationship with indices of copepod FD @ the global scale 

### Latest update: 07/08/23

### ------------------------------------------------------------------------------------------------------------------------------------------------------

library("ncdf4")
library("raster")
library("marmap")
library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("viridis")
library("ggpubr")

world <- map_data("world") 
world2 <- map_data("world2") 

setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/Clements&al.2023")

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### Open .nc file, examine variables 
nc <- nc_open("Clements&al._2023_Euphotic_Export.nc")
nc
# names(nc$var)
# Need to extract: 
variables <- c("Pred_BV","Pred_slope","Flux")
# Pred_BV = particle biovolume (ppm)
# Pred_slope = slope of particle size distribution (unitless)
# Flux =  carbon export flux from the seasonal euphotic zone (mgC/m2/d)

v <- "Flux" # for testing fun below

res <- lapply(variables, function(v) {
    
        # Message
        message(paste("Extracting ",v,"\n", sep = ""))
        
        lon <- ncvar_get(nc,"lon")
        lat <- ncvar_get(nc,"lat")
        var <- ncvar_get(nc,v)
        m.lon <- melt(lon)
        m.lat <- melt(lat)
        m.var <- melt(var)
        colnames(m.var) <- c("x","y","month","value")
        # unique(m.var$x); unique(m.var$y)
        # Adjust coordinates
        m.var$x2 <- m.var$x
        m.var[which(m.var$x > 180),"x2"] <- (m.var[which(m.var$x > 180),"x"]) - 360
        # unique(m.var$x2)
        m.var$y2 <- m.var$y - 90

        # Get closest 0.5 longitude and latitude
        longitudes <- seq(from = -179.5, to = 179.5, by = 1)
        latitudes <- seq(from = -89.5, to = 89.5, by = 1)
        
        # For each unique value of m.var$lon, find nearest value in 'longitudes' and attribute
        m.var$LONG <- NA
        # i <- 160.5000153
        for(i in unique(m.var$x2)) {
            # Find closest value to i in 'longitudes'
            index <- which.min(abs(longitudes - i))
            m.var[m.var$x2 == i,"LONG"] <- longitudes[index]
        } # eo for loop

        # For each unique value of m.var$lat, find nearest value in 'latitudes' and attribute
        m.var$LAT <- NA
        for(i in unique(m.var$y2)) {
            # Find closest value to i in 'longitudes'
            index <- which.min(abs(latitudes - i))
            m.var[m.var$y2 == i,"LAT"] <- latitudes[index]
        } # eo for loop
           
        # Add cell ID and compute average 
        m.var$ID <- paste(m.var$LONG, m.var$LAT, sep = "_")

        # Compute cliamtology within cell ID 
        clim <- data.frame( m.var %>% group_by(ID,month) %>% summarise(x = unique(LONG), y = unique(LAT), mean = mean(value, na.rm = T)) )
        # dim(clim) # good 
        
        # specify which variable
        clim$variable <- NA
        if(v == "Pred_BV") {
            clim$variable <- "Biovolume"
        } else if(v == "Pred_slope") {
            clim$variable <- "Slope"
        } else if(v == "Flux") {
            clim$variable <- "Export"
        } # eo else if loop
        
        # The .nc provided by the authors longitudes are actually shifted by 25°. Need to 'cut' the 'extra 25°' and put them back before 154.50
        clim$x <- (clim$x)+25
        clim$x2 <- clim$x
        clim[which(clim$x > 179.5),"x2"] <- (clim[which(clim$x > 179.5),"x"]) - 360
        # summary(clim)
        # ggplot() + geom_raster(aes(x = x2, y = y, fill = log10(mean)), data = na.omit(clim[clim$month == 9,])) +
#              geom_contour(colour = "grey70", binwidth = .1, size = 0.25, aes(x = x2, y = y, z = log10(mean)), data = na.omit(clim[clim$month == 9,])) +
#              scale_fill_viridis(name = "Particle biovolume\n(ppm)") +
#              geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "grey70", colour = "black", size = 0.3) +
#              coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#              panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#              scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#              scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
            
        # OK. Discard 'x2' 
        clim$x <- clim$x2
        clim <- subset(clim, select = -c(x2))
        # Correct ID
        clim$ID <- paste(clim$x, clim$y, sep = "_") # length(unique(clim$ID))
            
        # Return
        return(clim)
        
        message(paste("Returning ",v,"\n", sep = ""))
    
    } # eo fun
    
) # eo lapply vars

# Rbind
tab <- bind_rows(res)
# dim(tab); unique(tab$variable)
# summary(tab)
rm(res); gc()

### Draw maps
ggplot() + geom_raster(aes(x = x, y = y, fill = log10(mean)), data = na.omit(tab[tab$variable == "Biovolume",])) + 
    scale_fill_viridis(name = "Particle biovolume\n(ppm)") +
    geom_contour(colour = "grey30", binwidth = .1, size = 0.25, aes(x = x, y = y, z = log10(mean)), data = na.omit(tab[tab$variable == "Biovolume",])) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey80",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(month))

ggplot() + geom_raster(aes(x = x, y = y, fill = mean), data = na.omit(tab[tab$variable == "Slope",])) + 
    scale_fill_viridis(name = "PSD Slope") +
    geom_contour(colour = "grey30", binwidth = .2, size = 0.25, aes(x = x, y = y, z = mean), data = na.omit(tab[tab$variable == "Slope",])) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey80",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(month))

ggplot() + geom_raster(aes(x = x, y = y, fill = log10(mean)), data = na.omit(tab[tab$variable == "Export",])) + 
    scale_fill_viridis(name = "Carbon export\n(mgC/m2/d)") +
    geom_contour(colour = "grey30", binwidth = .2, size = 0.25, aes(x = x, y = y, z = log10(mean)), data = na.omit(tab[tab$variable == "Export",])) +
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey80",linetype = "dashed"), legend.position = "right") +
    scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL) +
    facet_wrap(.~ factor(month))
    
# Save 
save(x = tab, file = "table_clim_mon_part_BV_slope_carbon_Clements&al.2023_07_08_23.RData")    
    
nc_close(nc)

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------

##### 06/09/23: R script to extract the gridded satellite products from globcolour (CHLA and PFTs, 100km) and the standard VGPM NPP product and make monthly and annual climatologies from them on a 1°x1° cell grid © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - compute monthly and annual climatologies of CHLA concentration (+satellite PFTs and phyto size classes) --> 'globcolour_monthly_100km_CHL_REP'
# - compute monthly and annual climatologies of NPP derived from the standard VGPM algorithm (http://orca.science.oregonstate.edu/1080.by.2160.monthly.hdf.vgpm.m.chl.m.sst.php)

### Latest update: 06/09/23

library("tidyverse")
library("reshape2")
library("viridis")
library("parallel")
library("ncdf4")
library("raster")
library("marmap")

world <- map_data("world") # coastlines for maps
world2 <- map_data("world2") # coastlines for maps - Pacififc-centered

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### First, satellite CHLA (1997-2022) + contribution of phytoplankton PFTs

### To be consistent with NPP: keep 2003-2022 (complete years)
setwd("/net/kryo/work/datasets/gridded/ocean/2d/observation/chl/cmems_globcolour/globcolour_monthly_100km_CHL_REP/")#; dir()
files <- dir()[c(7:26)] # stick to complete years

# Examine one .nc file
nc <- nc_open(files[1])
vars <- names(nc$var)[c(1,4:12)]#; vars
nc_close(nc)

### For testing below 
# f <- files[1]
# v <- vars[1]

### For each file/year, extract all data and turn into monthyl 1°x1° data.frame
mclapply(files, function(f) {
                
            message(paste("\n",f,"\n", sep = ""))
            year <- substring(f,1,4) 
            
            for(v in vars) {
                
                setwd("/net/kryo/work/datasets/gridded/ocean/2d/observation/chl/cmems_globcolour/globcolour_monthly_100km_CHL_REP/")
                
                message(paste("Extracting ",v, sep = ""))
				
                # Get data from the nc file
				ras <- raster::stack(f, varname = v)
				# Turn into ddf and provide coordinates
				d <- as.data.frame(ras, xy = T)
                colnames(d)[c(3:length(d))] <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
                d$id <- paste(d$x,d$y, sep = "_") # length(unique(d$id))
                d <- d %>% relocate(id, .before = x)
                d <- d[order(d$id),]
                
                # Compute annual mean
                d$Annual <- rowMeans(x = as.matrix(d[,c(4:length(d))]), na.rm = TRUE)
                
                # Save ddf as data.frame in proper dir
                setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/CHLA_PFTs_globcolour_100km_2002-2022/")
                save(d, file = paste("table_clim_mon+ann_globcolour_100km_",year,"_",v,".Rdata", sep = ""))
                
                # Clean
                rm(d,ras); gc()
                    
            } # eo for loop - v in vars     
    
        }, mc.cores = 15
        
) # eo LAPPLY

setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/CHLA_PFTs_globcolour_100km_2002-2022/")
# For each vars, derive load all .Rdata and compute annual mean over the 2003-2022 period. Print in another .Rdata
mclapply(vars, function(v) {
    
        message(paste("\n","Extracting and computing annual climatology for ",v,"\n", sep = ""))
        files2load <- dir()[grep(paste("_",v, sep = ""),dir())]
        # f <- files2load[1]
        res <- lapply(files2load, function(f) { d <- get(load(f)) ; return(d[,c("id","x","y","Annual")]) } )
        # Rbind
        d <- bind_rows(res); rm(res); gc()
        # Compute annual clim
        clim <- data.frame( d %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), mean = mean(Annual, na.rm = T)) )
        
        # Save
        save(clim, file = paste("table_clim_ann_globcolour_100km_2003-2022_",v,".Rdata", sep = ""))
    
    }, mc.cores = length(vars)
    
) # eo mclapply - v in vars
# Good. Done.

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### For standard VGPM NPP product. NPP values are in units of mgC/m2/day

### NOTE: These NPP fields are not on the 100km regular grid:
longitudes <- seq(from = -179.5, to = 179.5, by = 1)
latitudes <- seq(from = -89.5, to = 89.5, by = 1)
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

setwd("/net/kryo/work/datasets/gridded/ocean/2d/observation/npp/xyz_files")
folders <- dir()[c(2:21)] # stick to complete years: 2003-2022
# fo <- folders[13]

### For each year/folder, extract the monthly NPP measurements, re-sample on 1x1 regular grid and compute average.
for(fo in folders) {
    
    setwd("/net/kryo/work/datasets/gridded/ocean/2d/observation/npp/xyz_files")
   
    year <- substring(fo,8,11) 
    message(paste("\n","Extracting and computing NPP climatologies for ",year,"\n", sep = ""))
    setwd(paste(getwd(),"/",fo, sep = ""))
    
    # vector of monthly files
    files <- dir() # index of files == index of month from Jan to Dec (e.g., 1st file = Jan, etc.)
    # f <- files[6]
    mclapply(files, function(f) {
        
                month <- which(files == f)
                message(paste("Month == ",month, sep = ""))
                # Read file
                d <- read.table(f, sep = "", head = T, dec = ".", skip = 1)
                colnames(d) <- c("x","y","NPP")
                # Replace -9999 by NA
                d <- na_if(d,-9999) # summary(d)
                
                # Attribute each cell to their closest long/lat coordinates on the classic 1x1 cell grid
                d$LONG <- NA
                # i <- 160.5000153
                for(i in unique(d$x)) {
                    # Find closest value to i in 'longitudes'
                    index <- which.min(abs(longitudes - i))
                    d[d$x == i,"LONG"] <- longitudes[index]
                } # eo for loop

                # Same with lat
                d$LAT <- NA
                for(i in unique(d$y)) {
                    index <- which.min(abs(latitudes - i))
                    d[d$y == i,"LAT"] <- latitudes[index]
                } # eo for loop
           
                # Add cell ID and compute average 
                d$id <- paste(d$LONG, d$LAT, sep = "_") # length(unique(d$id))

                # Compute cliamtology within cell ID 
                clim <- data.frame( d %>% group_by(id) %>% summarise(x = unique(LONG), y = unique(LAT), NPP = mean(NPP, na.rm = T)) )
                # Attribute month of interest
                clim$month <- months[month]
                # dim(clim); summary(clim)
                
                # Quick map to check
                # ggplot() + geom_raster(aes(x = x, y = y, fill = log10(NPP)), data = na.omit(clim)) +
#                     geom_contour(colour = "grey70", binwidth = .5, size = 0.25, aes(x = x, y = y, z = log10(NPP)), data = na.omit(clim)) +
#                     scale_fill_viridis(name = "NPP\n(logged)") +
#                     geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "grey70", colour = "black", size = 0.3) +
#                     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#                     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#                     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
        
                # Save 'clim' in proper dir
                setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/NPP_VGPM.standard_2002-2023")
                save(clim, file = paste("table_clim_mon_NPP_standard_VGPM_100km_",year,"_",months[month],".RData", sep = "") )
        
        }, mc.cores = length(files)
    
    ) # eo mclapply - f in files
    
} # eo for loop - F in FILES

### Good, now compute the 12 monthly climatologies as well as the annual climatology
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/NPP_VGPM.standard_2002-2023")

# For each m in months, load all clims and compute average 
# m <- "oct"
mclapply(months, function(m) {

         message(paste("\n","Computing monthly NPP climatology for ",m,"\n", sep = ""))
         files <- dir()[grepl(m,dir())]
         res <- lapply(files, function(f) { d <- get(load(f)) ; return(d)} ) # eo lapply
         # Rbind 
         d <- bind_rows(res) ; # dim(d); summary(d)
         rm(res); gc()
         
         # Compute monthly clim and save
         clim <- data.frame( d %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), NPP = mean(NPP, na.rm = T)) )
         
         save(clim, file = paste("table_clim_mon_NPP_standard_VGPM_2003-2022_",m,".RData", sep = "") )

     }, mc.cores = length(months) 
 
) # eo mclapply - m in months

### Now, load all files (those with '100km' in it), and compute annual mean
files <- dir()[grepl("100km",dir())]; files
res <- lapply(files, function(f) { d <- get(load(f)) ; return(d)} ) # eo lapply
d <- bind_rows(res)
rm(res); gc()

# Compute monthly clim and save
ann <- data.frame( d %>% group_by(id) %>% summarise(x = unique(x), y = unique(y), NPP = mean(NPP, na.rm = T)) )

# Quick map to check
# ggplot() + geom_raster(aes(x = x, y = y, fill = log10(NPP)), data = na.omit(ann)) +
#     geom_contour(colour = "grey70", binwidth = .33, size = 0.25, aes(x = x, y = y, z = log10(NPP)), data = na.omit(ann)) +
#     scale_fill_viridis(name = "Mean annual NPP\nlog10(mgC/m2/day)") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world[world$long <= 180,], fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

# Save
save(clim, file = paste("table_clim_ann_NPP_standard_VGPM_2003-2022.RData", sep = "") )

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
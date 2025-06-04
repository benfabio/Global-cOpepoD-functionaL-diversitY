
##### ATLANTECO SCRIPT 8.3 ----------------------------------------------------------------------------------------------------------------------------
##### 02/06/2025: R Script to create a standard netCDF file that will include all the global surface monthly estimates of copepod FD from Benedetti et al. (2025) for 'AtlantECO-MAPS-v2', based on Script#8.1 © Fabio Benedetti, ETH Zürich, IBP, UP Group.

# - Load a previous .nc file from 'AtlantECO-MAPS-v1' to create a first .nc file
# - Iteratively add the global (360x180) monthly (12) estimates of 8 copepod FD variables 
#   () + species richness from Benedetti et al. (2025) - Global Change Biology - https://onlinelibrary.wiley.com/doi/10.1111/gcb.70094 
# - Make one netCDF file for contemporary projections and a second one for the end-of-century (2100-2081) period

### Latest update: 04/06/2025 # Making the two netCDF files

library("tidyverse")
library("reshape2")
library("lubridate")
library("ncdf4")
library("raster")
library("terra")

### ----------------------------------------------------------------------------------------------------------------------------

### 02/05/2025: Go to the 'manuscript#2_Copepod_FD-EF' directory and get an overview of where ALL the monthly fields are
# setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)") # dir() 

### 1.1) Gather the data needed for the first .nc file (contemporary values)

## Species Richness & Faith index
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/future_projs/Faith/Faith")
files <- dir()[grepl("baseline",dir())] ; files
# f <- files[1]
res <- lapply(files, function(f) {
        # Useless message
        message(paste("Reading ",f, sep = ""))
        # Read f
        d <- get(load(f))
        # Add SDM and month from filename 
        d$SDM <- unlist(strsplit(f, "_"))[5]
        d$month <- NA
        if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "jan" ) {
            d$month <- 01
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "feb" ) {
            d$month <- 02
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "mar" ) {
            d$month <- 03
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "apr" ) {
            d$month <- 04
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "may" ) {
            d$month <- 05
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "jun" ) {
            d$month <- 06
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "jul" ) {
            d$month <- 07
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "aug" ) {
            d$month <- 08
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "sep" ) {
            d$month <- 09
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "oct" ) {
            d$month <- 10
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "nov" ) {
            d$month <- 11
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "dec" ) {
            d$month <- 12
        } # eo if else loop
        # Return
        return(d)
    } # eo FUN 
) # eo lapply - files
# Bind rows
table <- dplyr::bind_rows(res)

## First, add a cell id
table$ID <- paste(table$x, table$y, sep = "_")

## Compute Min/Mean/Max/Std per IDxmonth
# For SR
SR <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x),
            lat = unique(y),
            Min = min(SR, na.rm = TRUE),
            Mean = mean(SR, na.rm = TRUE),
            Median = median(SR, na.rm = TRUE),
            Max = max(SR, na.rm = TRUE),
            Stdev = sd(SR, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(SR)

# Save example ddf for chatGPT
# setwd("/Users/fabiobenedetti/Desktop/")
# save(x = t_SR, file = "ddf_for_ifa.Rdata")
# Go back to dir
# setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/future_projs/Faith/Faith")

# For Faith
Faith <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x),
            lat = unique(y),
            Min = min(Faith, na.rm = TRUE),
            Mean = mean(Faith, na.rm = TRUE),
            Median = median(Faith, na.rm = TRUE),
            Max = max(Faith, na.rm = TRUE),
            Stdev = sd(Faith, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(Faith)

rm(table,res); gc()


## SES
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/future_projs/Faith/SES")
files <- dir()[grepl("baseline",dir())] ; files
# f <- files[13]
res <- lapply(files, function(f) {
        # Useless message
        message(paste("Reading ",f, sep = ""))
        # Read f
        d <- get(load(f))
        # Add SDM and month from filename 
        d$SDM <- unlist(strsplit(f, "_"))[5]
        d$month <- NA
        if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "jan" ) {
            d$month <- 01
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "feb" ) {
            d$month <- 02
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "mar" ) {
            d$month <- 03
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "apr" ) {
            d$month <- 04
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "may" ) {
            d$month <- 05
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "jun" ) {
            d$month <- 06
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "jul" ) {
            d$month <- 07
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "aug" ) {
            d$month <- 08
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "sep" ) {
            d$month <- 09
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "oct" ) {
            d$month <- 10
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "nov" ) {
            d$month <- 11
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "dec" ) {
            d$month <- 12
        } # eo if else loop
        # Return
        return(d)
    } # eo FUN 
) # eo lapply - files
# Bind rows
table <- dplyr::bind_rows(res)
# dim(table); summary(table)
# Add cell id
table$ID <- paste(table$x, table$y, sep = "_")

SES <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x),
            lat = unique(y),
            Min = min(pd.obs.z, na.rm = TRUE),
            Mean = mean(pd.obs.z, na.rm = TRUE),
            Median = median(pd.obs.z, na.rm = TRUE),
            Max = max(pd.obs.z, na.rm = TRUE),
            Stdev = sd(pd.obs.z, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(SES)

rm(table,res); gc()


## FEve/FDis/FDiv
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/future_projs/db_FD/HSI_Gawdis_PCoA_Euclid")
files <- dir()[grepl("baseline",dir())] ; files
# f <- files[23]
res <- lapply(files, function(f) {
        # Useless message
        message(paste("Reading ",f, sep = ""))
        # Read f
        d <- get(load(f)) # str(d)
        # Add SDM and month from filename 
        d$SDM <- unlist(strsplit(f,"_"))[5]
        d$month <- NA
        if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "jan" ) {
            d$month <- 01
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "feb" ) {
            d$month <- 02
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "mar" ) {
            d$month <- 03
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "apr" ) {
            d$month <- 04
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "may" ) {
            d$month <- 05
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "jun" ) {
            d$month <- 06
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "jul" ) {
            d$month <- 07
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "aug" ) {
            d$month <- 08
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "sep" ) {
            d$month <- 09
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "oct" ) {
            d$month <- 10
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "nov" ) {
            d$month <- 11
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "dec" ) {
            d$month <- 12
        } # eo if else loop
        # Return
        return(d)
    } # eo FUN 
) # eo lapply - files
# Bind rows
table <- dplyr::bind_rows(res)
# dim(table); summary(table)
# Add cell id
table$ID <- paste(table$x, table$y, sep = "_")

# For FEve
FEve <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x),
            lat = unique(y),
            Min = min(FEve, na.rm = TRUE),
            Mean = mean(FEve, na.rm = TRUE),
            Median = median(FEve, na.rm = TRUE),
            Max = max(FEve, na.rm = TRUE),
            Stdev = sd(FEve, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(FEve)

# For FDis
FDis <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x),
            lat = unique(y),
            Min = min(FDis, na.rm = TRUE),
            Mean = mean(FDis, na.rm = TRUE),
            Median = median(FDis, na.rm = TRUE),
            Max = max(FDis, na.rm = TRUE),
            Stdev = sd(FDis, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(FDis)

# For FDiv
FDiv <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x),
            lat = unique(y),
            Min = min(FDiv, na.rm = TRUE),
            Mean = mean(FDiv, na.rm = TRUE),
            Median = median(FDiv, na.rm = TRUE),
            Max = max(FDiv, na.rm = TRUE),
            Stdev = sd(FDiv, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(FDiv)

rm(table,res); gc()


## Beta div indices (ßjac, ßjne & ßjtu)
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/future_projs/beta.div/")
files <- dir()[grepl("baseline",dir())] ; files
# f <- files[31]; f
res <- lapply(files, function(f) {
        # Useless message
        message(paste("Reading ",f, sep = ""))
        # Read f
        d <- get(load(f)) # str(d)
        # Add SDM and month from filename 
        d$SDM <- unlist(strsplit(f,"_"))[5]
        d$month <- NA
        if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "jan" ) {
            d$month <- 01
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "feb" ) {
            d$month <- 02
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "mar" ) {
            d$month <- 03
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "apr" ) {
            d$month <- 04
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "may" ) {
            d$month <- 05
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "jun" ) {
            d$month <- 06
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "jul" ) {
            d$month <- 07
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "aug" ) {
            d$month <- 08
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "sep" ) {
            d$month <- 09
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "oct" ) {
            d$month <- 10
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "nov" ) {
            d$month <- 11
        } else if( str_replace(unlist(strsplit(f,"_"))[6], ".Rdata", "") == "dec" ) {
            d$month <- 12
        } # eo if else loop
        # Return
        return(d)
    } # eo FUN 
) # eo lapply - files
# Bind rows
table <- dplyr::bind_rows(res)
# dim(table); summary(table)
# Add cell id
table$ID <- paste(table$x, table$y, sep = "_")

# For beta.jac
Trait_dissimilarity <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x),
            lat = unique(y),
            Min = min(beta.jac, na.rm = TRUE),
            Mean = mean(beta.jac, na.rm = TRUE),
            Median = median(beta.jac, na.rm = TRUE),
            Max = max(beta.jac, na.rm = TRUE),
            Stdev = sd(beta.jac, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(Trait_dissimilarity)

Trait_nestedness <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x),
            lat = unique(y),
            Min = min(beta.jne, na.rm = TRUE),
            Mean = mean(beta.jne, na.rm = TRUE),
            Median = median(beta.jne, na.rm = TRUE),
            Max = max(beta.jne, na.rm = TRUE),
            Stdev = sd(beta.jne, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(Trait_nestedness)

Trait_turnover <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x),
            lat = unique(y),
            Min = min(beta.jtu, na.rm = TRUE),
            Mean = mean(beta.jtu, na.rm = TRUE),
            Median = median(beta.jtu, na.rm = TRUE),
            Max = max(beta.jtu, na.rm = TRUE),
            Stdev = sd(beta.jtu, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(Trait_turnover)

rm(table,res); gc()


### -------------------------------------------------------

### 03/06/25: 1.2) Make the first NetCDF (contemporary projections)

### Examine the template netCDF file
setwd("/Users/fabiobenedetti/Desktop/")
# # Open the NetCDF file
# nc <- nc_open("AtlantECO-MAPS-v1_microbiome_monthly_surface_Calanoida_species_richness_2012-2031_Benedettietal.2021_20220919.nc")
# # Print summary of the file
# print(nc)
# # List all variable names
# names(nc$var)
# # List all dimension names and their lengths
# lapply(nc$dim, function(x) list(name = x$name, len = x$len))
# # Close the file
# nc_close(nc)


### Create new .nc file for SR first 

# === STEP 1: Load the template NetCDF file to get lon, lat, and month ===
template_nc <- nc_open("AtlantECO-MAPS-v1_microbiome_monthly_surface_Calanoida_species_richness_2012-2031_Benedettietal.2021_20220919.nc")
lon_vals <- ncvar_get(template_nc, "lon") # lon_vals
lat_vals <- ncvar_get(template_nc, "lat") # lat_vals 
month_vals <- ncvar_get(template_nc, "month") # month_vals
nc_close(template_nc)

# Names of the variables to put in the NetCDF
stat_names <- c("Mean","Median","Min","Max","Stdev")


# === STEP 2: Align data frames to full grid and  containing the FD values to 4D array [lon, lat, month, stat] ===

### Write FUN to convert df to 4D array 

# Check that our FD data.frames match the NetCDF in terms of dimensions
# length(unique(SR$lon)) # 360 - good
# length(unique(SR$lat)) # 161 ONLY - NOT GOOD
# summary(SR)  Min.   :-77.500 & Max.   : 82.500 

### -> Need to align your original data.frames to match the template NetCDF grid exactly BEFORE we do any array conversion

## Define the full grid (long lat month combinations)
full_grid <- expand.grid(
    lon = lon_vals,
    lat = lat_vals,
    month = month_vals
) 

## Let's fill in missing (lon, lat, month) combos with NAs for each FD indicator
# Merge all target ddf in a list
fd_names <- c("SR","Faith","SES","FDis","FDiv","FEve","Trait_dissimilarity","Trait_turnover","Trait_nestedness")
fd_df_list <- mget(fd_names) # class(fd_df_list)
names(fd_df_list) <- fd_names

## Align each data frames along full_grid
fd_df_list_aligned <- lapply(fd_df_list, function(df) {

    # Join full grid to fill in missing locations with NA
    full_df <- full_grid %>% left_join(df, by = c("lon","lat","month")) %>% arrange(lon,lat,month)

    # Ensure all stat columns are present and in correct order
    for(stat in stat_names) {
        if( !stat %in% names(full_df) ) {
            full_df[[stat]] <- NA  # Create missing stat column if needed
        } # eo if loop 
    } # eo for loop - stat
    
    return(full_df)
    
    } # eo FUN
     
) # eo lapply
# str(fd_df_list_aligned)

convert_df_to_array <- function(data, lon, lat, months) {
    
    array_out <- array(NA, dim = c(length(lon), length(lat), length(months), 5))

    i2 <- sapply(data$lon, function(x) which.min(abs(lon - x))) 
    j2 <- sapply(data$lat, function(y) which.min(abs(lat - y))) 
    l2 <- sapply(data$month, function(x) which(months == x)) 

    array_out[cbind(i2,j2,l2,1)] <- data$Mean
    array_out[cbind(i2,j2,l2,2)] <- data$Median
    array_out[cbind(i2,j2,l2,3)] <- data$Min
    array_out[cbind(i2,j2,l2,4)] <- data$Max
    array_out[cbind(i2,j2,l2,5)] <- data$Stdev
    
    return(array_out)
    
} # eo FUN - convert_df_to_array()

## Test convert_df_to_array()
#test_array <- convert_df_to_array(data = SR, lon = lon_vals, lat = lat_vals, months = months)
# str(test_array) ; summary(test_array)
# Seems to work now 

### Apply convert_df_to_array() to all ddf - takes a few minutes
fd_array_list <- lapply(fd_df_list_aligned, function(df) {
        convert_df_to_array(df, lon = lon_vals, lat = lat_vals, months = month_vals)
    } # eo FUN
) # eo lapply - df
## Checks
# str(fd_array_list)
# summary(fd_array_list[[2]]) ; summary(Faith) # Looks good


# === STEP 3: Define NetCDF dimensions (add FD_indicator) ===
lon_dim <- ncdim_def("lon", "degrees_east", vals = lon_vals)
lat_dim <- ncdim_def("lat", "degrees_north", vals = lat_vals)
month_dim <- ncdim_def("time", "months since 2012-01", vals = month_vals, unlim = TRUE)
fd_dim <- ncdim_def("FD indicator", "", vals = 1:9)  # 9 possible indicators


# === STEP 4: Define variables ===
make_var <- function(varname) {
    ncvar_def(
        name = varname,
        units = "",
        dim = list(lon_dim,lat_dim,month_dim,fd_dim),
        missval = -9999,
        prec = "float"
    ) # eo ncvar_def
} # eo FUN - make_var

var_list <- lapply(c("Mean","Median","Min","Max","Stdev"), make_var)


# === STEP 5: Create new NetCDF file ===
# new_nc <- nc_create("AtlantECO-MAPS-v2_microbiome_monthly_surface_copepod_FD_contemporary_Benedettietal.2025_20250604.nc", vars = var_list)


# # === STEP 6: Write data for the first FD indicator ===
# fd_index <- 1  # We're only adding one for now
#
# for(s in seq_along(stat_names)) {
#
#     message(paste("Filling in new .nc with ",stat_names[s]," values", sep = ""))
#
#     stat <- stat_names[s]
#     varname <- stat
#
#     ncvar_put(
#         new_nc,
#         varid = varname,
#         vals = fd_array[,,,s],
#         start = c(1, 1, 1, fd_index),
#         count = c(-1, -1, -1, 1)
#     )  # Write only for fd_index = 1
#
# } # eo for loop - s in 'stat_names'
# # Check
# new_nc
#
# # nc_close(new_nc)


### Create NetCDF File Including These Named Indicators

new_nc <- nc_create("AtlantECO-MAPS-v2_microbiome_monthly_surface_copepod_FD_contemporary_Benedettietal.2025_20250604.nc", vars = var_list, force_v4 = TRUE)

### Fill new_nc with the data from the arrays stored in a list
for(fd_index in c(1:9)) {
    
    message(paste("Filling ",fd_names[fd_index]," values in the new NetCDF", sep = ""))
    
    # Subset array of interest
    arr <- fd_array_list[[fd_index]]
    
    for(s in seq_along(stat_names)) {
        # s <- 3
        stat <- stat_names[s]
        
        ncvar_put(
            new_nc,
            varid = stat,
            vals = arr[,,,s],
            start = c(1,1,1,fd_index),
            count = c(-1,-1,-1,1)
        ) # eo - ncvar_put()
        
     } # eo 2nd for loop - s

} # eo 1st for loop - fd_index

print(new_nc)
names(new_nc$var)
lapply(new_nc$dim, function(x) list(name = x$name, len = x$len))


### Define Metadata for Each FD Indicator and putting them as attributes in new nc files

## Write metadata file
fd_metadata <- data.frame(
  
  name = fd_names,
  long_name = c(
    "Species Richness",
    "Faith's index for Functional Richness (Faith)",
    "Standardized Effect Size of Faith",
    "Functional Dispersion (FDis)",
    "Functional Divergence (FDiv)",
    "Functional Evenness (FEve)",
    "Trait Dissimilarity (Jaccard's dissimilarity index)",
    "Trait Turnover (turnover component of Jaccard's index)",
    "Trait Nestedness (nestedness component of Jaccard's index)"
  ),
  units = rep("unitless", 9),  # Adjust as needed
  description = c(
    "Potential number of species present per grid cell",
    "Assemblages with higher Faith values and more numerous the present species. Assemblages with are those where branches on the total volume filled by represent more distant functional dendrogram (i.e., more functional the assemblage)",
    "SES Faith values < 0 indicate that func-tional clustering (or functional convergence) occurs due to environmental filtering in the copepod assemblage whereas values > 0 indicate that there is functional overdispersion",
    "Assemblages with higher FDis values are those whose species are further away from each other and from the centroid in the functional space (i.e., more specialized species)",
    "Assemblages with higher FDiv values are characterized by higher HSI values at the vertices of their convex hull (i.e., more extreme traits values)",
    "Higher FEve values indicate that species in the assemblage display similar HSI at equal distances between nearest neighbors in the functional space. Lower FEve values indicate the coexistence of scattered clouds of functional units",
    "Trait dissimilarity values close to 1 indicate that two assemblages display functional dendrograms with very different number of branches that are non-overlapping",
    "Trait turnover values close to 1 indicate that total trait dissimilarity is driven by the replacement of branches",
    "Trait nestedness values close to 1 indicate that total trait dissimilarity is driven by different number of branches, whatever their identity"
  ), stringsAsFactors = FALSE
  
)

## Add FD indicator metadata to NetCDF
for(i in c(1:nrow(fd_metadata))) {
    
    # idx <- i
    
    ncatt_put(new_nc, varid = 0, 
            attname = paste0("FD_indicator_", i, "_name"), 
            attval = fd_metadata$name[i]
    )
  
    ncatt_put(new_nc, varid = 0, 
            attname = paste0("FD_indicator_", i, "_long_name"), 
            attval = fd_metadata$long_name[i]
    )
  
    ncatt_put(new_nc, varid = 0, 
            attname = paste0("FD_indicator_", i, "_description"), 
            attval = fd_metadata$description[i]
    )
  
    ncatt_put(new_nc, varid = 0, 
            attname = paste0("FD_indicator_", i, "_units"), 
            attval = fd_metadata$units[i]
    )
            
} # eo for loop

### Add global attributes to the new NetCDF
# if varid == 0, then a global attribute is written instead of a variable's attribute
ncatt_put(new_nc, varid = 0, attname = "Dataset", attval = "AtlantECO data - WP2 - Traditinal microscopy")
ncatt_put(new_nc, varid = 0,"Institution","ETH Zürich, D-USYS, IBP, UP group")
ncatt_put(new_nc, varid = 0,"Description","Mean/Median/Min/Max and Stdev of monthly projections of copepod species richness and eight functional diversity (FD) indicators for the contemporary global surface ocean stemming from three species distribution models (GLM, GAM and ANN). See full description in Benedetti et al. (2025): Emergent Relationships Between the Functional Diversity of Marine Planktonic Copepods and Ecosystem Functioning in the Global Ocean")
ncatt_put(new_nc, varid = 0,"DOI","doi.org/10.1111/gcb.70094")
ncatt_put(new_nc, varid = 0, "Funding statement","This project has received funding from the European Union's Horizon 2020 Research and Innovation Programme under grant agreement no. 862923 (AtlantECO) and under grant agreement no. 101059915 (BIOceans5D). This output reflects only the author's view, and the European Union cannot be held responsible for any use that may be made of the information contained therein.")
# Even add history
history <- paste("Fabio Benedetti (fabio.benedetti@unibe.ch); last update: ", date(), sep = ", ")
ncatt_put(new_nc, varid = 0, "Last update by:", history)


### Close nc - shpuld save the data 
nc_close(new_nc)

## Open the NetCDF file again to check
nc <- nc_open("AtlantECO-MAPS-v2_microbiome_monthly_surface_copepod_FD_contemporary_Benedettietal.2025_20250604.nc")
print(nc)
names(nc$var)
lapply(nc$dim, function(x) list(name = x$name, len = x$len))
nc_close(nc)


## Open as raster? 
# library("raster")
# Load a specific variable, e.g., "Mean"
# r <- raster::brick("AtlantECO-MAPS-v2_microbiome_monthly_surface_copepod_FD_contemporary_Benedettietal.2025_20250604.nc", varname = "Stdev")
# Plot quickly
# plot(r)
### Looks ok


### --------------------------------------------------------------------------------------------------------------

### 2.1) Gather the data needed for the second .nc file (end-of-century (2100-2081) period)

### BEWARE: FOR FUTURE PROJECTIONS, LONGITUDES ARE in 0°-360° -> Need to rotate to classic WGS84 

## SR & Faith
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/future_projs/Faith/Faith")
files <- dir()[grepl("2100-2000",dir())] ; files
# f <- files[13]
res <- lapply(files, function(f) {
        # Useless message
        message(paste("Reading ",f, sep = ""))
        # Read f
        d <- get(load(f))
        # Add ESM, SDM and month from filename 
        d$SDM <- unlist(strsplit(f, "_"))[6]
        d$ESM <- unlist(strsplit(f, "_"))[5]
        d$month <- NA
        if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "jan" ) {
            d$month <- 01
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "feb" ) {
            d$month <- 02
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "mar" ) {
            d$month <- 03
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "apr" ) {
            d$month <- 04
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "may" ) {
            d$month <- 05
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "jun" ) {
            d$month <- 06
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "jul" ) {
            d$month <- 07
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "aug" ) {
            d$month <- 08
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "sep" ) {
            d$month <- 09
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "oct" ) {
            d$month <- 10
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "nov" ) {
            d$month <- 11
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "dec" ) {
            d$month <- 12
        } # eo if else loop
        # Return
        return(d)
    } # eo FUN 
) # eo lapply - files
# Bind rows
table <- dplyr::bind_rows(res)
## Convert longitudes from [0, 360] to [-180, 180]
table$x2 <- ifelse(table$x > 180, table$x - 360, table$x)
## First, add a cell id
table$ID <- paste(table$x2, table$y, sep = "_")

## Compute Min/Mean/Max/Std per IDxmonth
# For SR
SR <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x2),
            lat = unique(y),
            Min = min(SR, na.rm = TRUE),
            Mean = mean(SR, na.rm = TRUE),
            Median = median(SR, na.rm = TRUE),
            Max = max(SR, na.rm = TRUE),
            Stdev = sd(SR, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(SR)

# Faith
Faith <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x2),
            lat = unique(y),
            Min = min(Faith, na.rm = TRUE),
            Mean = mean(Faith, na.rm = TRUE),
            Median = median(Faith, na.rm = TRUE),
            Max = max(Faith, na.rm = TRUE),
            Stdev = sd(Faith, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(Faith)

rm(table,res); gc()


## SES 
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/future_projs/Faith/SES")
files <- dir()[grepl("2100-2000",dir())] ; files
# f <- files[13]
res <- lapply(files, function(f) {
        # Useless message
        message(paste("Reading ",f, sep = ""))
        # Read f
        d <- get(load(f))
        # Add ESM, SDM and month from filename 
        d$SDM <- unlist(strsplit(f, "_"))[6]
        d$ESM <- unlist(strsplit(f, "_"))[5]
        d$month <- NA
        if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "jan" ) {
            d$month <- 01
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "feb" ) {
            d$month <- 02
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "mar" ) {
            d$month <- 03
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "apr" ) {
            d$month <- 04
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "may" ) {
            d$month <- 05
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "jun" ) {
            d$month <- 06
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "jul" ) {
            d$month <- 07
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "aug" ) {
            d$month <- 08
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "sep" ) {
            d$month <- 09
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "oct" ) {
            d$month <- 10
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "nov" ) {
            d$month <- 11
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "dec" ) {
            d$month <- 12
        } # eo if else loop
        # Return
        return(d)
    } # eo FUN 
) # eo lapply - files
# Bind rows
table <- dplyr::bind_rows(res)
## Convert longitudes from [0, 360] to [-180, 180]
table$x2 <- ifelse(table$x > 180, table$x - 360, table$x)
## First, add a cell id
table$ID <- paste(table$x2, table$y, sep = "_")

# For SES
SES <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x2),
            lat = unique(y),
            Min = min(pd.obs.z, na.rm = TRUE),
            Mean = mean(pd.obs.z, na.rm = TRUE),
            Median = median(pd.obs.z, na.rm = TRUE),
            Max = max(pd.obs.z, na.rm = TRUE),
            Stdev = sd(pd.obs.z, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
summary(SES)

rm(table,res); gc()


## FEve/FDis/FDiv
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/future_projs/db_FD/HSI_Gawdis_PCoA_Euclid")
files <- dir()[grepl("2100-2000",dir())] ; files
# f <- files[23]; f
res <- lapply(files, function(f) {
        # Useless message
        message(paste("Reading ",f, sep = ""))
        # Read f
        d <- get(load(f)) # str(d)
        # Add SDM and month from filename 
        d$SDM <- unlist(strsplit(f,"_"))[6]
        d$ESM <- unlist(strsplit(f,"_"))[5]
        d$month <- NA
        if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "jan" ) {
            d$month <- 01
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "feb" ) {
            d$month <- 02
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "mar" ) {
            d$month <- 03
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "apr" ) {
            d$month <- 04
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "may" ) {
            d$month <- 05
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "jun" ) {
            d$month <- 06
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "jul" ) {
            d$month <- 07
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "aug" ) {
            d$month <- 08
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "sep" ) {
            d$month <- 09
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "oct" ) {
            d$month <- 10
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "nov" ) {
            d$month <- 11
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "dec" ) {
            d$month <- 12
        } # eo if else loop
        # Return
        return(d)
    } # eo FUN 
) # eo lapply - files
# Bind rows
table <- dplyr::bind_rows(res)
## Convert longitudes from [0, 360] to [-180, 180]
table$x2 <- ifelse(table$x > 180, table$x - 360, table$x) # summary(table$x2)
# Add cell id
table$ID <- paste(table$x2, table$y, sep = "_")

# For FEve
FEve <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x2),
            lat = unique(y),
            Min = min(FEve, na.rm = TRUE),
            Mean = mean(FEve, na.rm = TRUE),
            Median = median(FEve, na.rm = TRUE),
            Max = max(FEve, na.rm = TRUE),
            Stdev = sd(FEve, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
summary(FEve)

# For FDis
FDis <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x2),
            lat = unique(y),
            Min = min(FDis, na.rm = TRUE),
            Mean = mean(FDis, na.rm = TRUE),
            Median = median(FDis, na.rm = TRUE),
            Max = max(FDis, na.rm = TRUE),
            Stdev = sd(FDis, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
summary(FDis)

# For FDiv
FDiv <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x2),
            lat = unique(y),
            Min = min(FDiv, na.rm = TRUE),
            Mean = mean(FDiv, na.rm = TRUE),
            Median = median(FDiv, na.rm = TRUE),
            Max = max(FDiv, na.rm = TRUE),
            Stdev = sd(FDiv, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
summary(FDiv)

rm(table,res); gc()


## Beta div indices (ßjac, ßjne & ßjtu)
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/ETHZ/OVERSEE/manuscript#2_Copepod_FD-EF/data (new)/future_projs/beta.div/")
files <- dir()[grepl("2100-2000",dir())] ; files
# f <- files[100]; f
res <- lapply(files, function(f) {
        # Useless message
        message(paste("Reading ",f, sep = ""))
        # Read f
        d <- get(load(f)) # str(d)
        # Add SDM and month from filename 
        d$SDM <- unlist(strsplit(f,"_"))[6]
        d$ESM <- unlist(strsplit(f,"_"))[5]
        d$month <- NA
        if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "jan" ) {
            d$month <- 01
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "feb" ) {
            d$month <- 02
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "mar" ) {
            d$month <- 03
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "apr" ) {
            d$month <- 04
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "may" ) {
            d$month <- 05
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "jun" ) {
            d$month <- 06
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "jul" ) {
            d$month <- 07
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "aug" ) {
            d$month <- 08
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "sep" ) {
            d$month <- 09
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "oct" ) {
            d$month <- 10
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "nov" ) {
            d$month <- 11
        } else if( str_replace(unlist(strsplit(f,"_"))[7], ".Rdata", "") == "dec" ) {
            d$month <- 12
        } # eo if else loop
        # Return
        return(d)
    } # eo FUN 
) # eo lapply - files
# Bind rows
table <- dplyr::bind_rows(res)
# dim(table); summary(table)
# Add cell id - THIS ONE ALREADY HAS CORRECT LONGITUDES
table$ID <- paste(table$x2, table$y, sep = "_")

# For overall trait dissimilarity (jaccard index)
Trait_dissimilarity <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x2),
            lat = unique(y),
            Min = min(beta.jac, na.rm = TRUE),
            Mean = mean(beta.jac, na.rm = TRUE),
            Median = median(beta.jac, na.rm = TRUE),
            Max = max(beta.jac, na.rm = TRUE),
            Stdev = sd(beta.jac, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(Trait_dissimilarity)

# For overall trait turnover
Trait_nestedness <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x2),
            lat = unique(y),
            Min = min(beta.jne, na.rm = TRUE),
            Mean = mean(beta.jne, na.rm = TRUE),
            Median = median(beta.jne, na.rm = TRUE),
            Max = max(beta.jne, na.rm = TRUE),
            Stdev = sd(beta.jne, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(Trait_nestedness)

# For overall trait nestedness
Trait_turnover <- data.frame(
    table %>% 
        group_by(ID,month) %>% 
        summarise(
            lon = unique(x2),
            lat = unique(y),
            Min = min(beta.jtu, na.rm = TRUE),
            Mean = mean(beta.jtu, na.rm = TRUE),
            Median = median(beta.jtu, na.rm = TRUE),
            Max = max(beta.jtu, na.rm = TRUE),
            Stdev = sd(beta.jtu, na.rm = TRUE)
        ) # eo summarise
) # eo ddf
# summary(Trait_turnover)

rm(table,res); gc()


### 2.2) Make the second NetCDF (future projections)

# === STEP 1: Load the template NetCDF file to get lon, lat, and month ===
setwd("/Users/fabiobenedetti/Desktop/")
template_nc <- nc_open("AtlantECO-MAPS-v1_microbiome_monthly_surface_Calanoida_species_richness_2012-2031_Benedettietal.2021_20220919.nc")
lon_vals <- ncvar_get(template_nc, "lon") # lon_vals
lat_vals <- ncvar_get(template_nc, "lat") # lat_vals 
month_vals <- ncvar_get(template_nc, "month") # month_vals
nc_close(template_nc)

# Names of the variables to put in the NetCDF
stat_names <- c("Mean","Median","Min","Max","Stdev")


# === STEP 2: Align data frames to full grid and  containing the FD values to 4D array [lon, lat, month, stat] ===

# Check that our FD data.frames match the NetCDF in terms of dimensions
# length(unique(SR$lon)) # 360 - good
# length(unique(SR$lat)) # 161 ONLY - NOT GOOD
# summary(SR)  Min.   :-77.500 & Max.   : 82.500 
### -> Again, need to align your original data.frames to match the template NetCDF grid exactly BEFORE we do any array conversion

## Define the full grid (long lat month combinations)
full_grid <- expand.grid(
    lon = lon_vals,
    lat = lat_vals,
    month = month_vals
) # eo expand.grid

## Let's fill in missing (lon, lat, month) combos with NAs for each FD indicator
fd_names <- c("SR","Faith","SES","FDis","FDiv","FEve","Trait_dissimilarity","Trait_turnover","Trait_nestedness")
fd_df_list <- mget(fd_names) # class(fd_df_list); str(fd_df_list)
names(fd_df_list) <- fd_names

## Align each data frames along full_grid
fd_df_list_aligned <- lapply(fd_df_list, function(df) {
    # Join full grid to fill in missing locations with NA
    full_df <- full_grid %>% left_join(df, by = c("lon","lat","month")) %>% arrange(lon,lat,month)
    # Ensure all stat columns are present and in correct order
    for(stat in stat_names) {
        if( !stat %in% names(full_df) ) {
            full_df[[stat]] <- NA
        } # eo if loop 
    } # eo for loop - stat
    return(full_df)
    } # eo FUN
) # eo lapply
# str(fd_df_list_aligned)


### Re-use the convert_df_to_array() FUN 
convert_df_to_array <- function(data, lon, lat, months) {
    
    array_out <- array(NA, dim = c(length(lon), length(lat), length(months), 5))

    i2 <- sapply(data$lon, function(x) which.min(abs(lon - x))) 
    j2 <- sapply(data$lat, function(y) which.min(abs(lat - y))) 
    l2 <- sapply(data$month, function(x) which(months == x)) 

    array_out[cbind(i2,j2,l2,1)] <- data$Mean
    array_out[cbind(i2,j2,l2,2)] <- data$Median
    array_out[cbind(i2,j2,l2,3)] <- data$Min
    array_out[cbind(i2,j2,l2,4)] <- data$Max
    array_out[cbind(i2,j2,l2,5)] <- data$Stdev
    
    return(array_out)
    
} # eo FUN - convert_df_to_array()

## Test convert_df_to_array()
# test_array <- convert_df_to_array(data = Faith, lon = lon_vals, lat = lat_vals, months = month_vals)
# str(test_array) ; summary(test_array)
# Seems to work now 


### Apply convert_df_to_array() to all ddf - takes a few minutes
fd_array_list <- lapply(fd_df_list_aligned, function(df) {
        convert_df_to_array(df, lon = lon_vals, lat = lat_vals, months = month_vals)
    } # eo FUN
) # eo lapply - df
## Checks
# str(fd_array_list) # Looks good
# Compare summary( fd_array_list[[1]] ) #against# summary(SR)
# Looks good


# === STEP 3: Define NetCDF dimensions (add FD_indicator) ===
lon_dim <- ncdim_def("lon", "degrees_east", vals = lon_vals)
lat_dim <- ncdim_def("lat", "degrees_north", vals = lat_vals)
month_dim <- ncdim_def("time", "months since 2012-01", vals = month_vals, unlim = TRUE)
fd_dim <- ncdim_def("FD indicator", "", vals = 1:9)


# === STEP 4: Define variables ===
make_var <- function(varname) {
    ncvar_def(
        name = varname,
        units = "",
        dim = list(lon_dim,lat_dim,month_dim,fd_dim),
        missval = -9999,
        prec = "float"
    ) # eo ncvar_def
} # eo FUN - make_var

var_list <- lapply(c("Mean","Median","Min","Max","Stdev"), make_var)


# === STEP 5: Create new NetCDF file ===

## Create it on the Desktop
setwd("/Users/fabiobenedetti/Desktop/")

# Create new .nc file
new_nc <- nc_create("AtlantECO-MAPS-v2_microbiome_monthly_surface_copepod_FD_2081-2100_Benedettietal.2025_20250604.nc", vars = var_list, force_v4 = TRUE)

### Fill new_nc with the data from the arrays stored in a list
for(fd_index in c(1:9)) {
    
    message(paste("Filling ",fd_names[fd_index]," values in the new NetCDF", sep = ""))
    
    # Subset array of interest
    arr <- fd_array_list[[fd_index]]
    
    for(s in seq_along(stat_names)) {

        stat <- stat_names[s]
        
        ncvar_put(
            new_nc,
            varid = stat,
            vals = arr[,,,s],
            start = c(1,1,1,fd_index),
            count = c(-1,-1,-1,1)
        ) # eo - ncvar_put()
        
     } # eo 2nd for loop - s

} # eo 1st for loop - fd_index
# print(new_nc)
# names(new_nc$var)
# lapply(new_nc$dim, function(x) list(name = x$name, len = x$len))


### Define Metadata for Each FD Indicator and putting them as attributes in new nc files

## Write metadata file
fd_metadata <- data.frame(
  
  name = fd_names,
  long_name = c(
    "Species Richness",
    "Faith's index for Functional Richness (Faith)",
    "Standardized Effect Size of Faith",
    "Functional Dispersion (FDis)",
    "Functional Divergence (FDiv)",
    "Functional Evenness (FEve)",
    "Trait Dissimilarity (Jaccard's dissimilarity index)",
    "Trait Turnover (turnover component of Jaccard's index)",
    "Trait Nestedness (nestedness component of Jaccard's index)"
  ),
  units = rep("unitless", 9),  # Adjust as needed
  description = c(
    "Potential number of species present per grid cell",
    "Assemblages with higher Faith values and more numerous the present species. Assemblages with are those where branches on the total volume filled by represent more distant functional dendrogram (i.e., more functional the assemblage)",
    "SES Faith values < 0 indicate that func-tional clustering (or functional convergence) occurs due to environmental filtering in the copepod assemblage whereas values > 0 indicate that there is functional overdispersion",
    "Assemblages with higher FDis values are those whose species are further away from each other and from the centroid in the functional space (i.e., more specialized species)",
    "Assemblages with higher FDiv values are characterized by higher HSI values at the vertices of their convex hull (i.e., more extreme traits values)",
    "Higher FEve values indicate that species in the assemblage display similar HSI at equal distances between nearest neighbors in the functional space. Lower FEve values indicate the coexistence of scattered clouds of functional units",
    "Trait dissimilarity values close to 1 indicate that two assemblages display functional dendrograms with very different number of branches that are non-overlapping",
    "Trait turnover values close to 1 indicate that total trait dissimilarity is driven by the replacement of branches",
    "Trait nestedness values close to 1 indicate that total trait dissimilarity is driven by different number of branches, whatever their identity"
  ), stringsAsFactors = FALSE
  
)

## Add FD indicator metadata to NetCDF
for(i in c(1:nrow(fd_metadata))) {
    
    # idx <- i
    
    ncatt_put(new_nc, varid = 0, 
            attname = paste0("FD_indicator_", i, "_name"), 
            attval = fd_metadata$name[i]
    )
  
    ncatt_put(new_nc, varid = 0, 
            attname = paste0("FD_indicator_", i, "_long_name"), 
            attval = fd_metadata$long_name[i]
    )
  
    ncatt_put(new_nc, varid = 0, 
            attname = paste0("FD_indicator_", i, "_description"), 
            attval = fd_metadata$description[i]
    )
  
    ncatt_put(new_nc, varid = 0, 
            attname = paste0("FD_indicator_", i, "_units"), 
            attval = fd_metadata$units[i]
    )
            
} # eo for loop


### Add global attributes to the new NetCDF
# if varid == 0, then a global attribute is written instead of a variable's attribute
ncatt_put(new_nc, varid = 0, attname = "Dataset", attval = "AtlantECO data - WP2 - Traditinal microscopy")
ncatt_put(new_nc, varid = 0,"Institution","ETH Zürich, D-USYS, IBP, UP group")
ncatt_put(new_nc, varid = 0,"Description","Mean/Median/Min/Max and Stdev of monthly projections of copepod species richness and eight functional diversity (FD) indicators for the future (2081-2100) global surface ocean predicted by an ensemble of three species distribution models (GLM, GAM and ANN) and five earth system models (CESM-BEC, CNRM-PISCES, GFDL-TOPAZ, IPSL-PISCES and MRI-NEMURO). See full description in Benedetti et al. (2025): Emergent Relationships Between the Functional Diversity of Marine Planktonic Copepods and Ecosystem Functioning in the Global Ocean")
ncatt_put(new_nc, varid = 0,"DOI","doi.org/10.1111/gcb.70094")
ncatt_put(new_nc, varid = 0, "Funding statement","This project has received funding from the European Union's Horizon 2020 Research and Innovation Programme under grant agreement no. 862923 (AtlantECO) and under grant agreement no. 101059915 (BIOceans5D). This output reflects only the author's view, and the European Union cannot be held responsible for any use that may be made of the information contained therein.")
# Even add history
history <- paste("Fabio Benedetti (fabio.benedetti@unibe.ch); last update: ", date(), sep = ", ")
ncatt_put(new_nc, varid = 0, "Last update by:", history)

### Close nc - should save the data 
nc_close(new_nc)

## Open the NetCDF file again to check
# nc <- nc_open("AtlantECO-MAPS-v2_microbiome_monthly_surface_copepod_FD_2081-2100_Benedettietal.2025_20250604.nc")
# print(nc)
# names(nc$var)
# lapply(nc$dim, function(x) list(name = x$name, len = x$len))
# nc_close(nc)


## Open as raster? 
library("raster")
# varname for the var of interest (Mean/median/etc.) and level for the month
r <- raster::brick("AtlantECO-MAPS-v2_microbiome_monthly_surface_copepod_FD_2081-2100_Benedettietal.2025_20250604.nc", varname = "Max", level = 4)
plot(r)
### Looks ok


### ----------------------------------------------------------------------------------------------------------------------------
### ----------------------------------------------------------------------------------------------------------------------------
### ----------------------------------------------------------------------------------------------------------------------------

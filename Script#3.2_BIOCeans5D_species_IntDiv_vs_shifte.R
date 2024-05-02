
##### Global cOpepoD functionaL diversitY (GODLY) project 

##### 20/11/23: R script to compare the functional distinctiveness index (IntDi) to the strength of the species' spatial shifts induced by Cl. Ch. © Fabio Benedetti, ETH Zürich, IBP, UP Group.

### Aims to:
# - explore covariance between IntDi and traits (continuous and categorical)
# - explore covariance between IntDi and amplitude of spatial shifts (results from 'RSCRIPTBATCH_BIOceans5D_PA_shifte.R') --> are the most functionally distinct the most affected by Cl. Change

### Latest update: 21/11/23

### ------------------------------------------------------------------------------------------------------------------------------------------------------

library("gawdis")
library("FD")
library("ape")
library("vegan")
library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("viridis")
library("ggpubr")
library("parallel")
library("xlsx")
library("readxl")
library("gtools")
library("naniar")

setwd("/net/kryo/work/fabioben/GODLY/data")

### ------------------------------------------------------------------------------------------------------------------------------------------------------

### 1°) Load and format functional traits table 
traits <- read.csv("traits_table_Benedetti2023.csv", h = T, sep = ";", dec = ",")
colnames(traits)[8] <- "Body.length" # size vector that we will keep 
# Replace "" by NA
traits <- traits %>% replace_with_na_all(condition = ~.x == "")
# Convert Feeding.mode, Trophic.group and Spawning.mode to factors
traits$Spawning.mode <- as.factor(traits$Spawning.mode)
traits$Trophic.group <- as.factor(traits$Trophic.group)
traits$Feeding.mode <- as.factor(traits$Feeding.mode)
# Count NA
names <- colnames(traits)[c(8:11,16)] ; names
traits$na_count <- apply(traits[,names], 1, function(x) sum(is.na(x)))
# Drop species with missing body length and more than two missing traits
traits_red <- traits[!is.na(traits$Body.length),]
traits_red2 <- traits_red[traits_red$na_count < 2,]
# Convert to 1 and 0 for Gower (as to be factors for FAMD+Eucli though)
traits_red2$Myelination <- as.integer(as.logical(traits_red2$Myelination))
traits_red2$Omnivore <- as.integer(as.logical(traits_red2$Omnivore))
traits_red2$Carnivore <- as.integer(as.logical(traits_red2$Carnivore))
traits_red2$Herbivore <- as.integer(as.logical(traits_red2$Herbivore))
traits_red2$Detritivore <- as.integer(as.logical(traits_red2$Detritivore))
traits_red2$Current <- as.integer(as.logical(traits_red2$Current))
traits_red2$Cruise <- as.integer(as.logical(traits_red2$Cruise))
traits_red2$Ambush <- as.integer(as.logical(traits_red2$Ambush))
traits_red2 <- data.frame(traits_red2)


### 2°) Load the IntDi values per species that you computed before

list_mat <- get(load("list_gawdis_mat_IntDi_24_07_23.Rdata"))

### Compute average matrix
# ?Reduce: ‘Reduce’ uses a binary function to successively combine the elements of a given vector and a possibly given initial value
mean_matrix_gaw <- Reduce("+",list_mat)/length(list_mat) # fast
# dim(mean_matrix_gaw) ; str(mean_matrix_gaw) # species pairwise dist matrix
# summary(mean_matrix_gaw)

### Compute average IntDi values
IntDi <- apply(mean_matrix_gaw, 1, sum, na.rm = T) / (nrow(mean_matrix_gaw) + 1)
IntDi <- cbind(as.data.frame(IntDi), traits_red2$Species)
colnames(IntDi)[2] <- "Species"
# dim(IntDi); head(IntDi); str(IntDi)
# summary(IntDi)
# save(IntDi, "IntDi_gaw_with_NA.Rdata")
# IntDi[IntDi$IntDi < .2,] # Why 0? 
# IntDi[IntDi$IntDi > .4,] # 
IntDi[order(IntDi$IntDi, decreasing = T),]

### Plot
IntDi$Species <- factor(IntDi$Species, levels = IntDi[order(IntDi$IntDi),"Species"]) # 
ggplot(IntDi[IntDi$IntDi > 0,], aes(Species, IntDi)) + geom_point() + ylab("Integrated distinctiveness (IntDi)") +
    scale_x_discrete(name = "Species", labels = NULL) + theme_classic() + theme(axis.ticks = element_blank()) 



### 3°) Explore relationships between IntDi and traits

# Attribute IntDi to traits table
# IntDi[150:155,"Species"] ; traits_red2[150:155,"Species"]
traits_red2$IntDi <- IntDi$IntDi

### ~ body length
ggplot() + geom_point(aes(x = Body.length, y = IntDi), data = traits_red2[traits_red2$IntDi > 0,]) +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Maximum body length (mm)") +
    theme_classic()
    
# cor(traits_red2$Body.length, traits_red2$IntDi, method = "spearman")
# 0.0328 # NS 

### ~ FG
ggplot(aes(y = IntDi, x = factor(FG)), data = traits_red2[traits_red2$IntDi > 0 & !is.na(traits_red2$FG),]) +
    geom_boxplot(colour = "black", fill = "gray") + ylab("Integrated distinctiveness (IntDi)") + xlab("FG (Benedetti et al., 2022)") +
    theme_classic()
### --> FGs that show highest IntDiv: 2,5,1,4...then 3,6.
unique(traits_red2[traits_red2$FG == 2,"Species"]) # Oncaeidae (only detritivores...)
unique(traits_red2[traits_red2$FG == 5,"Species"]) # Very large carnivores
unique(traits_red2[traits_red2$FG == 1,"Species"]) # Clausocalanus because of spawning...
unique(traits_red2[traits_red2$FG == 4,"Species"]) # Small ambush carnivores (Corycaeus)
# Most specialized taxa. Makes sense! Not your average copepod
### --> FGs that show LOWEST IntDiv: 9,10,11
unique(traits_red2[traits_red2$FG == 9,"Species"]) # small-medium sized omnivores
unique(traits_red2[traits_red2$FG == 10,"Species"]) # oithona
unique(traits_red2[traits_red2$FG == 11,"Species"]) # Acartia/Centropages

### ~ Order/Families
ggplot(aes(y = IntDi, x = factor(Family)), data = traits_red2[traits_red2$IntDi > 0,]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Family") + theme_classic() + theme(axis.text.x=element_text(angle = -90, hjust = 0))
# --> Euchaeta, Sapphirina, Oncaea

ggplot(aes(y = IntDi, x = factor(Order)), data = traits_red2[traits_red2$IntDi > 0,]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Order") + theme_classic() + theme(axis.text.x=element_text(angle = -90, hjust = 0))
# --> Calanoida LESS distinct than Cyclopoida (as expected)

### ~ Myelination
ggplot(aes(y = IntDi, x = factor(Myelination)), data = traits_red2[traits_red2$IntDi > 0,]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Myelination (1/0)") + theme_classic() 
# n.s. 

### ~ Spawning strat
ggplot(aes(y = IntDi, x = factor(Spawning.mode)), data = traits_red2[traits_red2$IntDi > 0,]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Spawning mode") + theme_classic() 
# --> sac-spawners are more distinct

### ~ Trophic group
ggplot(aes(y = IntDi, x = factor(Trophic.group)), data = traits_red2[traits_red2$IntDi > 0,]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Trophic group") + theme_classic() 
# --> Detritivore and pure carnivores > distinct

### ~ Feeding.mode
ggplot(aes(y = IntDi, x = factor(Feeding.mode)), data = traits_red2[traits_red2$IntDi > 0,]) + geom_boxplot(colour = "black", fill = "gray") +
    ylab("Integrated distinctiveness (IntDi)") + xlab("Feeding mode") + theme_classic() 
# --> cruise-feeders and more distinct


### 4°) ~ 4 PCoA axes! (corr, plots)
mat.gow <- gawdis(x = traits_red2[,c(8,9,10,12:15,17:19)], groups = c(1,2,3,4,4,4,4,5,5,5), fuzzy = c(4,5))
pcoa2 <- wcmdscale(d = mat.gow, eig = T)
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

# Test correlations
# round(cor(IntDi[IntDi$IntDi > 0,"IntDi"], IntDi[IntDi$IntDi > 0,"PCoA1"], method = "spearman"),3) # -0.279
# round(cor(IntDi[IntDi$IntDi > 0,"IntDi"], IntDi[IntDi$IntDi > 0,"PCoA2"], method = "spearman"),3) # 0.669
# round(cor(IntDi[IntDi$IntDi > 0,"IntDi"], IntDi[IntDi$IntDi > 0,"PCoA3"], method = "spearman"),3) # -0.053 n.s.
# round(cor(IntDi[IntDi$IntDi > 0,"IntDi"], IntDi[IntDi$IntDi > 0,"PCoA4"], method = "spearman"),3) # 0.261
# PCoA 2 seems to be the main axis of variation

# Plot in PCoA space
p1 <- ggplot(aes(x = PCoA1, y = PCoA2), data = IntDi[IntDi$IntDi > 0,]) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..)), fill = "grey") + 
    geom_point(aes(fill = IntDi), pch = 21, colour = "black") +
    scale_fill_viridis(name = "IntDiv", option = "D", direction = 1) + 
    scale_alpha_continuous(range = c(0, 1)) + 
    theme_bw() + xlab(scores1) + ylab(scores2)+ guides(alpha = "none")

p2 <- ggplot(aes(x = PCoA3, y = PCoA4), data = IntDi[IntDi$IntDi > 0,]) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..)), fill = "grey") + 
    geom_point(aes(fill = IntDi), pch = 21, colour = "black") + 
    scale_fill_viridis(name = "IntDiv", option = "D", direction = 1) + 
    scale_alpha_continuous(range = c(0, 1)) + 
    theme_bw() + xlab(scores3) + ylab(scores4)+ guides(alpha = "none")

ggarrange(p1,p2, labels = letters, align = "hv", ncol = 2, nrow = 1)


### Show taxa with highest PCoA2 scores
IntDi[order(IntDi$PCoA2, decreasing = T),"Species"]


### ----------------------------------------------------------------------

### 21/11/23: 5°) Combine with ensemble estimates of shifts in spatial centroids (between future and contemp. communities) and compare to IntDi index 
### Aims to answer the following question: Are more functionally distinct species/FGs disproportionally affected by Cl. Ch. ? Relate to Benedetti et al. (2018b)
### You may also plot the shift metrics into the PCoA space etc. 

library("geosphere")

### But first, load centroids position and compute shifts in mean annual centroids positions for each SDMxESM combin
setwd("/net/kryo/work/fabioben/GODLY/data/shift_metrics")
res <- lapply(dir()[grep("contemp",dir())], function(f) { d <- get(load(f)) ; return(d) })
base <- dplyr::bind_rows(res); rm(res); gc()
# head(base) 

# For the 'global' centroids, compute each future shift relative to the SDM-specific contemporary positions
# f <- dir()[grepl("_fut",dir())][1]
res <- lapply(dir()[grepl("fut",dir())], function(f) {
        
        message(paste("\n","Compute shifts in mean annual centroids for ",f,"\n", sep = ""))
        d <- get(load(f)) 
        # Choose reference of interest wihtin 'base'
        centroids.base <- base[which(base$SDM == unique(d$SDM) & base$region == "Global"),]
        
        # For each species, compute spatial shifts
        spp <- unique(d$species) 
        # s <- "Acartia_Acartiura_clausi" # for testing
        res2 <- mclapply(spp, function(s) {
            
                    message(paste("Computing shifts for ",s, sep = ""))
                    ref <- centroids.base[centroids.base$species == s,]
                    fut <- d[d$species == s & d$region == "Global",]
                        
                    dist <- distm(ref[,c("lon","lat")], fut[,c("lon","lat")], fun = distHaversine)
                    # Convert to km
                    dist.km <- dist/1000
                    
                    return(data.frame(species = s, shift.km = dist.km))
            
                }, mc.cores = 20
                
        ) # eo mclapply
        shifts <- dplyr::bind_rows(res2); rm(res2); gc()
        shifts$SDM <- unique(d$SDM)
        shifts$ESM <- unique(d$ESM)
        
        # Return
        return(shifts)
    }
         
) # eo lapply
shifts <- dplyr::bind_rows(res)
rm(res); gc()
# dim(shifts); summary(shifts)

### Summarize ensemble prediction
ens <- data.frame( shifts %>% group_by(species) %>% summarise(mean.shift = mean(shift.km), sd = sd(shift.km)) )
# summary(ens)
ens$species <- as.character(ens$species)

### Supply to 'traits_red2' and analyze covariance with IntDi
# unique(ens$species); unique(traits_red2$Species) # Uneven species list, restrict to commons
traits_red2$mean.shift <- NA
traits_red2$sd.shift <- NA

for(s in unique(ens$species)) {
    message(s)
    traits_red2[traits_red2$Species == s,"mean.shift"] <- ens[ens$species == s,"mean.shift"]
    traits_red2[traits_red2$Species == s,"sd.shift"] <- ens[ens$species == s,"sd"]
} # eo for loop - s in species

# Test covariance between IntDi and mean shifts amplitude
# ggplot() + geom_point(aes(x = traits_red2$IntDi, y = traits_red2$mean.shift), data = traits_red2) + theme_bw()
# No relationship
df <- traits_red2 %>% filter(!(is.na(mean.shift)))
# cor(df$IntDi, df$mean.shift)
# summary(lm(IntDi ~ mean.shift, data = df))

traits_red2$Species <- factor(traits_red2$Species, levels = traits_red2[order(traits_red2$IntDi),"Species"]) # 

# ggplot(traits_red2[traits_red2$IntDi > 0,], aes(x = Species, y = IntDi, fill = mean.shift, size = mean.shift)) +
#     geom_point(pch = 21, colour = "black") + ylab("Integrated distinctiveness (IntDi)") +
#     scale_fill_distiller(palette = "YlOrRd") + scale_x_discrete(name = "Species", labels = NULL) +
#     theme_classic() + theme(axis.ticks = element_blank())

# No relationship. Plot distrib of shifts per FG (boxplots)
ggplot(aes(y = mean.shift, x = factor(FG)), data = traits_red2[!is.na(traits_red2$FG),]) +
    geom_boxplot(colour = "black", fill = "gray") + ylab("Centroid shift (km)") + xlab("FG (Benedetti et al., 2023)") +
    theme_classic()
### --> all medians kind of overlap with IQR; no strong variance in shifts per known FG
    
### Explore distrbution of mean shifts in PCoA space
IntDi$mean.shift <- NA
IntDi$sd.shift <- NA
for(s in unique(ens$species)) {
    message(s)
    IntDi[IntDi$Species == s,"mean.shift"] <- ens[ens$species == s,"mean.shift"]
    IntDi[IntDi$Species == s,"sd.shift"] <- ens[ens$species == s,"sd"]
} # eo for loop - s in species
# summary(IntDi) # remove NAs
df <- IntDi %>% filter(!(is.na(mean.shift)))

# round(cor(df[df$IntDi > 0,"mean.shift"], df[df$IntDi > 0,"PCoA1"], method = "spearman"),3) # 0.076 n.s.
# round(cor(df[df$IntDi > 0,"mean.shift"], df[df$IntDi > 0,"PCoA2"], method = "spearman"),3) # -0.007 n.s.
# round(cor(df[df$IntDi > 0,"mean.shift"], df[df$IntDi > 0,"PCoA3"], method = "spearman"),3) # 0.073 n.s.
# round(cor(df[df$IntDi > 0,"mean.shift"], df[df$IntDi > 0,"PCoA4"], method = "spearman"),3) # -0.101 forget it

p1 <- ggplot(aes(x = PCoA1, y = PCoA2), data = df[df$IntDi > 0,]) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..)), fill = "grey") + 
    geom_point(aes(fill = mean.shift, size = mean.shift), pch = 21, colour = "black") +
    scale_fill_viridis(name = "Shift in mean annual\ncentroids (km)", option = "D", direction = 1) + 
    scale_alpha_continuous(range = c(0, 1)) + 
    theme_bw() + xlab(scores1) + ylab(scores2)+ guides(alpha = "none")

p2 <- ggplot(aes(x = PCoA3, y = PCoA4), data = df[df$IntDi > 0,]) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    stat_density_2d(geom = "polygon", aes(alpha = (..level..)), fill = "grey") + 
    geom_point(aes(fill = mean.shift, size = mean.shift), pch = 21, colour = "black") +
    scale_fill_viridis(name = "Shift in mean annual\ncentroids (km)", option = "D", direction = 1) + 
    scale_alpha_continuous(range = c(0, 1)) + 
    theme_bw() + xlab(scores3) + ylab(scores4)+ guides(alpha = "none")

ggarrange(p1,p2, labels = letters, align = "hv", ncol = 2, nrow = 1)
### No real pattern, as expected

### Examine distrbution of shifts...seems like there is an outlier
ggplot(df[df$IntDi > 0,], aes(x = mean.shift)) + geom_histogram(binwidth=50, colour="black", fill="white")
# Who are those 3-5 species with shifts > 1000km?
ens[order(ens$mean.shift, decreasing = T),"species"]
#   [1] "Spinocalanus_elongatus"             "Parvocalanus_latus"                
#   [3] "Acartia_Odontacartia_amboinensis"   "Scolecithricella_orientalis"       
#   [5] "Atrophia_glacialis"

### Examine covariance between strength of mean shift and sd
ggplot(aes(x = mean.shift, y = sd), data = ens) + geom_point() + geom_smooth(method = "lm", se = T) + theme_classic()
### Obvious and very strong linear trend between amplitude of mean shift and its uncertainty
### Re-run tests above but only with those species that shift by < 1000 km
df2 <- df[df$IntDi > 0 & df$mean.shift < 1000,]
# dim(df2)
# cor(df2$mean.shift, df2$IntDi, method = "spearman") # -0.050 n.s.

# Per PCoA axes:
# round(cor(df2[,"mean.shift"], df2[,"PCoA1"], method = "spearman"),3) # 0.069 n.s.
# round(cor(df2[,"mean.shift"], df2[,"PCoA2"], method = "spearman"),3) # -0.004 n.s.
# round(cor(df2[,"mean.shift"], df2[,"PCoA3"], method = "spearman"),3) # 0.062 n.s.
# round(cor(df2[,"mean.shift"], df2[,"PCoA4"], method = "spearman"),3) # -0.103 n.s.

# ggplot(df2, aes(x = Species, y = IntDi, fill = mean.shift, size = mean.shift)) +
#     geom_point(pch = 21, colour = "black") + ylab("Integrated distinctiveness (IntDi)") +
#     scale_fill_distiller(palette = "YlOrRd", direction = 1) +
#     scale_x_discrete(name = "Species", labels = NULL) +
#     theme_classic() + theme(axis.ticks = element_blank())

### Nope. That's not it. No pattern.
### CCL: cl ch. does not impact most distinct FG and species disproportionally.

### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------------------------------------------------

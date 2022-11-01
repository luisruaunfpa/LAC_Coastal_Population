# UNFPA LACRO - LAC Coastal Population Assessment #
# Luis de la Rua - July 2022 #

# OVERVIEW: Estimate coastal population in LAC region using population grids and 1, 5, 10km distance to coast buffers
# Different approach - Split analysis by country pieces to better handle data processing load.

# RESULT: Population count and % by country and 1,5,10Km buffers in a single table
# OUTPUT LAYERS: Buffers by country, Projected population grids into equal area projection 54034.

#for Southern Lat countries check difference between running analysis on 4326 or 54034

# LIBRARIES ---------------------------------------------------------------
library(sf)
library(sp)
library(maptools)
library(raster)
library(spData)
library(rgeos)
library(dplyr)
library(rgdal)
library(magrittr)
library(terra)
library(data.table)
library(leaflet)

#SET PARAMETERS-----------------------------------------------------------------
#working directory (this hast to be changed in case we are working from shared folder)
wd <- "G:/.shortcut-targets-by-id/1i7Dq5-HfwQYp5ApFZZAw3pZb8GM0S76S/UNFPA Luis de la Rua/Coastal_Population_LAC" #Luis local  drive
setwd(wd)

# blocks viewer screens and forces the maps to be displayed on browser
options(viewer = NULL)
# Set Equal Area projection
ea_crs<- crs('+proj=cea +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')


### ABW  -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'ABW' ###NEEDS TO BE CHANGED

# admin boundary
#ab_path <- (paste0("layers/ad_bound/",country,"/","AIA_adm0.shp"))  ###NEEDS TO BE CHANGED

#

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')

## WITH SP LIBRARY

# Project option points to raster > project points > raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize
pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### AIA  -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'AIA' ###NEEDS TO BE CHANGED

# admin boundary
# ab_path <- (paste0("layers/ad_bound/",country,"/","AIA_adm0.shp"))  ###NEEDS TO BE CHANGED
# 
# ab <- st_read(ab_path) #dont knonw yet if it is necessary. (We need it to create Buffers but this is done in QGIS)
# ab_ea <- st_transform(ab,ea_crs)

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v1.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')

## WITH SP LIBRARY

# Project option points to raster > project points > raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize
pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### ARG 1km res-----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'ARG' ###NEEDS TO BE CHANGED

# admin boundary
# ab_path <- (paste0("layers/ad_bound/",country,"/","AIA_adm0.shp"))  ###NEEDS TO BE CHANGED
# 
# ab <- st_read(ab_path) #dont knonw yet if it is necessary. (We need it to create Buffers but this is done in QGIS)
# ab_ea <- st_transform(ab,ea_crs)

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v1.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project points > raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
# For Big Countries we change raster resoultion to 1spkm
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 1000
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize
pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
#writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### ATG  -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'ATG' ###NEEDS TO BE CHANGED

# admin boundary
#ab_path <- (paste0("layers/ad_bound/",country,"/","AIA_adm0.shp"))  ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg,'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project points > raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize
pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### BHS  -----------------------------------------------------------------------

# BRING DATA IN ## need to change in every Country section
country <- 'BHS' ### NEEDS TO BE CHANGED

# admin boundary
ab_path <- (paste0("layers/ad_bound/",country,"/","BHS_adm0.shp"))  ###NEEDS TO BE CHANGED

ab <- st_read(ab_path) #dont knonw yet if it is necessary. (We need it to create Buffers but this is done in QGIS)
ab_ea <- st_transform(ab,ea_crs)

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))   ###NEEDS TO BE CHANGED Version numbers tend to differ
pg_path
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')

## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize
pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### BLZ  -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'BLZ' ###NEEDS TO BE CHANGED

# admin boundary
ab_path <- (paste0("layers/ad_bound/",country,"/","BLZ_adm0.shp"))

ab <- st_read(ab_path) #dont knonw yet if it is necessary. (We need it to create Buffers but this is done in QGIS)
ab_ea <- st_transform(ab,ea_crs)

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_BLZ_v1.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')

## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### BMU  -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'BMU' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_BMU_v3.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### BOL NO SEA BOUNDARIES

### BRA 1km res-----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'BRA' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg #calculate total population from whole raster

# for raster analysis use raster clipped by b10 extent
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2_cliped.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)

## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 1000 #big country 1km resolution otherwise too long process
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### BRB  -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'BRB' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### BOL NO SEA BOUNDARIES
### CRI  -----------------------------------------------------------------------

# BRING DATA IN ## need to change in every Country section
country <- 'CRI' ### NEEDS TO BE CHANGED

# admin boundary
ab_path <- (paste0("layers/ad_bound/",country,"/","cri_admbnda_adm0_2021_WGS_84.shp"))  ###NEEDS TO BE CHANGED

ab <- st_read(ab_path) #dont knonw yet if it is necessary. (We need it to create Buffers but this is done in QGIS)
ab_ea <- st_transform(ab,ea_crs)

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))   ###NEEDS TO BE CHANGED Version numbers tend to differ
pg_path
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')

## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize
pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### BOL NO SEA BOUNDARIES
### CHL 1km res-----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'CHL' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 1000 #big country 1km resolution otherwise too long process
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### COL 1km res-----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'COL' ###NEEDS TO BE CHANGED

#bring buffer using st (see possibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 1000 #big country 1km resolution otherwise too long process
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### CUB -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'CUB' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))


### CUW -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'CUW' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))


### CYM -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'CYM' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))


### DMA -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'DMA' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))


### DOM -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'DOM' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))


### ECU 1km res-----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'ECU' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 1000 #big country 1km resolution otherwise too long process
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### GLP 1km res-----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'GLP' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_2020.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 1000 #big country 1km resolution otherwise too long process
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### GRD -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'GRD' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))


### GTM -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'GTM' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))


### GUF -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'GUF' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","guf_ppp_2020.tif"))   ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))


### GUY -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'GUY' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v1.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### HND -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'HND' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### HTI -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'HTI' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","hti_ppp_2020.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### JAM -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'JAM' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### KNA -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'KNA' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### LCA -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'LCA' ###NEEDS TO BE CHANGED

# admin boundary

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### MEX 1km res-----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'MEX' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 1000 #big country 1km resolution otherwise too long process
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### MSR -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'MSR' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=7,byrow=T)
colnames(output) <- c('country','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
# BRING DATA IN ## need to change in every Country section
country <- 'MSR' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))


### MTQ -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'MTQ' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_2020.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### NIC -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'NIC' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### PAN -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'PAN' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### PER -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'PER' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))

### PRI -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'PRI' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_2020.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### SLV -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'SLV' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### SUR -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'SUR' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### TCA -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'TCA' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### TTO -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'TTO' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### URY -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'URY' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### VCT -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'VCT' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v1.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### VEN 1km res-----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'VEN' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v2.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 1000 #big country 1km resolution otherwise too long process
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
### VGB -----------------------------------------------------------------------
# BRING DATA IN ## need to change in every Country section
country <- 'VGB' ###NEEDS TO BE CHANGED

#bring buffer using st (see posibility of generating buffer from line on 4326) here below so far buffer generated in QGIS
b1_path <- (paste0('layers/buffer/',country,"/",country,"_buff1_54034.gpkg"))    ###NEEDS TO BE CHANGED
b5_path <- (paste0('layers/buffer/',country,"/",country,"_buff5_54034.gpkg"))    ###NEEDS TO BE CHANGED
b10_path <- (paste0('layers/buffer/',country,"/",country,"_buff10_54034.gpkg"))  ###NEEDS TO BE CHANGED

b1 <- st_read(b1_path)
b5 <- st_read(b5_path)
b10 <- st_read(b10_path)

# Pop grid
pg_path <- (paste0("layers/pop_grid/",country,"/","ppp_",country,"_v3.tif"))  ###NEEDS TO BE CHANGED
pg  <- raster(pg_path)
sum_pg <- cellStats(pg, 'sum')
sum_pg
## WITH SP LIBRARY

# Project option points to raster > project >raster to points 

# Raster to point
pg_pts <- rasterToPoints(pg,spatial=T)
names(pg_pts)[1]="pop" ### Assigning same variable name to all pts layers

# Reproject Equal Area CRS
pg_pts_ea <- spTransform(pg_pts,ea_crs)

# create blank raster to embed point data 
pg_ea <- raster()
extent(pg_ea)<- extent(pg_pts_ea) # same as point layer
res(pg_ea) <- 100
crs(pg_ea) <- crs(pg_pts_ea)
pg_ea

# Rasterize

pg_ea <- rasterize(pg_pts_ea,pg_ea,'pop',fun=sum)
writeRaster(pg_ea,paste0("layers/pop_grid/",country,"/",country,"_pg_ea.tif"),overwrite=T) ### OUTPUT

# Check sum population for both points and raster
sum(pg_pts_ea$pop) 
cellStats(pg_ea,'sum')

# Mask raster using different buffer vectors
b1mask <- mask(pg_ea,b1)
b5mask <- mask(pg_ea,b5)
b10mask <- mask(pg_ea,b10)

# Calculating population within buffers
pop_b1 <- cellStats(b1mask,'sum')
pop_b5 <- cellStats(b5mask,'sum')
pop_b10 <- cellStats(b10mask,'sum')

# Calculating percentages
per_pop_b1 <- (pop_b1/sum_pg)
per_pop_b5 <- (pop_b5/sum_pg)
per_pop_b10 <- (pop_b10/sum_pg)

output <- matrix(c(country,sum_pg,pop_b1,pop_b5,pop_b10,per_pop_b1,per_pop_b5,per_pop_b10),ncol=8,byrow=T)
colnames(output) <- c('country','tpop','pop_b1','pop_b5','pop_b10','per_pop_b1','per_pop_b5','per_pop_b10')
write.csv(output,paste0('result/',country,'_output.csv'))
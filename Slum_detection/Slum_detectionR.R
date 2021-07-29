#########################
# Extracting slums from satellite imagery. R
# Francina Cabrera, 2021

# Description:

# Data:

# Code source: 
# https://pages.cms.hu-berlin.de/EOL/gcg_eo/index.html
# https://andrewmaclachlan.github.io/CASA0005repo/advanced-raster-analysis.html
# http://www.qgistutorials.com/en/docs/3/raster_mosaicing_and_clipping.html

################################################################################

# Load all the required packages
#install.packages('')

library(sp)
library(raster)
library(rgdal)
library(sf)
library(here)
library(rgeos)
library(rasterVis)
library(ggplot2)
library(stringr)
library(fs)
library(tidyverse)
library(maptools)
library(tmap)
library(GGally)
library(hexbin)
library(reshape2)
library(randomForest)
library(caTools)
library(googledrive)
library(utils)

################################################################################
# 1. Load and prepare the data
# a.1) Distrito Nacional's boundary
#dn_boundary <- st_read(here::here("data","geo","DN_boundary.shp")) # Test
dn_boundary <- st_read(here("data","geo","DN_boundary.shp")) %>% 
  # Project to EPSG:32619. This is the projected coordinate system for the Dominican Republic.
  st_transform(.,32619)
# Quick look of the feature
qtm(dn_boundary)

# a.2) Circumscription 3
c3_boundary <- st_read(here("data","geo","Circ3.shp")) %>% 
  # Project to EPSG:32619. This is the projected coordinate system for the Dominican Republic.
  st_transform(.,32619)
# Quick look of the feature
qtm(c3_boundary)

# b) PlanetScope scenes from Distrito Nacional taken on the 23/02/2021.
# The images will be downloaded from Google Drive.
# First scene:
# Create a temporary file that will be substituted with the downloaded raster.
temp49 <- tempfile(fileext = ".tif")
# Download the raster with the file's id from Google Drive.
dl49 <- drive_download(as_id("1ftQXEszRYiHfiCDkm1kJ3parha85569w"),
                     path = temp49, overwrite = TRUE)
# Create the raster stack
raster49 <- stack(temp49)
# Check the coordinate system is WGS 84 / UTM zone 19N,
# if not use: %>% projectRaster(., crs=32619)
crs(raster49) 

# Second scene:
# Create a temporary file that will be substituted with the downloaded raster.
temp50 <- tempfile(fileext = ".tif")
# Download the raster with the file's id from Google Drive.
dl50 <- drive_download(as_id("1But-SQaZkdazRBNMUKXjocLyiD7GQkJL"),
                       path = temp50, overwrite = TRUE)
# Create the raster stack
raster50 <- stack(temp50)
# Check the coordinate system is WGS 84 / UTM zone 19N,
# if not use: %>% projectRaster(., crs=32619)
crs(raster50) 

# Quick look at the images
plot(raster49, main="PlanetScope Distrito Nacional 20210223_123349_104e_3B_AnalyticMS_SR")
plot(raster50, main="PlanetScope Distrito Nacional 20210223_123350_104e_3B_AnalyticMS_SR")

# Merge the images to create a mosaic using the function mosaic()
raster_mosaic <- mosaic(raster49, raster50, fun=mean, tolerance=0.05, filename="raster_mosaic", overwrite=TRUE)

# Quick look at the mosaics
plot(raster_mosaic, main="Mosaic of PlanetScope scenes")

# Crop the mosaic to the Distrito Nacional's boundary with crop()
mosaic_DN <- crop(raster_mosaic, dn_boundary, filename="mosaicDN_crop", overwrite=TRUE) %>% 
  mask(., dn_boundary)

# Crop the mosaic to the Circunscription 3 (C3)'s boundary with crop()
mosaic_C3 <- crop(raster_mosaic, c3_boundary, filename="mosaicC3_crop", overwrite=TRUE) %>% 
  mask(., c3_boundary)

# Quick look at the cropped raster
plot(mosaic_DN, main="Mosaic cropped to Distrito Nacional's boundary")
plot(mosaic_C3, main="Mosaic cropped to Circunscription 3's boundary")

# Information on the Distrito Nacional's raster stack
extent(mosaic_DN) # extent
ncell(mosaic_DN) # number of cells
dim(mosaic_DN) # number of rows, columns, layers
nlayers(mosaic_DN) # number of layers
res(mosaic_DN) # xres, yres

# Information on the Circunscription 3 (C3)'s raster stack 
extent(mosaic_C3) # extent
ncell(mosaic_C3) # number of cells
dim(mosaic_C3) # number of rows, columns, layers
nlayers(mosaic_C3) # number of layers
res(mosaic_C3) # xres, yres

# Name the Bands based on where they sample the electromagnetic spectrum
# DN's mosaic
names(mosaic_DN) <- c('blue', 'green', 'red', 'NIR') 
# C3' mosaic
names(mosaic_C3) <- c('blue', 'green', 'red', 'NIR')

# Plot the data in true colours and false composite.
# true colour composite
rC_rgb <- stack(mosaic_C3$red, mosaic_C3$green, mosaic_C3$blue) %>% 
  plotRGB(.,axes=TRUE, stretch="lin")
# false colour composite
rC_false <- stack(mosaic_C3$NIR, mosaic_C3$red, mosaic_C3$green) %>% 
  plotRGB(.,axes=TRUE, stretch="lin")

# Check the similarity between bands
# set the plot window size (2 by 2)
par(mfrow = c(2,2))
# plot the bands - Mosaic C
plot(mosaic_C3$blue, main = "Blue")
plot(mosaic_C3$green, main = "Green")
plot(mosaic_C3$red, main = "Red")
plot(mosaic_C3$NIR, main = "NIR")

# Statistics of these bands
# Mosaic C3
mosaic_C3 %>%
  as.data.frame(., na.rm=TRUE) %>%
  sample_n(., 100) %>%
  ggpairs(.,axisLabels="none")

# Training data 2
# The training data was collected digitally using a very high resolution (VHR) imagery as a reference: the Google maps terrain in QGIS.
# The location of the informal settlements in the Distrito Nacional was obtained from: http://adn.gob.do/joomlatools-files/docman-files/borrador_plan_est/Borrador%20Plan%20Estrategico%20del%20Distrito%20Nacional%202030%20%20%20V.%2028%20JUL%202020.pdf (page 75)
# Following the classification for informal and formal settlements from: https://ieeexplore.ieee.org/document/6236225,
# the classification follows six classes:
# 1 - Informal settlements Type I: defined roads
# 2 - Informal settlements Type II: undefined roads
# 3 - Informal settlements Type II: undefined roads and located in hazardous areas
# 4 - Formal settlements (all built up areas)
# 5 - Non-settlemets (vegetation, bare ground, water, grassland, parking lots, sport courts)
# 6 - Roads

# Read training points, the following code assumes that it contains only the class attribute
# in readOGR, dsn specifies the path to the folder containing the file (may not end with /), 
# layer specifies the name of the shapefile without extension (.shp)
trainC3 <- st_read(here("data", "training_data","TrainingData_C3.shp")) %>% 
  # Project to EPSG:32619.
  st_transform(.,32619)
# Quick look of the feature
qtm(trainC3)
crs(trainC3) # Check crs

























##

# temp10 <- tempfile(fileext = ".shp")
# download.file("https://github.com/francicabrera/DI_Francina2021/blob/main/Slum_detection/data/geo/DN_boundary.shp",
#               destfile = "/Users/francinacabrera/Documents/test.shp")
# 
# dn_bound <- st_read(here("Users","francinacabrera","Documents","test.shp"))
# 
# 
# 
# 
# https://raw.githubusercontent.com/knightlab-analyses/qurro-mackerel-analysis/master/AnalysisOutput/qurro-plot.qzv
# 
# 
# dn_bound <- st_read(here("var","folders","l2","z4t104ms2ln5m4dw5m6dycrm0000gn","T","RtmpLQcmD9","filed24f32c9914f.shp"))
# 
# https://raw.githubusercontent.com/francicabrera/DI_Francina2021/main/Slum_detection/data/geo/DN_boundary.shp?raw=true
# 
# 
# test <- readOGR("https://github.com/francicabrera/DI_Francina2021/blob/main/Slum_detection/data/geo/DN_boundary.shp?raw=true")
# 
# 
# test <- readOGR("https://github.com/francicabrera/DI_Francina2021/blob/main/Slum_detection/data/geo/DN_boundary.shp")
# 
# 
# 
# ##
# plot(raster01, main="PlanetScope Distrito Nacional 20210223_123349_104e_3B_AnalyticMS_SR")
# 
# 
# listfiles<-dir_info(here::here("CASA","Dissertation","DIS_FRANCINA_2021","Data","raster")) %>%
#   dplyr::filter(str_detect(path, ".tif")) %>%
#   dplyr::select(path)%>%
#   dplyr::pull()%>%
#   #print out the .gz file
#   as.character()%>%
#   utils::untar(exdir=here::here("CASA","Dissertation","DIS_FRANCINA_2021","Data","raster"))
# 
# 
# 
# 
# 
# 
# salida <- utils::untar(exdir=here::here("prac7_data", "exampleGoogleDrivedata"))
# 
# out <- unzip(temp, exdir = tempdir()) %>% 
#   stack()
# 
# bank <- read.csv(out[14], sep = ";")
# 
# 
# library(curl)
# id <- "1ftQXEszRYiHfiCDkm1kJ3parha85569w"
# URL <- sprintf("https://drive.google.com/drive/folders/%s", id)
# con <- curl(URL)
# 
# myUrl <- "https://drive.google.com/file/d/1ftQXEszRYiHfiCDkm1kJ3parha85569w/view?usp=sharing"
# test <- raster(url(myUrl))
# 
# 
# ###
# 
# o<-drive_download("https://drive.google.com/file/d/1ftQXEszRYiHfiCDkm1kJ3parha85569w/view?usp=sharing",
#                   path="CASA/Dissertation/DIS_FRANCINA_2021/Data/raster/20210223_123349_104e_3B_AnalyticMS_SR.tif", 
#                   overwrite=T)
# 
# ####
# dl2 <- drive_download(as_id("1ftQXEszRYiHfiCDkm1kJ3parha85569w"),
#                      path = "CASA/Dissertation/DIS_FRANCINA_2021/Data/raster/20210223_123349_104e_3B_AnalyticMS_SR.tif", overwrite = TRUE)
# 
# #########
# 
# # 1. Load and prepare the data
# # a.1) Distrito Nacional's boundary
# #dn_boundary <- st_read(here::here("data","geo","DN_boundary.shp")) # Test
# dn_boundary <- st_read(here("data","geo","DN_boundary.shp")) %>% 
#   # Project to EPSG:32619. This is the projected coordinate system for the Dominican Republic.
#   st_transform(.,32619)
# # Quick look of the feature
# qtm(dn_boundary)
# 
# # a.2) Circumscription 3
# c3_boundary <- st_read(here("data","geo","Circ3.shp")) %>% 
#   # Project to EPSG:32619. This is the projected coordinate system for the Dominican Republic.
#   st_transform(.,32619)
# # Quick look of the feature
# qtm(c3_boundary)
# 
# # b) PlanetScope scenes from Distrito Nacional taken on the 23/02/2021 (Option A)
# o <- drive_download("https://drive.google.com/file/d/1ftQXEszRYiHfiCDkm1kJ3parha85569w/view?usp=sharing",
#                     path = "CASA/Dissertation/DIS_FRANCINA_2021/Data/raster/20210223_123349_104e_3B_AnalyticMS_SR.tif", overwrite = T) 
# raster49 <- dir_info(here::here("CASA","Dissertation","DIS_FRANCINA_2021","Data","raster")) %>% 
#   dplyr::filter(str_detect(path, "20210223_123349_104e_3B_AnalyticMS_SR.tif")) %>% 
#   dplyr::select(path) %>% 
#   pull() %>% 
#   as.character() %>% 
#   # Load the raster layers into a stack
#   stack()
# 
# o <- drive_download("https://drive.google.com/file/d/16MK4QWC33wKaraciE8GNXIQOG1eElW1a/view?usp=sharing",
#                     path = "CASA/Dissertation/DIS_FRANCINA_2021/Slum_detection/data/SantoDominfo_FEB2021_psscene4band_analytic_sr_udm2/files/PSScene4Band/20210223_123349_104e/analytic_sr_udm2", overwrite = T) 
# raster49 <- dir_info(here::here("CASA","Dissertation","DIS_FRANCINA_2021","Slum_detection","data","SantoDominfo_FEB2021_psscene4band_analytic_sr_udm2","files","PSScene4Band","20210223_123349_104e","analytic_sr_udm2")) %>% 
#   dplyr::filter(str_detect(path, "20210223_123349_104e_3B_AnalyticMS_SR.tif")) %>% 
#   dplyr::select(path) %>% 
#   pull() %>% 
#   as.character() %>% 
#   # Load the raster layers into a stack
#   stack()
# 
# 
# 
# 
# o<-drive_download("https://drive.google.com/open?id=1MV7ym_LW3Pz3MxHrk-qErN1c_nR0NWXy",
#                   path="prac7_data/exampleGoogleDrivedata/LC08_L1TP_203023_20190513_20190521_01_T1.tar.gz", 
#                   overwrite=T)
# 
# listfiles<-dir_info(here::here("prac7_data", "exampleGoogleDrivedata")) %>%
#   dplyr::filter(str_detect(path, ".gz")) %>%
#   dplyr::select(path)%>%
#   dplyr::pull()%>%
#   #print out the .gz file
#   print()%>%
#   as.character()%>%
#   utils::untar(exdir=here::here("prac7_data", "exampleGoogleDrivedata"))
# 
# # List your raster files excluding band 8 using the patter argument
# listlandsat<-dir_info(here::here("prac7_data", "Lsatdata"))%>%
#   dplyr::filter(str_detect(path, "[B123456790].TIF")) %>%
#   dplyr::select(path)%>%
#   pull()%>%
#   as.character()%>%
#   # Load our raster layers into a stack
#   stack()
# 
#   
# 
# 
# 
# 
# 

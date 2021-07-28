#########################
# Extracting slums from satellite imagery. R
# Francina Cabrera, 2021
# Description:

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

# b) PlanetScope scenes from Distrito Nacional taken on the 23/02/2021 (Option A)
o <- drive_download("",
                    path = "CASA/Dissertation/DIS_FRANCINA_2021/Data/raster", overwrite = T) 
raster49 <- dir_info(here::here("CASA","Dissertation","DIS_FRANCINA_2021","Data","raster")) %>% 
  dplyr::filter(str_detect(path, "20210223_123349_104e_3B_AnalyticMS_SR.tif")) %>% 
  dplyr::select(path) %>% 
  pull() %>% 
  as.character() %>% 
  # Load the raster layers into a stack
  stack()

o <- drive_download("https://drive.google.com/file/d/16MK4QWC33wKaraciE8GNXIQOG1eElW1a/view?usp=sharing",
                    path = "CASA/Dissertation/DIS_FRANCINA_2021/Slum_detection/data/SantoDominfo_FEB2021_psscene4band_analytic_sr_udm2/files/PSScene4Band/20210223_123349_104e/analytic_sr_udm2", overwrite = T) 
raster49 <- dir_info(here::here("CASA","Dissertation","DIS_FRANCINA_2021","Slum_detection","data","SantoDominfo_FEB2021_psscene4band_analytic_sr_udm2","files","PSScene4Band","20210223_123349_104e","analytic_sr_udm2")) %>% 
  dplyr::filter(str_detect(path, "20210223_123349_104e_3B_AnalyticMS_SR.tif")) %>% 
  dplyr::select(path) %>% 
  pull() %>% 
  as.character() %>% 
  # Load the raster layers into a stack
  stack()




  






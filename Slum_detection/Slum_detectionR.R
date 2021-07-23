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
install.packages('')

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
test <- st_read("https://github.com/francicabrera/GISassessment/blob/456328fe369c37a3eb910b7f49687379cc047818/Data/geo/ResponseZones/SALUD_DN.shp")
  
  

#https://github.com/francicabrera/GISassessment/blob/456328fe369c37a3eb910b7f49687379cc047818/Data/geo/ResponseZones/SALUD_DN.shp


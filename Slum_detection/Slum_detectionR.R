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
library(e1071)
library(corrr)

################################################################################
# 1. Load and prepare the data
# # a.1) Distrito Nacional's boundary
# #dn_boundary <- st_read(here::here("data","geo","DN_boundary.shp")) # Test
# dn_boundary <- st_read(here("data","geo","DN_boundary.shp")) %>% 
#   # Project to EPSG:32619. This is the projected coordinate system for the Dominican Republic.
#   st_transform(.,32619)
# # Quick look of the feature
# qtm(dn_boundary)

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

# Quick look at the images by the different bands
plot(raster49, main="PlanetScope Distrito Nacional 20210223_123349_104e_3B_AnalyticMS_SR")
plot(raster50, main="PlanetScope Distrito Nacional 20210223_123350_104e_3B_AnalyticMS_SR")

# Merge the images to create a mosaic using the function mosaic()
raster_mosaic <- mosaic(raster49, raster50, fun=mean, tolerance=0.05, filename="raster_mosaic", overwrite=TRUE)

# Quick look at the mosaics
plot(raster_mosaic, main="Mosaic of PlanetScope scenes")

# # Crop the mosaic to the Distrito Nacional's boundary with crop()
# mosaic_DN <- crop(raster_mosaic, dn_boundary, filename="mosaicDN_crop", overwrite=TRUE) %>% 
#   mask(., dn_boundary)

# Crop the mosaic to the Circunscription 3 (C3)'s boundary with crop()
mosaic_C3 <- crop(raster_mosaic, c3_boundary, filename="mosaicC3_crop", overwrite=TRUE) %>% 
  mask(., c3_boundary)

# Quick look at the cropped raster
#plot(mosaic_DN, main="Mosaic cropped to Distrito Nacional's boundary")
plot(mosaic_C3, main="Mosaic cropped to Circunscription 3's boundary")

# # Information on the Distrito Nacional's raster stack
# extent(mosaic_DN) # extent
# ncell(mosaic_DN) # number of cells
# dim(mosaic_DN) # number of rows, columns, layers
# nlayers(mosaic_DN) # number of layers
# res(mosaic_DN) # xres, yres

# Information on the Circunscription 3 (C3)'s raster stack 
extent(mosaic_C3) # extent
ncell(mosaic_C3) # number of cells
dim(mosaic_C3) # number of rows, columns, layers
nlayers(mosaic_C3) # number of layers
res(mosaic_C3) # xres, yres

# Name the Bands based on where they sample the electromagnetic spectrum
# DN's mosaic
# names(mosaic_DN) <- c('blue', 'green', 'red', 'NIR') 
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

# Cloud Assessment 
# PlanetScope provides Usable Data Masks (UDM) to assess the data quality of their images
# See more here: https://developers.planet.com/docs/data/udm-2/
# Let's verify the appearance of pixels covered with clouds or shadow.
# Create a temporary file that will be substituted with the downloaded UDM raster.
temp49_mask <- tempfile(fileext = ".tif")
# Download the UDM raster with the file's id from Google Drive.
dl49_mask <- drive_download(as_id("18r619B2o98q1JHxtkSyPFArvVKqKcpev"),
                       path = temp49_mask, overwrite = TRUE)
# Create the raster stack
raster49_mask <- stack(temp49_mask)
# Check the coordinate system is WGS 84 / UTM zone 19N,
# if not use: %>% projectRaster(., crs=32619)
crs(raster49_mask)
# Plot the raster to visually verify the appearance on pixels with clouds or shadows.
plot(raster49_mask) 

# Find the most frequent values using freq().
# Consider the meaning of 0 as “FALSE” or “NO”, 1 equals “TRUE” or “YES”.
# See more here: https://pages.cms.hu-berlin.de/EOL/gcg_eo/02_data_quality.html
freq(raster49_mask) 

# Repeat the same process with the UDM raster from the second scene.
# Create a temporary file that will be substituted with the downloaded raster.
temp50_mask <- tempfile(fileext = ".tif")
# Download the UDM raster with the file's id from Google Drive.
dl50_mask <- drive_download(as_id("18HnyLNzhswsV2S3b1EokS5OHXKsesZ99"),
                            path = temp50_mask, overwrite = TRUE)
# Create the raster stack
raster50_mask <- stack(temp50_mask)
# Check the coordinate system is WGS 84 / UTM zone 19N,
# if not use: %>% projectRaster(., crs=32619)
crs(raster50_mask)
# Plot the raster to visually verify the appearance on pixels with clouds or shadows.
plot(raster50_mask)
# Find the most frequent values using freq().
freq(raster50_mask) 




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

# Extract image values at training point locations
trainC3.sr <- raster::extract(mosaic_C3, trainC3, sp=T)

# Convert to data.frame 
trainC3.df <- as.data.frame(trainC3.sr)
trainC3.df$classID <- as.factor(trainC3.df$classID) # convert classID into factor

# Create boxplots of reflectance grouped by land cover class
# Melt dataframe containing point id, classID, and 6 spectral bands
spectra.df <- melt(trainC3.df, id.vars='classID', 
                    measure.vars=c('blue', 'green', 'red', 'NIR'))

# Create boxplots of spectral bands per class
ggplot(spectra.df, aes(x=variable, y=value, color=classID)) +
  geom_boxplot() +
  theme_bw()

# Create 2D scatterplot of image data and locations of training points

# Convert image to data.frame and remove missing values
sr.C3.val <- data.frame(getValues(mosaic_C3))
sr.C3.val <- na.omit(sr.C3.val)

# Randomly sub-sample 1000 to speed up visualisation
sr.C3.val <- sr.C3.val[sample(nrow(sr.C3.val), 1000),]  

# Specify which bands to use for the x and y axis of the plot - Red band
xband <- "red"
yband <- "NIR"

# Create plot of band value density and training data
ggplot() + 
  geom_hex(data = sr.C3.val, aes(x = get(xband), y = get(yband)), bins = 100) + 
  geom_point(data = trainC3.df, aes(x = get(xband), y = get(yband), color=classID, shape=classID), 
             size = 2, inherit.aes = FALSE, alpha=1) + 
  scale_fill_gradientn(colours=c("black","white"), na.value=NA) + 
  scale_x_continuous(xband, limits=c(-10, quantile(sr.C3.val[xband], 0.98, na.rm=T))) +
  scale_y_continuous(yband, limits=c(-10, quantile(sr.C3.val[yband], 0.98, na.rm=T))) +
  scale_color_manual(values=c("red", "blue", "green", "purple","black","brown")) +
  theme_bw()

# Specify which bands to use for the x and y axis of the plot - Blue band
xband2 <- "blue"

# Create plot of band value density and training data
ggplot() + 
  geom_hex(data = sr.C3.val, aes(x = get(xband2), y = get(yband)), bins = 100) + 
  geom_point(data = trainC3.df, aes(x = get(xband2), y = get(yband), color=classID, shape=classID), 
             size = 2, inherit.aes = FALSE, alpha=1) + 
  scale_fill_gradientn(colours=c("black","white"), na.value=NA) + 
  scale_x_continuous(xband2, limits=c(-10, quantile(sr.C3.val[xband2], 0.98, na.rm=T))) +
  scale_y_continuous(yband, limits=c(-10, quantile(sr.C3.val[yband], 0.98, na.rm=T))) +
  scale_color_manual(values=c("red", "blue", "green", "purple","black","brown")) +
  theme_bw()

# Specify which bands to use for the x and y axis of the plot - Green band
xband3 <- "green"

# Create plot of band value density and training data
ggplot() + 
  geom_hex(data = sr.C3.val, aes(x = get(xband3), y = get(yband)), bins = 100) + 
  geom_point(data = trainC3.df, aes(x = get(xband3), y = get(yband), color=classID, shape=classID), 
             size = 2, inherit.aes = FALSE, alpha=1) + 
  scale_fill_gradientn(colours=c("black","white"), na.value=NA) + 
  scale_x_continuous(xband3, limits=c(-10, quantile(sr.C3.val[xband3], 0.98, na.rm=T))) +
  scale_y_continuous(yband, limits=c(-10, quantile(sr.C3.val[yband], 0.98, na.rm=T))) +
  scale_color_manual(values=c("red", "blue", "green", "purple","black","brown")) +
  theme_bw()


# Machine learning model: Random Forest
# The code for this section can be found in this source: https://pages.cms.hu-berlin.de/EOL/gcg_eo/05_machine_learning.html
# Training data (dataframe): trainC3.df
# stack: mosaic_C3

# the randomForest() function expects the dependent variable to be of type factor.
# Use as.factor() for conversion of the classID column.
trainC3.df$classID <- as.factor(trainC3.df$classID)
str(trainC3.df) #allows you to see the classes of the variables (all numeric)

# Include only useful predictors in the model.
trainC3_df <- select(trainC3.df, -c("id", "class_name","circ_num","coords.x1","coords.x2")) 
# The RF algorithm cannot deal with NoData (NA) values. Remove NAs from the data.frame.
is.na(trainC3_df) # check for null values
sum(is.na(trainC3_df)) # check how many null values there are
trainC3_df <- na.omit(trainC3_df)
sum(is.na(trainC3_df))

# Train a randomForest() classification model with the data.frame created in the prior step. 
# The code for this section can be found in this source: https://towardsdatascience.com/random-forest-in-r-f66adf80ec9
# 1. Set a portion of the data aside for testing
sample <- sample.split(trainC3_df$classID, SplitRatio = .80)
train <- subset(trainC3_df, sample == TRUE)
test <- subset(trainC3_df, sample == FALSE)
dim(train)
dim(test)

# 2. Model
RF_modelC3 <- randomForest(classID ~ ., data = train, ntree=750, mtry=2)

# 3. Model performance
print(RF_modelC3)
RF_modelC3$err.rate[,1] # OOB estimate of  error rate

# 4. Variable importance
varImpPlot(RF_modelC3, sort=TRUE, main='Variable importance')

# Red band
partialPlot(RF_modelC3, pred.data=train, x.var = 'red', which.class = '1',  plot = TRUE)
partialPlot(RF_modelC3, pred.data=train, x.var = 'red', which.class = '2',  plot = TRUE)
partialPlot(RF_modelC3, pred.data=train, x.var = 'red', which.class = '3',  plot = TRUE)

# 4. Perform a classification of the image stack using the predict() function. 
# The code for this section can be found in this source: https://pages.cms.hu-berlin.de/EOL/gcg_eo/05_machine_learning.html
# Run predict() to store RF predictions
map <- predict(mosaic_C3, RF_modelC3)

map2 <- RF_modelC3 %>% 
  predict(mosaic_C3, test)

# Plot raster
plot(map)
#plot(map2)

# Write classification to disk
writeRaster(map, filename="predicted_map", datatype="INT1S", overwrite=T)

# 5. Calculate class probabilities for each pixel.
# Run predict() to store RF probabilities for class 1-6
RF_modelC3_p <- predict(mosaic_C3, RF_modelC3, type = "prob", index=c(1:6))

# Plot raster
plot(RF_modelC3_p$layer.2)

# Scale probabilities to integer values 0-100 and write to disk
writeRaster(RF_modelC3_pb*100, filename = 'prob_map2', datatype="INT1S", overwrite=T)


# Extract image values at training point locations
predictC3.sr <- raster::extract(map, trainC3, sp=T)

# Convert to data.frame and convert classID into factor
predictC3.df <- as.data.frame(predictC3.sr)



#Automated hyperparameter optimization

# Define accuracy from 5-fold cross-validation as optimization measure
cv <- tune.control(cross = 4) 

# Use tune.randomForest to assess the optimal combination of ntree and mtry
C3rf.tune500 <- tune.randomForest(classID~., data = train, ntree=c(500), mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 54.29%
C3rf.tune750 <- tune.randomForest(classID~., data = train, ntree=c(750), mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 53.46
C3rf.tune1000 <- tune.randomForest(classID~., data = train, ntree=c(1000), mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 54.11
C3rf.tune1250 <- tune.randomForest(classID~., data = train, ntree=c(1250), mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 54.21
C3rf.tune1500 <- tune.randomForest(classID~., data = train, ntree=c(1500), mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 54.42
C3rf.tune1750 <- tune.randomForest(classID~., data = train, ntree=c(1750), mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 53.65
C3rf.tune2000 <- tune.randomForest(classID~., data = train, ntree=c(2000), mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 53.71

# Store the best model in a new object for further use
C3rf.best <- C3rf.tune$best.model

# Is the parametrization and/or different from your previous model?
print(C3rf.best)


# How accurate is our model?
# The code for this section can be found in this source: https://andrewmaclachlan.github.io/CASA0005repo/advanced-r-maup-and-more-regression.html#cross-validation

# We can run a correlation between the data we left out and the predicted data to assess the accuracy.
actuals_preds <- data.frame(cbind(actuals=trainC3.df$classID,
                                  predicteds=predictC3.df$layer))

actuals_preds_correlation <- actuals_preds %>%
  correlate() %>%
  print()

# We can also use min-max accuracy to see how close the actual and predicted values are, using the equation:
# MinMaxAccuracy = mean (min(actuals, pedicteds) / max(actuals, pedicteds))

min_max_accuracy <- mean(apply(actuals_preds, 1, min, na.rm=TRUE) /
                           apply(actuals_preds, 1, max, na.rm=TRUE))  
min_max_accuracy

# Kappa statistics

# Cross validation












# Create a stratified reference sample using sampleStratified()

strat_smpl <- sampleStratified(map, size = 25, sp = TRUE, na.rm = TRUE)

writeOGR()

plot(strat_smpl)


# 











# ggplot(spectra.df, aes(x=RF_modelC3$err.rate, y=RF_modelC3$ntree, color='red')) +
#   geom_line() +
#   theme_bw()
# 
# 





















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

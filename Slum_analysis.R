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
# install.packages('caTools')
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
raster49 <- stack(here("data","SantoDominfo_FEB2021_psscene4band_analytic_sr_udm2","files","PSScene4Band","20210223_123349_104e","analytic_sr_udm2","20210223_123349_104e_3B_AnalyticMS_SR.tif"))
raster50 <- stack(here("data","SantoDominfo_FEB2021_psscene4band_analytic_sr_udm2","files","PSScene4Band","20210223_123350_104e","analytic_sr_udm2","20210223_123350_104e_3B_AnalyticMS_SR.tif"))

# Check the coordinate system
crs(raster49)
crs(raster50)

# Quick look at the images
plot(raster49, main="PlanetScope Distrito Nacional 20210223_123349_104e_3B_AnalyticMS_SR")
plot(raster50, main="PlanetScope Distrito Nacional 20210223_123350_104e_3B_AnalyticMS_SR")


# c) PlanetScope tiles from Distrito Nacional taken on the 14/01/2021 (Option B)
raster7510 <- stack(here("data","orthotics_Jan2021_psorthotile_analytic_5b_udm2","files","PSOrthoTile","4070244_1947510_2021-01-14_2307","analytic_5b_udm2","4070244_1947510_2021-01-14_2307_BGREN_Analytic.tif"))
raster7610 <- stack(here("data","orthotics_Jan2021_psorthotile_analytic_5b_udm2","files/PSOrthoTile","4070244_1947610_2021-01-14_2307","analytic_5b_udm2","4070244_1947610_2021-01-14_2307_BGREN_Analytic.tif"))
raster7611 <- stack(here("data","orthotics_Jan2021_psorthotile_analytic_5b_udm2","files","PSOrthoTile","4070244_1947611_2021-01-14_2307","analytic_5b_udm2","4070244_1947611_2021-01-14_2307_BGREN_Analytic.tif"))

# Check the coordinate system
crs(raster7510)
crs(raster7610)
crs(raster7611)

# Quick look at the images
plot(raster7510, main="PlanetScope Distrito Nacional 4070244_1947510_2021-01-14_2307_BGREN_Analytic")
plot(raster7610, main="PlanetScope Distrito Nacional 4070244_1947610_2021-01-14_2307_BGREN_Analytic")
plot(raster7611, main="PlanetScope Distrito Nacional 4070244_1947611_2021-01-14_2307_BGREN_Analytic")


# 2. Merge the images to create a mosaic using the function mosaic()
raster_mosaicA <- mosaic(raster49, raster50, fun=mean, tolerance=0.05, filename="mosaicA", overwrite=TRUE)

raster_mosaicB <- mosaic(raster7510, raster7610, raster7611, fun=mean, filename="mosaicB", tolerance=0.05, overwrite=TRUE)

# Quick look at the mosaics
plot(raster_mosaicA, main="Mosaic A. PlanetScope scenes")
plot(raster_mosaicB, main="Mosaic B. PlanetScope tiles")

# Alternatively, create the mosaic using the merge() function, it will take the first image as reference for the values. Uncomment:
# raster_mosaicA2 <- merge(raster49, raster50, tolerance=0.05, filename="mosaicA2", overlap=TRUE, overwrite=TRUE)
# raster_mosaicB2 <- mosaic(raster7510, raster7610, raster7611, fun=mean, tolerance=0.05, filename='mosaicB2', overwrite=TRUE)
# Quick look at the mosaics
# plot(raster_mosaicA2, main="Mosaic A2. PlanetScope scenes")
# plot(raster_mosaicB2, main="Mosaic B2. PlanetScope tiles")


# 3. Crop the mosaic to the Distrito Nacional's boundary with crop()
mosaicA_DN <- crop(raster_mosaicA, dn_boundary, filename="mosaicA_crop", overwrite=TRUE) %>% 
  mask(., dn_boundary)

mosaicB_DN <- crop(raster_mosaicB, dn_boundary, filename="mosaicB_crop", overwrite=TRUE) %>% 
  mask(., dn_boundary)

# Quick look at the cropped raster
plot(mosaicA_DN, main="Mosaic (A) cropped to Distrito Nacional's boundary")
plot(mosaicB_DN, main="Mosaic (B) cropped to Distrito Nacional's boundary")

# Information on the raster stack A
crs(mosaicA_DN) # coordinate system
extent(mosaicA_DN) # extent
ncell(mosaicA_DN) # number of cells
dim(mosaicA_DN) # number of rows, columns, layers
nlayers(mosaicA_DN) # number of layers
res(mosaicA_DN) # xres, yres

# Information on the raster stack B
crs(mosaicB_DN) # coordinate system
extent(mosaicB_DN) # extent
ncell(mosaicB_DN) # number of cells
dim(mosaicB_DN) # number of rows, columns, layers
nlayers(mosaicB_DN) # number of layers
res(mosaicB_DN) # xres, yres


# 4. Alternatively, use this mosaic created in QGIS from the PlanetScope scenes
# following the tutorial from this link: http://www.qgistutorials.com/en/docs/3/raster_mosaicing_and_clipping.html
mosaicC_DN <- stack(here("data","Mosaic","mosaic_DNclipped2.tif"))
# Quick plot
plot(mosaicC_DN, main="Distrito Nacional 23 Feb 2021 (mosaic from QGIS)")

# Information on the raster stack 
crs(mosaicC_DN) # coordinate system
extent(mosaicC_DN) # extent
ncell(mosaicC_DN) # number of cells
dim(mosaicC_DN) # number of rows, columns, layers
nlayers(mosaicC_DN) # number of layers
res(mosaicC_DN) # xres, yres

# # 4.1 Mosaic with Additional no data value = 0
# raster_dnNA0 <- raster(here("data","Mosaic","mosaic_DNclipped.tif"))
# # Quick plot
# plot(raster_dnNA0, main="Distrito Nacional 23 Feb 2021 (mosaic from QGIS - Additional no data value = 0)")
# # Check the coordinate system
# crs(raster_dnNA0)

# 5. C3 raster
mosaicC_C3 <- crop(mosaicC_DN, c3_boundary, filename="mosaicC_C1", overwrite=TRUE) %>% 
  mask(., c3_boundary)

# Information on the raster stack 
crs(mosaicC_C3) # coordinate system
extent(mosaicC_C3) # extent
ncell(mosaicC_C3) # number of cells
dim(mosaicC_C3) # number of rows, columns, layers
nlayers(mosaicC_C3) # number of layers
res(mosaicC_C3) # xres, yres


# 6. Name the Bands based on where they sample the electromagentic spectrum
names(mosaicA_DN) <- c('blue', 'green', 'red', 'NIR') 
names(mosaicB_DN) <- c('blue', 'green', 'red', 'red-edge', 'NIR') 
names(mosaicC_DN) <- c('blue', 'green', 'red', 'NIR') 
names(mosaicC_C3) <- c('blue', 'green', 'red', 'NIR') 

# # 7. Extract bands from the mosaic (C)
# r1C <- mosaicC_DN[[1]] # Band 1: BLUE
# # Quick plot
# plot(r1, main="Band 1 - Blue")
# 
# r2C <- mosaicC_DN[[2]] # Band 2: GREEN
# # Quick plot
# plot(r2, main="Band 2 - Green")
# 
# r3C <- mosaicC_DN[[3]] # Band 3: RED
# # Quick plot
# plot(r3, main="Band 3 - Red")
# 
# r4C <- mosaicC_DN[[4]] # Band 4: NIR
# # Quick plot
# plot(r4, main="Band 4 - NIR")

# 8. Plot the data in true colours and false composite.
# true colour composite
rC_rgb <- stack(mosaicC_DN$red, mosaicC_DN$green, mosaicC_DN$blue) %>% 
  plotRGB(.,axes=TRUE, stretch="lin")
# false colour composite
rC_false <- stack(mosaicC_DN$NIR, mosaicC_DN$red, mosaicC_DN$green) %>% 
  plotRGB(.,axes=TRUE, stretch="lin")

# 9. Check the similarity between bands
# set the plot window size (2 by 2)
par(mfrow = c(2,2))
# plot the bands - Mosaic C
plot(mosaicC_DN$blue, main = "Blue")
plot(mosaicC_DN$green, main = "Green")
plot(mosaicC_DN$red, main = "Red")
plot(mosaicC_DN$NIR, main = "NIR")

# set the plot window size (2 by 2)
par(mfrow = c(2,2))
# plot the bands - Mosaic C3
plot(mosaicC_C3$blue, main = "Blue")
plot(mosaicC_C3$green, main = "Green")
plot(mosaicC_C3$red, main = "Red")
plot(mosaicC_C3$NIR, main = "NIR")



# Statistics of these bands
pairs(mosaicC_DN[[1:4]])

# Mosaic A
mosaicA_DN %>%
  as.data.frame(., na.rm=TRUE) %>%
  sample_n(., 100) %>%
  ggpairs(.,axisLabels="none")

# Mosaic B
mosaicB_DN %>%
  as.data.frame(., na.rm=TRUE) %>%
  sample_n(., 100) %>%
  ggpairs(.,axisLabels="none")

# Mosaic C
mosaicC_DN %>%
  as.data.frame(., na.rm=TRUE) %>%
  sample_n(., 100) %>%
  ggpairs(.,axisLabels="none")

# Mosaic C3
mosaicC_C3 %>%
  as.data.frame(., na.rm=TRUE) %>%
  sample_n(., 100) %>%
  ggpairs(.,axisLabels="none")

# Training data - Circumscription 3 (C3)
# The training data was collected digitally using a very high resolution (VHR) imagery as a reference: the Google maps terrain in QGIS.
# following the classification for informal and formal settlements from: https://ieeexplore.ieee.org/document/6236225.
# The location of the informal settlements in the Distrito Nacional was obtained from: http://adn.gob.do/joomlatools-files/docman-files/borrador_plan_est/Borrador%20Plan%20Estrategico%20del%20Distrito%20Nacional%202030%20%20%20V.%2028%20JUL%202020.pdf (page 75)
# The code for this section can be found in this source: https://pages.cms.hu-berlin.de/EOL/gcg_eo/03_training_data.html

# # Load the shapefile with the training data for the model.
# training_dataC3 <- st_read(here("data","Training_data","training_dataC3.shp")) %>% 
#   # Project to EPSG:32619. This is the projected coordinate system for the Dominican Republic.
#   st_transform(.,32619)
# # Quick look of the feature
# qtm(training_dataC3)

# Read training points, the following code assumes that it contains only the class attribute
# in readOGR, dsn specifies the path to the folder containing the file (may not end with /), 
# layer specifies the name of the shapefile without extension (.shp)
trainC3 <- readOGR(dsn='data/Training_data', layer='training_dataC3')

qtm(trainC3) # Plot
crs(trainC3) # Check crs

# Extract image values at training point locations
trainC3.sr <- raster::extract(mosaicC_C3, trainC3, sp=T)

# Convert to data.frame and convert classID into factor
trainC3.df <- as.data.frame(trainC3.sr)
trainC3.df$classID <- as.factor(trainC3.df$classID)

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
sr.C3.val <- data.frame(getValues(mosaicC_C3))
sr.C3.val <- na.omit(sr.C3.val)

# Randomly sub-sample 100,000 to speed up visualisation
sr.C3.val <- sr.C3.val[sample(nrow(sr.C3.val), 100000),]  

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
  scale_color_manual(values=c("red", "blue", "green", "purple")) +
  theme_bw()

# Specify which bands to use for the x and y axis of the plot - Blue band
xband2 <- "blue"
yband2 <- "NIR"

# Create plot of band value density and training data
ggplot() + 
  geom_hex(data = sr.C3.val, aes(x = get(xband2), y = get(yband2)), bins = 100) + 
  geom_point(data = trainC3.df, aes(x = get(xband2), y = get(yband2), color=classID, shape=classID), 
             size = 2, inherit.aes = FALSE, alpha=1) + 
  
  scale_fill_gradientn(colours=c("black","white"), na.value=NA) + 
  scale_x_continuous(xband2, limits=c(-10, quantile(sr.C3.val[xband2], 0.98, na.rm=T))) +
  scale_y_continuous(yband2, limits=c(-10, quantile(sr.C3.val[yband2], 0.98, na.rm=T))) +
  scale_color_manual(values=c("red", "blue", "green", "purple")) +
  theme_bw()

# Specify which bands to use for the x and y axis of the plot - Green band
xband3 <- "green"
yband3 <- "NIR"

# Create plot of band value density and training data
ggplot() + 
  geom_hex(data = sr.C3.val, aes(x = get(xband3), y = get(yband3)), bins = 100) + 
  geom_point(data = trainC3.df, aes(x = get(xband3), y = get(yband3), color=classID, shape=classID), 
             size = 2, inherit.aes = FALSE, alpha=1) + 
  
  scale_fill_gradientn(colours=c("black","white"), na.value=NA) + 
  scale_x_continuous(xband3, limits=c(-10, quantile(sr.C3.val[xband3], 0.98, na.rm=T))) +
  scale_y_continuous(yband3, limits=c(-10, quantile(sr.C3.val[yband3], 0.98, na.rm=T))) +
  scale_color_manual(values=c("red", "blue", "green", "purple")) +
  theme_bw()


# Machine learning model: Random Forest
# The code for this section can be found in this source: https://pages.cms.hu-berlin.de/EOL/gcg_eo/05_machine_learning.html
# Training data (dataframe): trainC3.df
# stack: mosaicC_C3

# the randomForest() function expects the dependent variable to be of type factor.
# Use as.factor() for conversion of the classID column.
trainC3.df$classID <- as.factor(trainC3.df$classID)
str(trainC3.df) #allows you to see the classes of the variables (all numeric)

# Also, include only useful predictors in the model.
trainC3_df <- select(trainC3.df, -c("Id", "class_name","circ_num","coords.x1","coords.x2")) 
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

RF_modelC3 <- randomForest(classID ~ ., data = train, ntree=800)
##### I am not sure on how many decision trees should I specify?

# 3. Model performance
print(RF_modelC3)
RF_modelC3$err.rate[,1] # OOB estimate of  error rate: 48.09%

# 4. Variable importance
varImpPlot(RF_modelC3, sort=TRUE, main='Variable importance')

# Red band
partialPlot(RF_modelC3, pred.data=train, x.var = 'red', which.class = '1',  plot = TRUE)
partialPlot(RF_modelC3, pred.data=train, x.var = 'red', which.class = '2',  plot = TRUE)
partialPlot(RF_modelC3, pred.data=train, x.var = 'red', which.class = '3',  plot = TRUE)

# 4. Perform a classification of the image stack using the predict() function. 
# The code for this section can be found in this source: https://pages.cms.hu-berlin.de/EOL/gcg_eo/05_machine_learning.html
# Run predict() to store RF predictions
map <- predict(mosaicC_C3, RF_modelC3)
map2 <- predict(mosaicC_DN, RF_modelC3)

# Plot raster
plot(map)
plot(map2)

# Write classification to disk
writeRaster(map, filename="predicted_map1", datatype="INT1S", overwrite=T)

# 5. Calculate class probabilities for each pixel.
# Run predict() to store RF probabilities for class 1-3
RF_modelC3_p <- predict(mosaicC_C3, RF_modelC3, type = "prob", index=c(1:3))

# Plot raster
plot(RF_modelC3_p)

# Scale probabilities to integer values 0-100 and write to disk
writeRaster(RF_modelC3_p*100, filename = 'prob_map1', datatype="INT1S", overwrite=T)









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
install.packages('glcm')

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
library(diffeR)
library(caret)
library(corrr)
library(grid)
library(RStoolbox)
library(ggsn)
library(ggspatial)
library(grid)
library(glcm)


################################################################################
# 1. Load and prepare the data
# Provinces
provinces <- st_read(here("data","geo", "PROVCenso2010.shp")) %>% 
  # Project to EPSG:32619. This is the projected coordinate system for the Dominican Republic.
  st_transform(.,32619)
qtm(provinces)

# Distrito Nacional's boundary
dn_boundary <- st_read(here("data","geo","DN_boundary.shp")) %>%
  # Project to EPSG:32619. This is the projected coordinate system for the Dominican Republic.
  st_transform(.,32619)
# Quick look of the feature
qtm(dn_boundary)

# Distrito Nacional's circumscriptions
dn_circ <- st_read(here("data","geo","DN_circumscriptions.shp")) %>%
  # Project to EPSG:32619. This is the projected coordinate system for the Dominican Republic.
  st_transform(.,32619) 
# Quick look of the feature
qtm(dn_circ)

# Circumscription 3
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
# Or, 0: no cloud, 1: cloud; 0: no shadow, 1: shadow.
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


# Merge the two original scenes to create a mosaic using the function mosaic()
raster_mosaic <- mosaic(raster49,
                        raster50, 
                        fun=mean, 
                        tolerance=0.05, 
                        filename="raster_mosaic", 
                        overwrite=TRUE)
# Quick look at the mosaic
plot(raster_mosaic, main="Mosaic of PlanetScope scenes")

# Crop the mosaic to the Distrito Nacional's boundary with crop()
mosaic_DN <- crop(raster_mosaic, 
                  dn_boundary, 
                  filename="mosaicDN_crop", 
                  overwrite=TRUE) %>%
  mask(., dn_boundary)

# Crop the mosaic to the Circumscription 3 (C3)'s boundary with crop()
mosaic_C3 <- crop(raster_mosaic, 
                  c3_boundary, 
                  filename="mosaicC3_crop", 
                  overwrite=TRUE) %>% 
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
# Mosaic
names(raster_mosaic) <- c('blue', 'green', 'red', 'NIR') 
# DN's mosaic
names(mosaic_DN) <- c('blue', 'green', 'red', 'NIR') 
# C3' mosaic
names(mosaic_C3) <- c('blue', 'green', 'red', 'NIR')

# Plot the data in true colours and false composite.
# DN
# true colour composite
rC_rgb <- stack(raster_mosaic$red, raster_mosaic$green, raster_mosaic$blue) %>% 
  plotRGB(.,axes=TRUE, stretch="lin")
# false colour composite
rC_false <- stack(raster_mosaic$NIR, raster_mosaic$red, raster_mosaic$green) %>% 
  plotRGB(.,axes=TRUE, stretch="lin")

# C3
# true colour composite
rC3_rgb <- stack(mosaic_C3$red, mosaic_C3$green, mosaic_C3$blue) %>% 
  plotRGB(.,axes=TRUE, stretch="lin")
# false colour composite
rC3_false <- stack(mosaic_C3$NIR, mosaic_C3$red, mosaic_C3$green) %>% 
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

########

# NDVI
# Let’s make a function called NDVIfun
NDVIfun <- function(NIR, Red) {
  NDVI <- (NIR - Red) / (NIR + Red)
  return(NDVI)
}

# Call the function
source('NDVIfun')

# Calculate NDVI
ndvi <- NDVIfun(mosaic_C3$NIR, mosaic_C3$red)

# plot
ndvi %>%
  plot(.,col = rev(terrain.colors(10)), main = "NDVI")

# Let's look at the histogram for this dataset
ndvi %>%
  hist(., breaks = 40, main = "NDVI Histogram", xlim = c(-.3,.8))

veg <- ndvi %>%
  reclassify(., cbind(-Inf, 0.3, NA))

veg %>%
  plot(.,main = 'Possible Veg cover')



#######
# GLCM
# Source: https://zia207.github.io/geospatial-r-github.io/texture-analysis.html#texture-analysis
# Calculate using default 90 degree shift textures_shift1 <- glcm(raster(L5TSR_1986, layer=1)) plot(textures_shift1)
# Calculate over all directions
# Red band
texturesRed <- glcm(raster(mosaic_C3, layer=1), 
                         window = c(25, 25), 
                         statistics = "variance",
                         shift=list(c(2,-2), c(1,-2), c(1,-1), c(2,-1), c(1,0), 
                                    c(2,0), c(0,1), c(1,1), c(2,1), c(0,2), c(1,2), c(2,2)))
plot(texturesRed)

# Green band
texturesGreen <- glcm(raster(mosaic_C3, layer=2), 
                         window = c(25, 25), 
                         statistics = "variance",
                         shift=list(c(2,-2), c(1,-2), c(1,-1), c(2,-1), c(1,0), 
                                    c(2,0), c(0,1), c(1,1), c(2,1), c(0,2), c(1,2), c(2,2)))
plot(texturesGreen)

# Blue band
texturesBlue <- glcm(raster(mosaic_C3, layer=3), 
                          window = c(25, 25), 
                          statistics = "variance",
                          shift=list(c(2,-2), c(1,-2), c(1,-1), c(2,-1), c(1,0), 
                                     c(2,0), c(0,1), c(1,1), c(2,1), c(0,2), c(1,2), c(2,2)))
plot(texturesBlue)

# NIR band
texturesNIR <- glcm(raster(mosaic_C3, layer=4), 
                          window = c(25, 25), 
                          statistics = "variance",
                          shift=list(c(2,-2), c(1,-2), c(1,-1), c(2,-1), c(1,0), 
                                     c(2,0), c(0,1), c(1,1), c(2,1), c(0,2), c(1,2), c(2,2)))
plot(texturesNIR)

# Principal Component Analysis (PCA) of Texture Bands
# To reduce redundancy of texture bands and to determine the appropriate texture features,
# we will apply PCA to all texture images.
# The function rasterPCA() will calculate the PCA of our raster stack and will return a raster brick with multiple layers of PCA scores
# Source code: https://zia207.github.io/geospatial-r-github.io/texture-analysis.html#texture-analysis

# Stack the glcm layers of all bands
mosaic_C3glcm <- stack(texturesRed, texturesGreen, texturesBlue, texturesNIR)

# Develop a PCA model
mosaic_C3PCA < -scale(mosaic_C3glcm)        # scale the data
mosaic_C3PCA[is.na(mosaic_C3PCA)] <- 0  # define zero  all miising vaalues
mosaic_C3PCAmodel <- rasterPCA(mosaic_C3PCA, nComp=4)
summary(mosaic_C3PCAmodel$model)

# Extract first 2 PCs from the model
# Since, the first two PCs account for more of the variability of these textures (see standard deviation)
# we will extract these 2 components.
# The values from these raster layers will be used as features for our classification.
PC1 <- mosaic_C3PCAmodel$map$PC1
PC2 <- mosaic_C3PCAmodel$map$PC2

# Stack into one raster.
PC_glcm <- stack(PC1, PC2)




# Training data 
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
# Check how many data points are per class
trainC3 %>% group_by(classID) %>% count()

# Extract image values at training point locations
trainC3.sr <- raster::extract(mosaic_C3, trainC3, sp=T)
# GLCM
trainC3.sr2 <- raster::extract(PC_glcm, trainC3, sp=T)

# Convert to data.frame 
trainC3.df <- as.data.frame(trainC3.sr)
# GLCM
trainC3.df2 <- as.data.frame(trainC3.sr2)

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
sample <- sample.split(trainC3.df$classID, SplitRatio = .80)
train <- subset(trainC3.df, sample == TRUE)
test <- subset(trainC3.df, sample == FALSE)

# Check how many data points are per class
train %>% group_by(classID) %>% count()
test %>% group_by(classID) %>% count()

# Check the dimension of the objects
dim(train)
dim(test)

# 2. Model
# Automated hyperparameter optimization
# This part of the process uses the package: "e1071"
# Number of trees (ntree): it is unnecessary to tune in the ntree, instead it is recommended
# to set it to a large number and compare across multiple runs of the model.
# Read more here: https://stats.stackexchange.com/questions/348245/do-we-have-to-tune-the-number-of-trees-in-a-random-forest
# Number of variables (mtry): set the number of k-folds that will run to find the optimal parameter.
# Read more here: https://stats.stackexchange.com/questions/348245/do-we-have-to-tune-the-number-of-trees-in-a-random-forest
# Other sources: https://www.youtube.com/watch?v=v5Bmz2eMd7M

# Define accuracy from 5-fold cross-validation as optimization measure
cv <- tune.control(cross = 5) 

# Use tune.randomForest to assess the optimal combination of ntree and mtry
RF_modelC3.tune500 <- tune.randomForest(classID~., data = train, ntree=500, mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 54.29%
RF_modelC3.tune750 <- tune.randomForest(classID~., data = train, ntree=750, mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 53.46
RF_modelC3.tune1000 <- tune.randomForest(classID~., data = train, ntree=1000, mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 54.11
RF_modelC3.tune1250 <- tune.randomForest(classID~., data = train, ntree=1250, mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 54.21
RF_modelC3.tune1500 <- tune.randomForest(classID~., data = train, ntree=1500, mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 54.42
RF_modelC3.tune1750 <- tune.randomForest(classID~., data = train, ntree=1750, mtry=c(2:10), tunecontrol = cv)
#OOB estimate of  error rate: 53.65

# Store the best model in a new object for further use
RF_modelC3 <- RF_modelC3.tune750$best.model

####
# Model 2
cv2 <- trainControl(method = "cv",
                    number = 10,
                    savePredictions = TRUE)
rfGrid <- expand.grid(mtry = (2:4))
RF_modelC32 <- train(classID~., data = train,
                     method = "rf",
                     trControl = cv2,
                     verbose = TRUE,
                     tuneGrid = rfGrid,
                     importance = TRUE)
RF_modelC32
plot(RF_modelC32)
varImp(RF_modelC32, scale = FALSE)


predicted_class_test <- predict(RF_modelC32, test)
predicted_class_test
test[1]
confusionMatrix(predicted_class_test, test[1])

classified <- predict(mosaic_C3,
                      RF_modelC32)
classified
par(mfrow=c(1,2))
plot(classified)
plotRGB(mosaic_C3,"red","green","blue","NIR",stretch="lin")

###

# 3. Model performance
plot(RF_modelC3)
print(RF_modelC3)
RF_modelC3$err.rate[,1] # OOB estimate of  error rate

# 4. Variable importance
varImpPlot(RF_modelC3, sort=TRUE, main='Variable importance')
varImp(RF_modelC3, scale = FALSE)

# Red band
partialPlot(RF_modelC3, pred.data=train, x.var = 'red', which.class = '1',  plot = TRUE)
partialPlot(RF_modelC3, pred.data=train, x.var = 'red', which.class = '2',  plot = TRUE)
partialPlot(RF_modelC3, pred.data=train, x.var = 'red', which.class = '3',  plot = TRUE)

# 5. Accuracy assessment with test data
predicted_class_test <- predict(RF_modelC3, test)
predicted_class_test
test[,4]
confusionMatrix(predicted_class_test, test)

# 6. Perform a classification of the image stack using the predict() function. 
# The code for this section can be found in this source: https://pages.cms.hu-berlin.de/EOL/gcg_eo/05_machine_learning.html
# Run predict() to store RF predictions
map <- predict(mosaic_C3, RF_modelC3)


# Plot raster
plot(map)
freq(map)

# Write classification to disk
writeRaster(map, filename="predicted_map", datatype="INT1S", overwrite=T)

# 6. Calculate class probabilities for each pixel.
# Run predict() to store RF probabilities for class 1-6
RF_modelC3_p <- predict(mosaic_C3, RF_modelC3, type = "prob", index=c(1:6))

# Plot raster of class: informal settlement Type II
plot(RF_modelC3_p$layer.2)

freq(RF_modelC3_p)

# Scale probabilities to integer values 0-100 and write to disk
writeRaster(RF_modelC3_pb*100, filename = 'prob_map2', datatype="INT1S", overwrite=T)






# Extract image values at training point locations
predictC3.sr <- raster::extract(map, trainC3, sp=T)

# Convert to data.frame and convert classID into factor
predictC3.df <- as.data.frame(predictC3.sr)










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

min_max_accuracy <- mean(apply(actuals_preds$, 1, min, na.rm=TRUE) /
                           apply(actuals_preds$, 1, max, na.rm=TRUE))  
min_max_accuracy

# Kappa statistics


# Cross validation












# Create a stratified reference sample using sampleStratified()

strat_smpl <- sampleStratified(map, size = 50, sp = TRUE, na.rm = TRUE)

writeOGR()

plot(strat_smpl)


# 





# Other Method
# https://www.youtube.com/watch?v=ww8KWgT98Hw
set.seed(42) # allows reproducibility to a process with a random component (random forest)
# Set the train control parameters. Method: cross validation, 5 k-folds, verboseIter = True (to print everything)
trainctrl <- trainControl(method = "cv", number = 5, verboseIter = TRUE)

rf.model2 <- train(classID~., data=train, method= "ranger",
                   tuneLength = 10,
                   preProcess = c("center", "scale"),
                   trControl = trainctrl,
                   metric = "Kappa")





################################################################################
# Plots
# Map 1 - Location of Study Area

# Labels
# Distrito Nacional
# choose a point on the surface of each geometry
DN_points <- sf::st_point_on_surface(dn_boundary)
# retrieve the coordinates
DN_coords <- as.data.frame(sf::st_coordinates(DN_points))
DN_coords$Label <- dn_boundary$TOPONIMIA

# C3
# choose a point on the surface of each geometry
DNc_points <- sf::st_point_on_surface(dn_circ)
# retrieve the coordinates
DNc_coords <- as.data.frame(sf::st_coordinates(DNc_points))
DNc_coords$Label <- dn_circ$TOPONIMIA


# Map of Study Area
Map1_main <- ggplot() +
  ggRGB(raster_mosaic,
        r = 3,
        g = 2,
        b = 1,
      stretch = "lin",
      ggLayer = TRUE) +
  geom_sf(data = c3_boundary,
          fill = alpha("#fed98e",0.4)) +
  theme_minimal() +
  labs(x="", 
       y="") +
  geom_sf(data = dn_circ,
          fill = NA,
          lwd = 0.3,
          colour = "#cc4c02") +
  geom_text(data = DNc_coords,
             aes(X,Y,
                 label= Label),
             colour = "white",
            size = 4) +
  theme(legend.position="bottom") +
  scalebar(dn_circ, dist = 5, dist_unit = "km", location = "bottomright",
           transform = TRUE) +
  # annotation_scale(plot_unit = "km",
  #                  aes(location = "br")) +
  north(data=dn_circ,
        location= "bottomright",
        symbol = 10)

# Inset map
Map1_loc <- ggplot() +
  geom_sf(data = provinces,
          fill = "#cccccc",
          lwd = 0.3,
          colour = "white") +
  geom_sf(data = dn_boundary,
          fill = "#fe9929",
          lwd = 0.2) +
  geom_text(data = DN_coords,
            aes(X,Y,
                label= Label),
            colour = "black",
            size = 4,
            nudge_x = 1,
            vjust = -1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.text=element_blank(), axis.ticks=element_blank()) +
  labs(x="", y="") 

# Combine the main map with the inset map.
Map1_main +
  annotation_custom(ggplotGrob(Map1_loc), 
                    ymin = -1, ymax=1, xmin=1, xmax=1)





  
  Map2 <- ggplot() +
  geom_bin2d(data = CasesHRTA1_pd,
             aes(X, Y),
             binwidth = c(0.005, 0.005)) + 
  geom_sf(data = RoadsA1,
          col="#636363",
          size=0.07) +
  theme_minimal() +
  scale_fill_distiller(palette = "GnBu",
                       direction = 1, 
                       name= "Count of cases with High RTs")+
  theme(legend.position="bottom")+
  labs(x="",
       y="")+
  geom_sf(data = EmerUnitsA1_pd,
          col="#8856a7",
          size=1,)+ 
  geom_text(data = EmerUnitsA1_pd,
            aes(X, Y, label = Ficha),
            colour = "#404040",
            size=3,
            vjust = -1)+
  annotation_scale(bar_cols=c("#f0f0f0","#636363"))+
  north(data=CasesHRTA1_pd,
        location= "bottomright",
        symbol = 3)


  
  
  

  tm_layout(inner.margin=c(0.1,0.04,0.04,0.04),
            legend.outside=TRUE,
            legend.outside.position = "right",
            legend.text.size = 0.5,
            legend.height = 0.5)




























################################################################################
################################################################################
################################################################################




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
# MAP LAYOUT
# # true colour composite
# rC_rgb <- stack(raster_mosaic$red, raster_mosaic$green, raster_mosaic$blue) %>% 
#   plotRGB(.,axes = TRUE, 
#           stretch = "lin", 
#           main = "Distrito Nacional") 
# 
# 
# 
#   tm_shape(dn_boundary) +
#     tm_fill(col = NA) +
#     tm_borders(col = "white",
#                lwd = 0.5) +
#     tm_shape(c3_boundary) +
#     tm_polygons(col = "red",
#                 border.col = "white",
#                 lwd=0.5) +
#     tm_shape(dn_circ) +
#     tm_polygons(col = "gray",
#                 border.col = "white",
#                 lwd=0.5)
#   
#   
# 
#   
# tmap_save(rC_rgb,
#           insets_tm = inset,
#           insets_vp=viewport(0.35, 0.22, width = 0.15, height = 0.15),
#           filename="Map1.png",
#           dpi=600)
#   
#   
# 
# 
# Map1 <-  tm_shape(raster_mosaic) +
#   tm_rgb(
#     r = 1,
#     g = 2,
#     b = 3,
#     alpha = NA,
#     saturation = 1,
#     interpolate = TRUE,
#     max.value = 255)
#   
# 
# 
# 
# 

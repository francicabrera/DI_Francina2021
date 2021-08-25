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
# install.packages(' ')

install.packages("remotes")
remotes::install_github("wilsontom/modelmisc")

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
library(RColorBrewer)
library(cowplot) # to arrange inset maps
library(glcm) # to create GLCM model
library(imager) # to convert RGB raster in grayscale
library(wilsontom/modelmisc) # for kappa accuracy
# To install this package, run this code:
# install.packages("remotes")
# remotes::install_github("wilsontom/modelmisc")

################################################################################
# 1. Load and prepare the data
# Dominican Republic Provinces (for visualisation only)
provinces <- st_read(here("data","geo", "PROVCenso2010.shp")) %>% 
  # Project to EPSG:32619. This is the projected coordinate system for the Dominican Republic.
  st_transform(.,32619)
qtm(provinces)

# Distrito Nacional's boundary (for visualisation only)
dn_boundary <- st_read(here("data","geo","DN_boundary.shp")) %>%
  # Project to EPSG:32619. This is the projected coordinate system for the Dominican Republic.
  st_transform(.,32619)
# Quick look of the feature
qtm(dn_boundary)

# Distrito Nacional's circumscriptions (for visualisation only)
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
freq(raster49_mask$cloud) 

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
freq(raster50_mask$cloud) 

# As both raster show no signal of clouds, the mask wont be used. 
# Therefore, we will remain with the original scenes.

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
# This raster will be used for visualisation only.
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

# Information on the Circumscription 3 (C3)'s raster stack 
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
trainC3.sp <- raster::extract(mosaic_C3, trainC3, sp=T)

# Convert to data.frame 
trainC3.df <- as.data.frame(trainC3.sp)

# Create boxplots of reflectance grouped by land cover class
# Melt dataframe containing point id, classID, and 6 spectral bands
spectra.df <- melt(trainC3.df, id.vars='classID', 
                   measure.vars=c('blue', 'green', 'red', 'NIR'))

# Create boxplots of spectral bands per class
ggplot(spectra.df, aes(x=variable, y=value, color=classID)) +
  geom_boxplot() +
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
# The code for this section can be found in this source: 
# https://stackoverflow.com/questions/20776887/stratified-splitting-the-data/30181396
# https://pages.cms.hu-berlin.de/EOL/gcg_eo/06_accuracy_assessment.html
# createDataPartition does a stratified random split of the data.
# Stratified random sampling helps us to validate the map using a sufficient amount of samples for each for the four classes
# Set a portion of the data aside for testing
sample <- createDataPartition(trainC3_df$classID, p = .80, list = FALSE)
train <- trainC3_df[sample,]
test <- trainC3_df[-sample,]

# Check how many data points are per class
train %>% group_by(classID) %>% count()
test %>% group_by(classID) %>% count()

# Check the dimension of the objects
dim(train)
dim(test)


# Model
# Random Forest. Automated hyperparameter optimisation.
# This part of the process uses the package: "e1071". We will optimise the mtry parameter.
# Number of trees (ntree): it is unnecessary to tune in the ntree, instead it is recommended
# to set it to a large number and compare across multiple runs of the model.
# Read more here: https://stats.stackexchange.com/questions/348245/do-we-have-to-tune-the-number-of-trees-in-a-random-forest
# Number of variables (mtry): set the number of k-folds that will run to find the optimal parameter.
# Read more here: https://stats.stackexchange.com/questions/348245/do-we-have-to-tune-the-number-of-trees-in-a-random-forest
# Other sources: https://www.youtube.com/watch?v=v5Bmz2eMd7M

# Define accuracy from 10-fold cross-validation as optimization measure
cv <- tune.control(cross = 10) 

# Use tune.randomForest to assess the optimal combination of ntree and mtry
RF_modelC3.tune <- tune.randomForest(classID~., 
                                     data = train, 
                                     ntree=750, 
                                     mtry=c(2:10), 
                                     importance = TRUE,
                                     tunecontrol = cv)
#OOB estimate of  error rate: 53.46

# Store the best model in a new object for further use
RF_modelC3 <- RF_modelC3.tune$best.model
print(RF_modelC3)

# Check at which number of trees the model stabilises. 
plot(RF_modelC3) 

# Error rate
RF_modelC3$err.rate[,1] # OOB estimate of  error rate

#Variable importance
# Check which variables provide more information
varImpPlot(RF_modelC3, sort=TRUE, main='Variable importance')
varImp(RF_modelC3, scale = FALSE)

# Model Performance
# Extract the predictions from the model, using the test dataset
predClass <- predict(RF_modelC3, test)
predClass

# Use as.factor() for conversion of the classID column in the test dataset
test_acc <- as.factor(test$classID) 

# Assess the accuracy from the confusion matrix
confusionMatrix(predClass, test_acc)

# Producer and user's accuracy
# Code source: https://blogs.fu-berlin.de/reseda/accuracy-statistics-in-r/
# Build a contingency table of the counts at each combination of factor levels.
# Rename the vectors accordingly.
accmat <- table("pred" = predClass, "ref" = test_acc)
accmat
# User's accuracy
UA <- diag(accmat) / rowSums(accmat) * 100
UA
# Producer's accuracy
PA <- diag(accmat) / colSums(accmat) * 100
PA
# Overall accuracy
OA <- sum(diag(accmat)) / sum(accmat) * 100
OA

# Perform a classification of the image stack using the predict() function. 
# Run predict() to store RF predictions
map <- predict(mosaic_C3, RF_modelC3)

mapttest <- map
mapttest[mapttest == 0] <- NA
trim(mapttest)
crop(mapttest, c3_boundary)

freq(mapttest)
# Plot raster
plot(map)
freq(map)

# Write classification to disk
writeRaster(map, filename="predicted_map", datatype="INT1S", overwrite=T)

# # Calculate class probabilities for each pixel.
# # Run predict() to store RF probabilities for class 1-6
# RF_modelC3_p <- predict(mosaic_C3, RF_modelC3, type = "prob", index=c(1:6))
# 
# # Plot raster of class: informal settlement Type II
# plot(RF_modelC3_p)
# freq(RF_modelC3_p)

# Area Adjusted Accuracy assessment
# Calculated as in Olofsson et al. 2014
# Code source: https://blogs.fu-berlin.de/reseda/area-adjusted-accuracies/
# We will use the previous accuracy matrix: 
confmat <- accmat
# get number of pixels per class and convert in km²
imgVal <- as.factor(getValues(map))
# Remove NA values from the classified dataset 
imgVal <- na.omit(imgVal)
# Extract the number of classes
nclass <- length(unique(train$classID))
# Calculate the total area per class
maparea <- sapply(1:nclass, function(x) sum(imgVal == x))
# Transform area in km2
maparea <- maparea * res(map)[1] ^ 2 / 1000000

# Set confidence interval
# 1.96 for a 95% confidence interval, or 1.645 for a 90% confidence interval.
# See more here: http://reddcr.go.cr/sites/default/files/centro-de-documentacion/olofsson_et_al._2014_-_good_practices_for_estimating_area_and_assessing_accuracy_of_land_change.pdf
conf <- 1.96

# total  map area
A <- sum(maparea)
# proportion of area mapped as class i
W_i <- maparea / A
# number of reference points per class
n_i <- rowSums(confmat) 
# population error matrix
p <- W_i * confmat / n_i
p[is.na(p)] <- 0

# area estimation
p_area <- colSums(p) * A
# area estimation confidence interval
p_area_CI <- conf * A * sqrt(colSums((W_i * p - p ^ 2) / (n_i - 1)))

# overall accuracy
OA <- sum(diag(p))
# producers accuracy
PA <- diag(p) / colSums(p)
# users accuracy
UA <- diag(p) / rowSums(p)

# overall accuracy confidence interval 
OA_CI <- conf * sqrt(sum(W_i ^ 2 * UA * (1 - UA) / (n_i - 1)))
# user accuracy confidence interval 
UA_CI <- conf * sqrt(UA * (1 - UA) / (n_i - 1)) 
# producer accuracy confidence interval 
N_j <- sapply(1:nclass, function(x) sum(maparea / n_i * confmat[ , x]) )
tmp <- sapply(1:nclass, function(x) sum(maparea[-x] ^ 2 * confmat[-x, x] / n_i[-x] * ( 1 - confmat[-x, x] / n_i[-x]) / (n_i[-x] - 1)) )
PA_CI <- conf * sqrt(1 / N_j ^ 2 * (maparea ^ 2 * ( 1 - PA ) ^ 2 * UA * (1 - UA) / (n_i - 1) + PA ^ 2 * tmp))

# gather results
result <- matrix(c(p_area, p_area_CI, PA * 100, PA_CI * 100, UA * 100, UA_CI * 100, c(OA * 100, rep(NA, nclass-1)), c(OA_CI * 100, rep(NA, nclass-1))), nrow = nclass)
result <- round(result, digits = 2) 
rownames(result) <- levels(as.factor(train$classID))
colnames(result) <- c("km²", "km²±", "PA", "PA±", "UA", "UA±", "OA", "OA±")
class(result) <- "table"
result


# Now, we will add texture features to see if it improves the slum detection. 
# GLCM - Gray level co-occurrence matrix
# This code follows the methodology of: https://ieeexplore.ieee.org/document/7447704
# Source: https://zia207.github.io/geospatial-r-github.io/texture-analysis.html#texture-analysis,
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
mosaic_C3glcm2 <- scale(mosaic_C3glcm)        # scale the data
mosaic_C3glcm2[is.na(mosaic_C3glcm2)] <- 0  # define zero  all missing values
mosaic_C3PCA <- rasterPCA(mosaic_C3glcm2, nComp=4)
summary(mosaic_C3PCA$model)

# Extract first 2 PCs from the model
# Since, the first two PCs account for more of the variability of these textures (see standard deviation)
# we will extract these 2 components.
# The values from these raster layers will be used as features for our classification.
PCA1 <- mosaic_C3PCA$map$PC1
PCA2 <- mosaic_C3PCA$map$PC2

# Stack into one raster.
PC_glcm <- stack(PCA1, PCA2)

# Plot
plot(PC_glcm)

# Stack with the original mosaic
mosaic_C3ST <- stack(mosaic_C3, PC_glcm)


# Random Forest model: Spectral + Texture values
# Extract image values at training point locations
trainC3.spGLCM <- raster::extract(PC_glcm, trainC3, sp=T)

# Convert to data.frame 
trainC3.dfGLCM <- as.data.frame(trainC3.spGLCM)

# Merge spectral and glcm values
trainC3.dfST <- merge(trainC3.df,trainC3.dfGLCM, by = "id")

# Rename column classID.x as classID
names(trainC3.dfST)[2] <- 'classID'

# the randomForest() function expects the dependent variable to be of type factor.
# Use as.factor() for conversion of the classID column.
trainC3.dfST$classID <- as.factor(trainC3.dfST$classID) 
str(trainC3.dfST) # allows you to see the classes of the variables (all numeric)

# Include only useful predictors in the model.
trainC3_dfST <- select(trainC3.dfST, -c("id", "class_name.x","circ_num.x","coords.x1.x","coords.x2.x",
                                    "classID.y","class_name.y","circ_num.y","coords.x1.y","coords.x2.y")) 
# The RF algorithm cannot deal with NoData (NA) values. Remove NAs from the data.frame.
is.na(trainC3_dfST) # check for null values
sum(is.na(trainC3_dfST)) # check how many null values there are
trainC3_dfST <- na.omit(trainC3_dfST)
sum(is.na(trainC3_dfST))

# Train a randomForest() classification model with the data.frame created in the prior step. 
# 1. Set a portion of the data aside for testing
sampleST <- createDataPartition(trainC3_dfST$classID, p = .80, list = FALSE)
trainST <- trainC3_dfST[sampleST,]
testST <- trainC3_dfST[-sampleST,]

# Check how many data points are per class
trainST %>% group_by(classID) %>% count()
testST %>% group_by(classID) %>% count()

# Check the dimension of the objects
dim(trainST)
dim(testST)

# RF Model
# Define accuracy from 10-fold cross-validation as optimization measure
cvST <- tune.control(cross = 10) 

# Use tune.randomForest to assess the optimal combination of ntree and mtry
RF_modelC3ST.tune <- tune.randomForest(classID~., 
                                     data = trainST, 
                                     ntree=750, 
                                     mtry=c(2:10), 
                                     importance = TRUE,
                                     tunecontrol = cvST)
#OOB estimate of  error rate: 48.16%

# Store the best model in a new object for further use
RF_modelC3ST <- RF_modelC3ST.tune$best.model
print(RF_modelC3ST)

# Check at which number of trees the model stabilises. 
plot(RF_modelC3ST) 

# Error rate
RF_modelC3ST$err.rate[,1] # OOB estimate of  error rate

# Variable importance
# Check which variables provide more information
varImpPlot(RF_modelC3ST, sort=TRUE, main='Variable importance')
varImp(RF_modelC3ST, scale = FALSE)

# Model Performance - Accuracy assessment
# Extract the predictions from the model, using the test dataset
predClassST <- predict(RF_modelC3ST, testST)
predClassST

# Use as.factor() for conversion of the classID column in the test dataset
test_accST <- as.factor(testST$classID) 

# Assess the overall accuracy and Kappa from the confusion matrix
confusionMatrix(predClassST, test_accST)

# Producer and user's accuracy
# Build a contingency table of the counts at each combination of factor levels.
# Rename the vectors accordingly.
accmat_ST <- table("pred" = predClassST, "ref" = test_accST)
accmat_ST
# User's accuracy
UA_ST <- diag(accmat_ST) / rowSums(accmat_ST) * 100
UA_ST
# Producer's accuracy
PA_ST <- diag(accmat_ST) / colSums(accmat_ST) * 100
PA_ST
# Overall accuracy
OA_ST <- sum(diag(accmat_ST)) / sum(accmat_ST) * 100
OA_ST

# Perform a classification of the image stack using the predict() function. 
# Run predict() to store RF predictions
mapST <- predict(mosaic_C3ST, RF_modelC3ST)

# Plot raster
plot(mapST)
freq(mapST)

# Write classification to disk
writeRaster(mapST, filename="predicted_mapST", datatype="INT1S", overwrite=T)

# # Calculate class probabilities for each pixel.
# # Run predict() to store RF probabilities for class 1-6
# RF_modelC3_pST <- predict(mosaic_C3ST, RF_modelC3ST, type = "prob", index=c(1:6))
# 
# # Plot raster of class: informal settlement Type II
# plot(RF_modelC3_pST)
# freq(RF_modelC3_pST)


# Area Adjusted Accuracy assessment
# We will use the previous accuracy matrix: 
confmat_ST <- accmat_ST
# get number of pixels per class and convert in km²
imgVal_ST <- as.factor(getValues(mapST))
# Remove NA values from the classified dataset 
imgVal_ST <- na.omit(imgVal_ST)
# Extract the number of classes
nclass_ST <- length(unique(trainST$classID))
# Calculate the total area per class
maparea_ST <- sapply(1:nclass_ST, function(x) sum(imgVal_ST == x))
# Transform area in km2
maparea_ST <- maparea_ST * res(mapST)[1] ^ 2 / 1000000

# Set confidence interval
# 1.96 for a 95% confidence interval, or 1.645 for a 90% confidence interval.
# See more here: http://reddcr.go.cr/sites/default/files/centro-de-documentacion/olofsson_et_al._2014_-_good_practices_for_estimating_area_and_assessing_accuracy_of_land_change.pdf
conf_ST <- 1.96

# total  map area
A_ST <- sum(maparea_ST)
# proportion of area mapped as class i
W_i_ST <- maparea_ST / A_ST
# number of reference points per class
n_i_ST <- rowSums(confmat_ST) 
# population error matrix 
p_ST <- W_i_ST * confmat_ST / n_i_ST
p_ST[is.na(p_ST)] <- 0

# area estimation
p_area_ST <- colSums(p_ST) * A_ST
# area estimation confidence interval 
p_area_CI_ST <- conf_ST * A_ST * sqrt(colSums((W_i_ST * p_ST - p_ST ^ 2) / (n_i_ST - 1)))

# overall accuracy 
OA_ST <- sum(diag(p_ST))
# producers accuracy 
PA_ST <- diag(p_ST) / colSums(p_ST)
# users accuracy 
UA_ST <- diag(p_ST) / rowSums(p_ST)

# overall accuracy confidence interval 
OA_CI_ST <- conf_ST * sqrt(sum(W_i_ST ^ 2 * UA_ST * (1 - UA_ST) / (n_i_ST - 1)))
# user accuracy confidence interval 
UA_CI_ST <- conf_ST * sqrt(UA_ST * (1 - UA_ST) / (n_i_ST - 1)) 
# producer accuracy confidence interval 
N_j_ST <- sapply(1:nclass_ST, function(x) sum(maparea_ST / n_i_ST * confmat_ST[ , x]) )
tmp_ST <- sapply(1:nclass_ST, function(x) sum(maparea_ST[-x] ^ 2 * confmat_ST[-x, x] / n_i_ST[-x] * ( 1 - confmat_ST[-x, x] / n_i_ST[-x]) / (n_i_ST[-x] - 1)) )
PA_CI_ST<- conf_ST * sqrt(1 / N_j_ST ^ 2 * (maparea_ST ^ 2 * ( 1 - PA_ST ) ^ 2 * UA_ST * (1 - UA_ST) / (n_i - 1) + PA_ST ^ 2 * tmp_ST))

# gather results
result_ST <- matrix(c(p_area_ST, p_area_CI_ST, PA_ST * 100, PA_CI_ST * 100, UA_ST * 100, UA_CI_ST * 100, c(OA_ST * 100, rep(NA, nclass_ST-1)), c(OA_CI_ST * 100, rep(NA, nclass_ST-1))), nrow = nclass_ST)
result_ST <- round(result_ST, digits = 2) 
rownames(result_ST) <- levels(as.factor(trainST$classID))
colnames(result_ST) <- c("km²", "km²±", "PA", "PA±", "UA", "UA±", "OA", "OA±")
class(result_ST) <- "table"
result_ST




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

# Get XY coordinates

c3_boundary <- c3_boundary %>% 
  cbind(st_coordinates(.)) # get X and Y coordinates

# Set the map's properties
# Map of Study Area
Map1_main <- ggplot() +
  ggRGB(raster_mosaic,
        r = 3,
        g = 2,
        b = 1,
      stretch = "lin",
      ggLayer = TRUE) +
  geom_sf(data = c3_boundary,
          aes(fill = TOPONIMIA)) +
  scale_fill_manual(name = 'Study Area', values = c(alpha('#d8b365',0.4))) +
  theme_minimal() +
  labs(x="", 
       y="") +
  geom_sf(data = dn_circ,
          fill = NA,
          lwd = 0.3,
          aes(color = TOPONIMIA)) +
  scale_color_manual(name = 'Circumscriptions', values = c("#cc4c02","#cc4c02","#cc4c02")) +
  geom_text(data = DNc_coords,
             aes(X,Y,
                 label= Label),
             colour = "white",
            size = 5) +
  theme(legend.position= c(1.25, 0.42)) +
  ggsn::scalebar(data = dn_circ,
                 dist = 2, 
                 dist_unit = "km",
                 transform = FALSE,
                 height = 0.01,
                 border.size = 0.5) +
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
            size = 2.5,
            nudge_x = 1,
            vjust = -1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.text=element_blank(), axis.ticks=element_blank()) +
  labs(x="", y="") 


# Combine the main map with the inset map.
Map1 <- ggdraw() +
  draw_plot(Map1_main) +
  draw_plot(Map1_loc, x = 0.714, y = 0.005, width = 0.3, height = 0.3)

Map1

# This will save the map in the same path as the project.
ggsave("Map1.png",
       plot=Map1,
       dpi = 600)



# Map 2 - Training data

# Load zoom polygon
zoomMap2 <- st_read(here("data","geo","zoom_Map2.shp")) %>% 
  # Project to EPSG:32619. This is the projected coordinate system for the Dominican Republic.
  st_transform(.,32619)

# # Get X and Y coordinates
# trainC3_plot <- trainC3 %>% 
#   cbind(st_coordinates(.)) 


# Set the map's properties
Map2 <- ggplot() +
  ggRGB(mosaic_C3,
        r = 3,
        g = 2,
        b = 1,
        stretch = "lin",
        ggLayer = TRUE) +
  geom_sf(data = trainC3,
          size = 0.4,
           aes(color = factor(classID))) +
  scale_color_manual(name = 'Classes',
                     values = c('#8c510a', '#d8b365', '#f6e8c3','#c7eae5','#5ab4ac','#01665e'),
                     labels = c('Informal Type I','Informal Type II', 'Informal Type III',
                                'Formal', 'No-settlement', 'Roads and asphalt')) +
  geom_sf(data = zoomMap2,
         color = 'red') +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key = element_rect(fill = NA),
        legend.position= c(1.2, 0.62),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.background = element_blank()) +
  labs(x="", 
       y="") +
  ggsn::scalebar(data = trainC3,
                 dist = 0.5, 
                 dist_unit = "km",
                 transform = FALSE,
                 height = 0.01,
                 border.size = 0.5,
                 location = "bottomleft") +
  north(data=trainC3,
        location= "topright",
        symbol = 10) 

Map2

# Save the map
ggsave("Map2.png",
       plot = Map2,
       dpi = 600)



# Map 3 - Classification maps

# Map classification #1
# convert raster in data frame
mapVal <-  rasterToPoints(map)    
mapdf <- data.frame(mapVal)
colnames(mapdf) <- c('Longitude', 'Latitude', 'classID')


# Set the map's properties
Map3_1 <- ggplot() +
  geom_raster(data = mapdf,
              aes(y=Latitude,
                  x=Longitude,
                  fill = factor(classID))) +
  coord_equal() +
  scale_fill_manual(name = 'Classes',
                     values = c('#8c510a', '#d8b365', '#f6e8c3','#c7eae5','#5ab4ac','#01665e'),
                     labels = c('Informal Type I','Informal Type II', 'Informal Type III',
                                'Formal', 'No-settlement', 'Roads and asphalt')) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  labs(x="", 
       y="") 

# Map classification #2
# convert raster in data frame 
mapValST <-  rasterToPoints(mapST)    
mapdfST <- data.frame(mapValST)
colnames(mapdfST) <- c('Longitude', 'Latitude', 'classID')

# Set the map's properties
Map3_2 <- ggplot() +
  geom_raster(data = mapdfST,
              aes(y=Latitude,
                  x=Longitude,
                  fill = factor(classID))) +
  coord_equal() +
  scale_fill_manual(name = 'Classes',
                    values = c('#8c510a', '#d8b365', '#f6e8c3','#c7eae5','#5ab4ac','#01665e'),
                    labels = c('Informal Type I','Informal Type II', 'Informal Type III',
                               'Formal', 'No-settlement', 'Roads and asphalt')) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key = element_rect(fill = NA),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank()) +
  labs(x="", 
       y="") +
  ggsn::scalebar(data = trainC3,
                 dist = 0.5, 
                 dist_unit = "km",
                 transform = FALSE,
                 height = 0.01,
                 border.size = 0.5,
                 location = "bottomright") +
  north(data=trainC3,
        location= "topright",
        symbol = 10) 

# Combine the main map with the inset map.
Map3 <- ggdraw() +
  draw_plot(Map3_1, x = -0.27, y = 0) +
  draw_plot(Map3_2, x = 0.21, y = - 0.009)

Map3

# Save the map
ggsave("Map3.png",
       plot=Map3,
       dpi = 600)




# Map 4 - Slum maps (by type)

# Map slums #1
# Subset dataframe 
mapdf_slums <- filter(mapdf, classID == '1' | classID == '2' | classID == '3')
  
# Set the map's properties
Map4_1 <- ggplot() +
  ggRGB(mosaic_C3,
        r = 3,
        g = 2,
        b = 1,
        stretch = "lin",
        ggLayer = TRUE) +
  geom_tile(data = mapdf_slums,
              aes(y=Latitude,
                  x=Longitude,
                  fill = factor(classID))) +
  coord_equal() +
  scale_fill_manual(name = 'Classes',
                    values = c('#8c510a', '#d8b365', '#f6e8c3'),
                    labels = c('Informal Type I','Informal Type II', 'Informal Type III')) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  labs(x="", 
       y="") 

# Map slums #2
# Subset dataframe 
mapdfST_slums <- filter(mapdfST, classID == '1' | classID == '2' | classID == '3')

# Set the map's properties
Map4_2 <- ggplot() +
  ggRGB(mosaic_C3,
        r = 3,
        g = 2,
        b = 1,
        stretch = "lin",
        ggLayer = TRUE) +
  geom_tile(data = mapdfST_slums,
            aes(y=Latitude,
                x=Longitude,
                fill = factor(classID))) +
  coord_equal() +
  scale_fill_manual(name = 'Classes',
                    values = c('#8c510a', '#d8b365', '#f6e8c3'),
                    labels = c('Informal Type I','Informal Type II', 'Informal Type III')) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key = element_rect(fill = NA),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank()) +
  labs(x="", 
       y="") +
  ggsn::scalebar(data = trainC3,
                 dist = 0.5, 
                 dist_unit = "km",
                 transform = FALSE,
                 height = 0.01,
                 border.size = 0.5,
                 location = "bottomright") +
  north(data=trainC3,
        location= "topright",
        symbol = 10) 

# Combine the main map with the inset map.
Map4 <- ggdraw() +
  draw_plot(Map4_1, x = -0.27, y = 0) +
  draw_plot(Map4_2, x = 0.21, y = - 0.009)

Map4

# Save the map
ggsave("Map4.png",
       plot=Map4,
       dpi = 600)




# Map 5 - Slums (single color)

# Map slums #1
# Set the map's properties
Map5_1 <- ggplot() +
  ggRGB(mosaic_C3,
        r = 3,
        g = 2,
        b = 1,
        stretch = "lin",
        ggLayer = TRUE) +
  geom_tile(data = mapdf_slums,
            aes(y=Latitude,
                x=Longitude,
                fill = factor(classID))) +
  coord_equal() +
  scale_fill_manual(name = 'Classes',
                    values = c('red', 'red', 'red'),
                    labels = c('Informal Type I','Informal Type II', 'Informal Type III')) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  labs(x="", 
       y="") 

# Map slums #2
# Set the map's properties
Map5_2 <- ggplot() +
  ggRGB(mosaic_C3,
        r = 3,
        g = 2,
        b = 1,
        stretch = "lin",
        ggLayer = TRUE) +
  geom_tile(data = mapdf_slums,
            aes(y=Latitude,
                x=Longitude,
                fill = factor(classID))) +
  coord_equal() +
  scale_fill_manual(name = 'Classes',
                    values = c('red', 'red', 'red'),
                    labels = c('Informal Type I','Informal Type II', 'Informal Type III')) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key = element_rect(fill = NA),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank()) +
  labs(x="", 
       y="") +
  ggsn::scalebar(data = trainC3,
                 dist = 0.5, 
                 dist_unit = "km",
                 transform = FALSE,
                 height = 0.01,
                 border.size = 0.5,
                 location = "bottomright") +
  north(data=trainC3,
        location= "topright",
        symbol = 10) 


# Combine the main map with the inset map.
Map5 <- ggdraw() +
  draw_plot(Map5_1, x = -0.27, y = 0) +
  draw_plot(Map5_2, x = 0.21, y = - 0.009)

Map5

# Save the map
ggsave("Map5.png",
       plot=Map5,
       dpi = 600)



# Graphic 1 - Scatterplots Predicted vs Observed

# Join predicted and observe into one data frame
comparisonsST <- merge(xx, testST, by="id")

# Statistics of these bands
# Mosaic C3
mosaic_C3 %>%
  as.data.frame(., na.rm=TRUE) %>%
  sample_n(., 100) %>%
  ggpairs(.,axisLabels="none")


predClassST, test_accST

xx <- as.data.frame(predClassST)

# We can run a correlation between the data we left out and the predicted data to assess the accuracy.
actuals_preds <- data.frame(cbind(actuals = testST$classID,
                                  predicteds = xx$predClassST))

actuals_preds_correlation <- actuals_preds %>%
  correlate() %>%
  print()




plot(xx$predClassST, testST$classID)


plot(computeddata$Temp, computeddata$NDBI)

ggplot(actuals_preds, 
       aes(x=actuals, y=predicteds)) +
  geom_hex(bins=100, na.rm=TRUE) +
  labs(fill = "Count per bin")+
  geom_smooth(method='lm', se=FALSE, size=0.6)+
  theme_bw()





















# Map  - Classification maps
# map and mapST
plot(map)
plot(mapST)

# convert raster in data frame
mapVal <-  rasterToPoints(map)    # as.data.frame(map$layer, xy = T)
mapdf <- data.frame(mapVal)
colnames(mapdf) <- c('Longitude', 'Latitude', 'classID')

ggplot(data=mapdf, aes(y=Latitude, x=Longitude)) +
  geom_raster(aes(fill=factor(classID)))+
  theme_minimal() +
  labs(x="", 
       y="")



map <- ratify(map)
rat <- levels(map)[[1]]
# Create new columns for land cover classes
rat$landcover <- c('Informal Type I','Informal Type II', 'Informal Type III',
                   'Formal', 'No-settlement', 'Roads and asphalt')
rat$class <- c('1', '2', '3', '4', '5', '6')
# set the map levels
levels(map) <- rat
# remove na values





# Do the same for the second map
mapST <- ratify(mapST)
levels(mapST) <- rat


Map <- 

gplot(map) +
  geom_tile(aes(fill = factor(value, labels = c('Informal Type I','Informal Type II', 'Informal Type III',
                                                'Formal', 'No-settlement', 'Roads and asphalt')))) +
  scale_fill_manual(values = c('#8c510a', '#d8b365', '#f6e8c3','#c7eae5','#5ab4ac','#01665e'), name = 'Classes') +
  coord_equal()  +
  theme_minimal() +
  labs(x="", 
       y="") 


+
  north(data=dn_circ,
        location= "bottomright",
        symbol = 10)





levelplot(map, col.regions=c('#8c510a', '#d8b365', '#f6e8c3','#c7eae5','#5ab4ac','#01665e')) 

levelplot(map, col.regions=c('#8c510a', '#d8b365', '#f6e8c3','#c7eae5','#5ab4ac','#01665e'))


#8c510a
#d8b365
#f6e8c3
#c7eae5
#5ab4ac
#01665e

Map2a <- ggplot() +
  geom_raster(map, aes(x = x, y = y))
  



# en theme, para quitar
panel.grid.major = element_blank()





# Map  - Informal settlements

slums <- map %>%
  reclassify(., cbind(-Inf, 0.3, NA))

slums %>%
  plot(.,main = 'Possible Veg cover')

























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
# # # Plot the data in true colours and false composite.
#   # DN
#   # true colour composite
#   rC_rgb <- stack(raster_mosaic$red, raster_mosaic$green, raster_mosaic$blue) %>% 
#     plotRGB(.,axes=TRUE, stretch="lin")
#   # false colour composite
#   rC_false <- stack(raster_mosaic$NIR, raster_mosaic$red, raster_mosaic$green) %>% 
#     plotRGB(.,axes=TRUE, stretch="lin")
# # 
# # 
#   # Create 2D scatterplot of image data and locations of training points
#   
#   # Convert image to data.frame and remove missing values
#   sr.C3.val <- data.frame(getValues(mosaic_C3))
#   sr.C3.val <- na.omit(sr.C3.val)
#   
#   # Randomly sub-sample 1000 to speed up visualisation
#   sr.C3.val <- sr.C3.val[sample(nrow(sr.C3.val), 1000),]  
#   
#   # Specify which bands to use for the x and y axis of the plot - Red band
#   xband <- "red"
#   yband <- "NIR"
#   
#   # Create plot of band value density and training data
#   ggplot() + 
#     geom_hex(data = sr.C3.val, aes(x = get(xband), y = get(yband)), bins = 100) + 
#     geom_point(data = trainC3.df, aes(x = get(xband), y = get(yband), color=classID, shape=classID), 
#                size = 2, inherit.aes = FALSE, alpha=1) + 
#     scale_fill_gradientn(colours=c("black","white"), na.value=NA) + 
#     scale_x_continuous(xband, limits=c(-10, quantile(sr.C3.val[xband], 0.98, na.rm=T))) +
#     scale_y_continuous(yband, limits=c(-10, quantile(sr.C3.val[yband], 0.98, na.rm=T))) +
#     scale_color_manual(values=c("red", "blue", "green", "purple","black","brown")) +
#     theme_bw()
#   
#   # Specify which bands to use for the x and y axis of the plot - Blue band
#   xband2 <- "blue"
#   
#   # Create plot of band value density and training data
#   ggplot() + 
#     geom_hex(data = sr.C3.val, aes(x = get(xband2), y = get(yband)), bins = 100) + 
#     geom_point(data = trainC3.df, aes(x = get(xband2), y = get(yband), color=classID, shape=classID), 
#                size = 2, inherit.aes = FALSE, alpha=1) + 
#     scale_fill_gradientn(colours=c("black","white"), na.value=NA) + 
#     scale_x_continuous(xband2, limits=c(-10, quantile(sr.C3.val[xband2], 0.98, na.rm=T))) +
#     scale_y_continuous(yband, limits=c(-10, quantile(sr.C3.val[yband], 0.98, na.rm=T))) +
#     scale_color_manual(values=c("red", "blue", "green", "purple","black","brown")) +
#     theme_bw()
#   
#   # Specify which bands to use for the x and y axis of the plot - Green band
#   xband3 <- "green"
#   
#   # Create plot of band value density and training data
#   ggplot() + 
#     geom_hex(data = sr.C3.val, aes(x = get(xband3), y = get(yband)), bins = 100) + 
#     geom_point(data = trainC3.df, aes(x = get(xband3), y = get(yband), color=classID, shape=classID), 
#                size = 2, inherit.aes = FALSE, alpha=1) + 
#     scale_fill_gradientn(colours=c("black","white"), na.value=NA) + 
#     scale_x_continuous(xband3, limits=c(-10, quantile(sr.C3.val[xband3], 0.98, na.rm=T))) +
#     scale_y_continuous(yband, limits=c(-10, quantile(sr.C3.val[yband], 0.98, na.rm=T))) +
#     scale_color_manual(values=c("red", "blue", "green", "purple","black","brown")) +
#     theme_bw()
  
  
  # Include only useful predictors in the model.
  trainC3_df <- select(trainC3.df, -c("id", "class_name.x","circ_num.x","coords.x1.x","coords.x2.x",
                                      "classID.y","class_name.y","circ_num.y","coords.x1.y","coords.x2.y")) 
  
  
  RF_modelC3.tune1000 <- tune.randomForest(classID~., data = train, ntree=1000, mtry=c(2:10), tunecontrol = cv)
  #OOB estimate of  error rate: 54.11
  RF_modelC3.tune1250 <- tune.randomForest(classID~., data = train, ntree=1250, mtry=c(2:10), tunecontrol = cv)
  #OOB estimate of  error rate: 54.21
  RF_modelC3.tune1500 <- tune.randomForest(classID~., data = train, ntree=1500, mtry=c(2:10), tunecontrol = cv)
  #OOB estimate of  error rate: 54.42
  RF_modelC3.tune1750 <- tune.randomForest(classID~., data = train, ntree=1750, mtry=c(2:10), tunecontrol = cv)
  #OOB estimate of  error rate: 53.65
#   
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
  
  
  # Train a randomForest() classification model with the data.frame created in the prior step. 
  # The code for this section can be found in this source: https://towardsdatascience.com/random-forest-in-r-f66adf80ec9
  # # 1. Set a portion of the data aside for testing
  # sampleST <- sample.split(trainC3_dfST$classID, SplitRatio = .80)
  # trainST <- subset(trainC3_dfST, sample == TRUE)
  # testST <- subset(trainC3_dfST, sample == FALSE)
  # 
  # # Check how many data points are per class
  # trainST %>% group_by(classID) %>% count()
  # testST %>% group_by(classID) %>% count()
  # 
  # # Check the dimension of the objects
  # dim(trainST)
  # dim(testST)
  
  # 
  # # Train a randomForest() classification model with the data.frame created in the prior step. 
  # # The code for this section can be found in this source: https://towardsdatascience.com/random-forest-in-r-f66adf80ec9
  # # 1. Set a portion of the data aside for testing
  # sample <- sample.split(trainC3_df$classID, SplitRatio = .80)
  # train <- subset(trainC3_df, sample == TRUE)
  # test <- subset(trainC3_df, sample == FALSE)
  # 
  # # Check how many data points are per class
  # train %>% group_by(classID) %>% count()
  # test %>% group_by(classID) %>% count()
  # 
  # # Check the dimension of the objects
  # dim(train)
  # dim(test)
  
  # 
  # # How accurate is our model?
  # # The code for this section can be found in this source: https://andrewmaclachlan.github.io/CASA0005repo/advanced-r-maup-and-more-regression.html#cross-validation
  # 
  # # We can run a correlation between the data we left out and the predicted data to assess the accuracy.
  # actuals_preds <- data.frame(cbind(actuals=trainC3.df$classID,
  #                                   predicteds=predictC3.df$layer))
  # 
  actuals_preds_correlation <- actuals_preds %>%
    correlate() %>%
    print()
  # 
  # # We can also use min-max accuracy to see how close the actual and predicted values are, using the equation:
  # # MinMaxAccuracy = mean (min(actuals, pedicteds) / max(actuals, pedicteds))
  # 
  # min_max_accuracy <- mean(apply(actuals_preds$, 1, min, na.rm=TRUE) /
  #                            apply(actuals_preds$, 1, max, na.rm=TRUE))  
  # min_max_accuracy
  # 
  # # Kappa statistics
  # 
  # 
  # # Cross validation
  
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
  
  ######
  
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
  
  
  fill = alpha('#d8b365',0.4),
  
  ggplot() +
    geom_sf(data = poly, aes(fill = "A")) +
    geom_sf(data = point, aes(colour = "B"), show.legend = "point") +
    geom_sf(data = line, aes(colour = "C"), show.legend = "line") +
    scale_fill_manual(values = c("A" = "yellow")) +
    scale_colour_manual(values = c("B" = "pink", "C" = "purple")) +
    theme_minimal()
  
  scale_color_manual(values = c("Women" = '#ff00ff','Men' = '#3399ff')) + 
    scale_shape_manual(values = c('Women' = 17, 'Men' = 16))
  
  
  
  Map1 <- Map1_main + annotation_custom(ggplotGrob(Map1_loc), 
                                        ymin = -1, ymax=1, xmin=1, xmax=1)
  plot(Map1)
  
library(terra)
library(randomForest)
library(caret)
library(sf)
library(plyr) 
library(dplyr)
library(tidyverse)
library(caret) 

# ---------------------------------------------------------------
# Load and align rasters
# ---------------------------------------------------------------
matching_tifs <- c("ETH_GlobalCanopyHeight_2020_10m_v1", "WSCI", 
                   "HLSLandsat", "HLSSentinel", "SentinelSAR","GEDI_biomass_filled")

rasterList <- lapply(matching_tifs, function(x)
  rast(paste0("/Users/abarenb1/Documents/PhD/Chapter3/InputsRandomForest/", x, ".tif")))

template <- rasterList[[1]]

rasterList_aligned <- lapply(seq_along(rasterList), function(i) {
  r <- rasterList[[i]]
  if (!compareGeom(r, template, stopOnError = FALSE)) {
    print(paste("Aligning", matching_tifs[i]))
    method <- ifelse(grepl("classified", matching_tifs[i]), "near", "bilinear")
    project(r, template, method = method)
  } else {
    r
  }
})

# ---------------------------------------------------------------
# Clean band names BEFORE stacking (fixes formula error too)
# ---------------------------------------------------------------
clean_names <- function(r, prefix = "") {
  n <- names(r)
  n <- gsub("[^a-zA-Z0-9_]", "_", n)   # Replace spaces/special chars with _
  n <- gsub("_+", "_", n)              # Collapse multiple underscores
  n <- gsub("^_|_$", "", n)            # Remove leading/trailing underscores
  if (prefix != "") n <- paste0(prefix, "_", n)
  names(r) <- n
  r
}

rasterList_aligned[[1]] <- clean_names(rasterList_aligned[[1]], "ETH")
rasterList_aligned[[2]] <- clean_names(rasterList_aligned[[2]], "WSCI")
rasterList_aligned[[3]] <- clean_names(rasterList_aligned[[3]], "Landsat")
rasterList_aligned[[4]] <- clean_names(rasterList_aligned[[4]], "Sentinel")
rasterList_aligned[[5]] <- clean_names(rasterList_aligned[[5]], "SAR")
rasterList_aligned[[6]] <- clean_names(rasterList_aligned[[6]], "Biomass")

# ---------------------------------------------------------------
# Stack rasters
# ---------------------------------------------------------------
rasterStack <- do.call(c, rasterList_aligned)
#rasterStack<-mask(rasterStack,rasterList_aligned[[6]])
bands <- names(rasterStack)
cat("Band names:\n"); print(bands)
cat("Total layers:", nlyr(rasterStack), "\n")

# ---------------------------------------------------------------
# Load GEDI data
# ---------------------------------------------------------------
L2A <- st_read("/Users/abarenb1/Documents/PhD/Chapter3/Variables/Bwindi_L2A.shp")
L2B <- st_read("/Users/abarenb1/Documents/PhD/Chapter3/Variables/Bwindi_L2B.shp")

L2B <- L2B %>% rename(fhd = fhd_normal)

fhd   <- L2B[, "fhd"]
cover <- L2B[, "cover"]
pai   <- L2B[, "pai"]
rh98  <- L2A[, "rh98"]
rh100 <- L2A[, "rh100"]

# Reproject to match raster CRS
fhd <- st_transform(fhd, crs = st_crs(crs(rasterStack)))

# ---------------------------------------------------------------
# 🔧 KEY FIX: Check NA coverage per layer at point locations
# ---------------------------------------------------------------
fhd_vect <- vect(fhd)

# Extract WITHOUT removing NAs first
samples_raw <- terra::extract(rasterStack, fhd_vect, df = TRUE)
samples_raw$fhd <- fhd$fhd

cat("\nNA count per band at GEDI point locations:\n")
na_counts <- sapply(samples_raw[, bands], function(x) sum(is.na(x)))
print(na_counts)

cat("\nTotal points before NA removal:", nrow(samples_raw), "\n")
cat("Complete cases (no NA in ANY band):", sum(complete.cases(samples_raw)), "\n")

# ---------------------------------------------------------------
# 🔧 OPTION A: Keep only bands with low NA counts
#    Drop bands where >50% of points are NA
# ---------------------------------------------------------------
na_frac <- na_counts / nrow(samples_raw)
bad_bands <- names(na_frac[na_frac > 0.5])

if (length(bad_bands) > 0) {
  cat("\n⚠️  Dropping high-NA bands (>50% NA at point locations):\n")
  print(bad_bands)
  bands_use <- setdiff(bands, bad_bands)
} else {
  bands_use <- bands
}

cat("\nBands used for Random Forest:", length(bands_use), "\n")

# ---------------------------------------------------------------
# 🔧 OPTION B: Use na.action in randomForest instead of na.omit
#    This keeps more samples by imputing or ignoring NAs per tree
# ---------------------------------------------------------------
samples <- samples_raw[, c("ID", bands_use, "fhd")]
samples <- na.omit(samples)   # Only remove rows with NA in KEPT bands

cat("\nSamples after NA removal:", nrow(samples), "\n")

# ---------------------------------------------------------------
# Train/test split
# ---------------------------------------------------------------
set.seed(123)
samples$random <- runif(nrow(samples))

split <- 0.8
training <- samples[samples$random < split, ]
testing  <- samples[samples$random >= split, ]

cat("Samples n =",  nrow(samples),  "\n")
cat("Training n =", nrow(training), "\n")
cat("Testing n =",  nrow(testing),  "\n")

#*************************************************************
# Train Random Forest Model
#*************************************************************
#cover$target_class <- cut(pai$pai, 
#                         breaks = quantile(pai$pai, probs = seq(0, 1, 0.10)),
#                         labels = c("1", "2", "3", "4","5","6","7","8","9","10"),
#                         include.lowest = TRUE)
#
#pai$target_class <- cut(pai$pai, 
#                        breaks = quantile(pai$pai, probs = c(0, 0.1, 0.9, 1.0)),
#                        labels = c("Low", "Medium", "High"),
#                        include.lowest = TRUE)
#
#class_counts <- table(pai$target_class)
#total_samples <- sum(class_counts)
#num_classes <- length(class_counts)
#class_weights <- total_samples / (num_classes * class_counts)
#names(class_weights) <- levels(pai$target_class)
#
## Prepare formula (all bands predicting pai)
formula <- as.formula(paste("fhd ~", paste(bands, collapse = " + ")))

# Train Random Forest
# ntree = 250 trees, mtry = 20 variables tried at each split
classifier <- randomForest(
  formula = formula,
  data = training,
  ntree = 500,
  mtry = 10,
  importance = TRUE,
  na.action = na.omit
)

print(classifier)

#############################################################
#k <- 5  # Number of folds
#folds <- createFolds(training$pai, k = k, list = TRUE, returnTrain = FALSE)
#
## Initialize vectors to store results
#rmse_values <- numeric(k)
#mae_values <- numeric(k)
#r2_values <- numeric(k)
#
## Perform k-fold cross-validation
#for(i in 1:k) {
#  # Split data into training and validation sets
#  train_fold <- training[-folds[[i]], ]
#  val_fold <- training[folds[[i]], ]
#  
#  # Train Random Forest on training fold
#  classifier_fold <- randomForest(
#    formula = formula,
#    data = train_fold,
#    ntree = 250,
#    mtry = 20,
#    importance = TRUE,
#    na.action = na.omit
#  )
#  
#  # Predict on validation fold
#  predictions <- predict(classifier_fold, newdata = val_fold)
#  
#  # Calculate performance metrics
#  rmse_values[i] <- sqrt(mean((val_fold$pai - predictions)^2, na.rm = TRUE))
#  mae_values[i] <- mean(abs(val_fold$pai - predictions), na.rm = TRUE)
#  
#  # Calculate R-squared
#  ss_res <- sum((val_fold$pai - predictions)^2, na.rm = TRUE)
#  ss_tot <- sum((val_fold$pai - mean(val_fold$pai, na.rm = TRUE))^2, na.rm = TRUE)
#  r2_values[i] <- 1 - (ss_res / ss_tot)
#  
#  cat(sprintf("Fold %d: RMSE = %.4f, MAE = %.4f, R² = %.4f\n", 
#              i, rmse_values[i], mae_values[i], r2_values[i]))
#}
#
## Print cross-validation results
#cat(sprintf("Mean R²: %.4f (+/- %.4f)\n", mean(r2_values), sd(r2_values)))
#cat("\n=== Cross-Validation Results ===\n")
#cat(sprintf("Mean RMSE: %.4f (+/- %.4f)\n", mean(rmse_values), sd(rmse_values)))
#cat(sprintf("Mean MAE: %.4f (+/- %.4f)\n", mean(mae_values), sd(mae_values)))
#
## Train final model on full training set
#classifier <- randomForest(
#  formula = formula,
#  data = training,
#  ntree = 250,
#  mtry = 20,
#  importance = TRUE,
#  na.action = na.omit
#)
#
#print(classifier)


#*************************************************************
# Test Model Accuracy - For Continuous Variables
#*************************************************************
library(pROC)
# Predict on testing data
predictions <- predict(classifier, testing)

# Calculate RMSE (Root Mean Square Error)
rmse <- sqrt(mean((predictions - testing$fhd)^2))
print(paste("RMSE:", round(rmse, 4)))

# Calculate R-squared
ss_res <- sum((testing$rh98 - predictions)^2)
ss_tot <- sum((testing$rh98 - mean(testing$rh98))^2)
r_squared <- 1 - (ss_res / ss_tot)
print(paste("R-squared:", round(r_squared, 4)))

# Calculate MAE (Mean Absolute Error)
mae <- mean(abs(predictions - testing$fhd))
print(paste("MAE:", round(mae, 4)))


# Plot observed vs predicted
plot(testing$fhd, predictions, 
     xlab = "Observed FHD", 
     ylab = "Predicted FHD",
     main = "Observed vs Predicted")
abline(0, 1, col = "red", lwd = 2)  # 1:1 line

# Calculate residuals
residuals <- testing$cover - predictions

# Plot residuals vs predicted
plot(predictions, residuals,
     xlab = "Predicted cover", 
     ylab = "Residuals",
     main = "Residuals")
abline(h=0, col = "red", lwd = 2)  # 1:1 line



residuals_sd<-sd(residuals)

# If testing is a matrix, convert to data frame
testing <- as.data.frame(testing)
training <- as.data.frame(training)

#*************************************************************
# Classify the entire image
#*************************************************************

# Predict across the entire raster
classifiedrf <- predict(rasterStack, classifier, type = "response")

#*************************************************************
# Visualize and Export Results
#*************************************************************

# Plot the classified raster
plot(classifiedrf, main = "fhd Classification")

# Export the classification
writeRaster(classifiedrf, 
            filename = "/Users/abarenb1/Documents/PhD/Chapter3/Variables/FHD_classified_2023_Mar12_2026.tif", 
            overwrite = TRUE)

# Optional: View variable importance
varImpPlot(classifier, main = "Variable Importance")
importance(classifier)

#*************************************************************
# Compare Lang Height to RH100
#*************************************************************
method <- "bilinear"

height<-rast("/Users/abarenb1/Documents/PhD/Chapter3/Variables/ETH_GlobalCanopyHeight_2020_10m_v1.tif")
# Project and resample to template
project(L2A, height, method = method)

rh98 <- L2A[, "rh98"]
rh100  <- L2A[, "rh100"]

samples <- terra::extract(height, rh98, df = TRUE, sp = TRUE)
samples$rh98 <- rh98$rh98

# Calculate R-squared
ss_res <- sum((samples$rh98 - samples$b1)^2)
ss_tot <- sum((samples$rh98 - mean(samples$rh98))^2)
r_squared <- 1 - (ss_res / ss_tot)
print(paste("R-squared:", round(r_squared, 4)))

plot(samples$b1,samples$rh98, 
     xlab = "RH98", 
     ylab = "Wall-to-Wall Height")
abline(0, 1, col = "red", lwd = 2)  # 1:1 line

#*************************************************************
# Saved Testing Of Raster
#*************************************************************

#####
#RH98
####
# Extract raster values at your sample points/polygons

method <- "bilinear"

# Reproject to match raster CRS
rh98 <- st_transform(rh98, crs = st_crs(crs(rasterStack)))

# ---------------------------------------------------------------
# 🔧 KEY FIX: Check NA coverage per layer at point locations
# ---------------------------------------------------------------
rh98_vect <- vect(rh98)

# Extract WITHOUT removing NAs first
samples_raw <- terra::extract(rasterStack, rh98_vect, df = TRUE)
samples_raw$rh98 <- rh98$rh98

cat("\nNA count per band at GEDI point locations:\n")
na_counts <- sapply(samples_raw[, bands], function(x) sum(is.na(x)))
print(na_counts)

cat("\nTotal points before NA removal:", nrow(samples_raw), "\n")
cat("Complete cases (no NA in ANY band):", sum(complete.cases(samples_raw)), "\n")

# ---------------------------------------------------------------
# 🔧 OPTION A: Keep only bands with low NA counts
#    Drop bands where >50% of points are NA
# ---------------------------------------------------------------
na_frac <- na_counts / nrow(samples_raw)
bad_bands <- names(na_frac[na_frac > 0.5])

if (length(bad_bands) > 0) {
  cat("\n⚠️  Dropping high-NA bands (>50% NA at point locations):\n")
  print(bad_bands)
  bands_use <- setdiff(bands, bad_bands)
} else {
  bands_use <- bands
}

cat("\nBands used for Random Forest:", length(bands_use), "\n")

# ---------------------------------------------------------------
# 🔧 OPTION B: Use na.action in randomForest instead of na.omit
#    This keeps more samples by imputing or ignoring NAs per tree
# ---------------------------------------------------------------
samples98 <- samples_raw[, c("ID", bands_use, "rh98")]
samples98 <- na.omit(samples98)   # Only remove rows with NA in KEPT bands

cat("\nSamples after NA removal:", nrow(samples98), "\n")

# ---------------------------------------------------------------
# Train/test split
# ---------------------------------------------------------------
set.seed(123)
samples98$random <- runif(nrow(samples98))

split <- 0.8
training98 <- samples98[samples98$random < split, ]
testing98  <- samples98[samples98$random >= split, ]

cat("Samples n =",  nrow(samples98),  "\n")
cat("Training n =", nrow(training98), "\n")
cat("Testing n =",  nrow(testing98),  "\n")

formula98 <- as.formula(paste("rh98 ~", paste(bands_use, collapse = " + ")))

# Train Random Forest
# ntree = 250 trees, mtry = 20 variables tried at each split
classifier98 <- randomForest(
  formula = formula98,
  data = training98,
  ntree = 500,
  mtry = 10,
  importance = TRUE,
  na.action = na.omit
)

predictions98 <- predict(classifier98, testing98)

p1<-ggplot()+geom_point(aes(testing98$rh98, predictions98,alpha = 1/20),show.legend = FALSE)+
  geom_abline(color="red",lwd=1)+labs(y = "Predicted",
                                      x = "Observed")+
  annotate("text", x=90, y=60, parse = TRUE,label= "R^2==0.57",
           color="blue",size=4)+
  theme_minimal()+theme(text = element_text(size = 12))

p1

classifiedrf <- predict(rasterStack, classifier98, type = "response")

p2<- plot(classifiedrf, main = "RH98 Classification")

library(patchwork)

p1 + p2 + plot_layout(ncol = 2, heights = c(1, 1))



###
#FHD
###

# Reproject to match raster CRS
fhd <- st_transform(fhd, crs = st_crs(crs(rasterStack)))

# ---------------------------------------------------------------
# 🔧 KEY FIX: Check NA coverage per layer at point locations
# ---------------------------------------------------------------
fhd_vect <- vect(fhd)

# Extract WITHOUT removing NAs first
samples_raw <- terra::extract(rasterStack, fhd_vect, df = TRUE)
samples_raw$fhd <- fhd$fhd

cat("\nNA count per band at GEDI point locations:\n")
na_counts <- sapply(samples_raw[, bands], function(x) sum(is.na(x)))
print(na_counts)

cat("\nTotal points before NA removal:", nrow(samples_raw), "\n")
cat("Complete cases (no NA in ANY band):", sum(complete.cases(samples_raw)), "\n")

# ---------------------------------------------------------------
# 🔧 OPTION A: Keep only bands with low NA counts
#    Drop bands where >50% of points are NA
# ---------------------------------------------------------------
na_frac <- na_counts / nrow(samples_raw)
bad_bands <- names(na_frac[na_frac > 0.5])

if (length(bad_bands) > 0) {
  cat("\n⚠️  Dropping high-NA bands (>50% NA at point locations):\n")
  print(bad_bands)
  bands_use <- setdiff(bands, bad_bands)
} else {
  bands_use <- bands
}

cat("\nBands used for Random Forest:", length(bands_use), "\n")

# ---------------------------------------------------------------
# 🔧 OPTION B: Use na.action in randomForest instead of na.omit
#    This keeps more samples by imputing or ignoring NAs per tree
# ---------------------------------------------------------------
samplesfhd <- samples_raw[, c("ID", bands_use, "fhd")]
samplesfhd <- na.omit(samplesfhd)   # Only remove rows with NA in KEPT bands

cat("\nSamples after NA removal:", nrow(samplesfhd), "\n")

# ---------------------------------------------------------------
# Train/test split
# ---------------------------------------------------------------
set.seed(123)
samplesfhd$random <- runif(nrow(samplesfhd))

split <- 0.8
trainingfhd <- samplesfhd[samplesfhd$random < split, ]
testingfhd  <- samplesfhd[samplesfhd$random >= split, ]

cat("Samples n =",  nrow(samplesfhd),  "\n")
cat("Training n =", nrow(trainingfhd), "\n")
cat("Testing n =",  nrow(testingfhd),  "\n")

formulafhd <- as.formula(paste("fhd ~", paste(bands_use, collapse = " + ")))

# Train Random Forest
# ntree = 250 trees, mtry = 20 variables tried at each split
classifierfhd <- randomForest(
  formula = formulafhd,
  data = trainingfhd,
  ntree = 500,
  mtry = 10,
  importance = TRUE,
  na.action = na.omit
)

predictionsfhd <- predict(classifierfhd, testingfhd)

p2<-ggplot()+geom_point(aes(testingfhd$fhd, predictionsfhd,alpha = 1/10),show.legend = FALSE)+
  geom_abline(color="red",lwd=1)+labs(y = "Predicted",
                                      x = "Observed")+
  annotate("text", x=2.9, y=3.5, parse = TRUE,label= "R^2==0.55",
           color="blue",size=4)+
  theme_minimal()+theme(text = element_text(size = 12))


p2

ss_res <- sum((testing98$rh98 - predictions98)^2)
ss_tot <- sum((testing98$rh98 - mean(testing98$rh98))^2)
r_squared <- 1 - (ss_res / ss_tot)
print(paste("R-squared:", round(r_squared, 4)))

ss_res <- sum((testingfhd$fhd - predictions)^2)
ss_tot <- sum((testingfhd$fhd - mean(testingfhd$fhd))^2)
r_squared <- 1 - (ss_res / ss_tot)
print(paste("R-squared:", round(r_squared, 4)))


library(patchwork)

p1 + p2 + plot_layout(ncol = 2, heights = c(1, 1))

##*************************************************************
## Plot Rasters of RH98 and FHD
##*************************************************************
rh98_Rast<-rast("/Users/abarenb1/Documents/PhD/Chapter3/Variables/RH98_classified_2023_Feb23_2026.tif")
fhd_Rast<-rast("/Users/abarenb1/Documents/PhD/Chapter3/Variables/FHD_classified_2023.tif")

rh98_df <-
  as.data.frame(rh98_Rast, xy = TRUE) %>%
  #--- remove cells with NA for any of the layers ---#
  na.omit() 

rh98Plot<-ggplot(data=rh98_df)+geom_raster(aes(x = x, y = y, fill = lyr1))+
  scale_fill_viridis_c(limits = c(2,50)) +theme_minimal()+theme(axis.title = element_blank())+
  labs(fill = "RH98")

rh98Plot

FHD_df <-
  as.data.frame(fhd_Rast, xy = TRUE) %>%
  #--- remove cells with NA for any of the layers ---#
  na.omit() 

FHDPlot<-ggplot(data=FHD_df)+geom_raster(aes(x = x, y = y, fill = lyr1))+
  scale_fill_viridis_c(option="cividis",limits = c(1.75,3.31)) +theme_minimal()+theme(axis.title = element_blank())+
  labs(fill = "FHD")

FHDPlot

p1+p2+rh98Plot + FHDPlot + plot_layout(ncol = 2, heights = c(1, 1))

##################################
#Compare values in and out of park
##################################
rh98_int<-rast("/Users/abarenb1/Documents/PhD/Chapter3/RH98_Masked.tif")
rh98_out<-rast("/Users/abarenb1/Documents/PhD/Chapter3/RH98_Difference.tif") 

fhd_int<-rast("/Users/abarenb1/Documents/PhD/Chapter3/FHD_Masked.tif") 
fhd_out<-rast("/Users/abarenb1/Documents/PhD/Chapter3/FHD_Difference.tif")

global(rh98_int, fun = 'mean',na.rm=TRUE)
global(rh98_int, fun = 'sd',na.rm=TRUE)

global(rh98_out, fun = 'mean',na.rm=TRUE)
global(rh98_out, fun = 'sd',na.rm=TRUE)

global(fhd_int, fun = 'mean',na.rm=TRUE)
global(fhd_int, fun = 'sd',na.rm=TRUE)

global(fhd_out, fun = 'mean',na.rm=TRUE)
global(fhd_out, fun = 'sd',na.rm=TRUE)


##*************************************************************
## Test with XGBoost
##*************************************************************
#library(xgboost)
#set.seed(123)  # For reproducibility
#
#samples<-samples[2:69]
#
##split into training (80%) and testing set (20%)
#parts = createDataPartition(samples$fhd, p = .8, list = F)
#train = samples[parts, ]
#test = samples[-parts, ]
#
##define predictor and response variables in training set
#train_x =as.matrix(train[, -13])
#train_y = train[,13]
#
##define predictor and response variables in testing set
#test_x = as.matrix(test[, -13])
#test_y = test[, 13]
#
##define final training and testing sets
#xgb_train = xgb.DMatrix(data = train_x, label = train_y)
#xgb_test = xgb.DMatrix(data = test_x, label = test_y)
#
##define watchlist
#watchlist = list(train=xgb_train, test=xgb_test)
#
## Model Building: XGBoost
#param_list = list(
#  objective = "reg:squarederror",
#  eta = 0.01,
#  gamma = 1,
#  max_depth = 6,
#  subsample = 0.8,
#  colsample_bytree = 0.5
#)
#
##fit XGBoost model and display training and testing data at each round
#model <- xgb.train(data = xgb_train, params=param_list, evals=watchlist, nrounds = 250)
#
#xgbpred <- function(model, data) {
#  predict(model, newdata=as.matrix(data))
#}
#
#p <- predict(rasterStack, model=model, fun=xgbpred)
#plot(p)
#
#fit<-predict(model, xgb_test)
#
## Plot observed vs predicted
#plot(test_x[,67], fit, 
#     xlab = "Observed fhd", 
#     ylab = "Predicted fhd",
#     main = "Observed vs Predicted")
#abline(0, 1, col = "red", lwd = 2)  # 1:1 line
#
## Calculate residuals
#residuals <- testing$fhd - predictions
#
## Plot residuals vs predicted
#plot(predictions, residuals,
#     xlab = "Predicted fhd", 
#     ylab = "Residuals",
#     main = "Residuals")
#abline(h=0, col = "red", lwd = 2)  # 1:1 line
#
#residuals_sd<-sd(residuals)

##*************************************************************
## Test with Cross Validation
##*************************************************************
#
#data <- samples
#
## in this cross validation example, we use the iris data set to 
## predict the Sepal Length from the other variables in the dataset 
## with the random forest model 
#
#k = 5 #Folds
#
## sample from 1 to k, nrow times (the number of observations in the data)
#data$id <- sample(1:k, nrow(data), replace = TRUE)
#list <- 1:k
#
## prediction and testset data frames that we add to with each iteration over
## the folds
#
#prediction <- data.frame()
#testsetCopy <- data.frame()
#
##Creating a progress bar to know the status of CV
#progress.bar <- create_progress_bar("text")
#progress.bar$init(k)
#
#formula <- as.formula(paste("fhd ~", paste(bands, collapse = " + ")))
#
#for (i in 1:k){
#  # remove rows with id i from dataframe to create training set
#  # select rows with id i to create test set
#  trainingset <- subset(data, id %in% list[-i])
#  testset <- subset(data, id %in% c(i))
#  
#  # run a random forest model
#  mymodel <- randomForest(formula = formula, data = trainingset, ntree = 100)
#  
#  # remove response column 1, Sepal.Length
#  temp <- as.data.frame(predict(mymodel, testset[,-1]))
#  # append this iteration's predictions to the end of the prediction data frame
#  prediction <- rbind(prediction, temp)
#  
#  # append this iteration's test set to the test set copy data frame
#  # keep only the Sepal Length Column
#  testsetCopy <- rbind(testsetCopy, as.data.frame(testset[,1]))
#  
#  progress.bar$step()
#}
#
## add predictions and actual Sepal Length values
#result <- cbind(prediction, testsetCopy[, 1])
#names(result) <- c("Predicted", "Actual")
#result$Difference <- abs(result$Actual - result$Predicted)
#
## As an example use Mean Absolute Error as Evalution 
#summary(result$Difference)
#
## Plot observed vs predicted
#plot(result$Actual, result$Predicted, 
#     xlab = "Observed fhd", 
#     ylab = "Predicted fhd",
#     main = "Observed vs Predicted")
#abline(0, 1, col = "red", lwd = 2)  # 1:1 line
#
## Calculate residuals
#residuals <- result$Actual - result$Predicted
#
## Plot residuals vs predicted
#plot(result$Predicted, residuals,
#     xlab = "Predicted fhd", 
#     ylab = "Residuals",
#     main = "Residuals")
#abline(h=0, col = "red", lwd = 2)  # 1:1 line
#
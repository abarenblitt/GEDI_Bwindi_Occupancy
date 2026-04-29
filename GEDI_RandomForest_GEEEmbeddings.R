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
# These are the rasters you are using as input bands to the Random Forest model
# ---------------------------------------------------------------
matching_tifs <- c("ETH_GlobalCanopyHeight_2020_10m_v1", #Global Canopy from Lang et al
                   "WSCI",                               #Waveform Structural Complexity Index
                   "HLSLandsat",             
                   "HLSSentinel", 
                   "SentinelSAR",
                   "GEDI_biomass_filled")                #Biomass from Armston personal comms

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
# Clean band names BEFORE stacking 
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
# Load GEDI data, these are the GEDI shots you will use for training and testing
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
# Check NA coverage per layer at point locations
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

#****************************************************************************
# Code below can be used to apply weights to values separated into quantiles
#****************************************************************************
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

#*************************************************************
# Train Random Forest Model
#*************************************************************
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

#****************************************************************************
# Alternatively, you can use a k-fold method to try to increase accuracy
#****************************************************************************
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


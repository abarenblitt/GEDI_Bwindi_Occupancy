library(terra)
library(randomForest)
library(caret)
library(sf)
library(dplyr)

#Merge output from GEE for Embeddings 2023

ras1<-rast("/Users/abarenb1/Documents/PhD/Chapter3/AlphaEmbeddings/Embeddings2023-0000000000-0000000000.tif")
ras2<-rast("/Users/abarenb1/Documents/PhD/Chapter3/AlphaEmbeddings/Embeddings2023-0000000000-0000003072.tif")
ras3<-rast("/Users/abarenb1/Documents/PhD/Chapter3/AlphaEmbeddings/Embeddings2023-0000003072-0000000000.tif")
ras4<-rast("/Users/abarenb1/Documents/PhD/Chapter3/AlphaEmbeddings/Embeddings2023-0000003072-0000003072.tif")



s<-sprc(ras1,ras2,ras3,ras4)
m<-merge(s)

writeRaster(m, filename="/Users/abarenb1/Documents/PhD/Chapter3/AlphaEmbeddings/Embeddings2023.tif")

#Pull full raster
embed<-rast("/Users/abarenb1/Documents/PhD/Chapter3/AlphaEmbeddings/Embeddings2023.tif")

#Pull GEDI Data
L2A<-st_read("/Users/abarenb1/Documents/PhD/Chapter3/Variables/Bwindi_L2A.shp")
L2B<-st_read("/Users/abarenb1/Documents/PhD/Chapter3/Variables/Bwindi_L2B.shp")

L2B <- L2B %>% rename_at('fhd_normal', ~'fhd')


fhd <- L2B[, "fhd"]
cover <- L2B[, "cover"]
pai <- L2B[, "pai"]
rh98 <- L2A[, "rh98"]

#RANDOM FOREST
################################################################
#Choose Variable


# Extract band names
bands <- names(embed)
print(bands)

# Extract raster values at your sample points/polygons
samples <- extract(embed, rh98, df = TRUE, sp = TRUE)
samples$rh98 <- rh98$rh98

# Remove any NA values
samples <- na.omit(samples)

# Add random column for train/test split
set.seed(123)  # For reproducibility
samples$random <- runif(nrow(samples))

#*************************************************************
# Split into training and testing
#*************************************************************

split <- 0.9
training <- samples[samples$random < split, ]
testing <- samples[samples$random >= split, ]

# Print sample sizes
cat('Samples n =', nrow(samples), '\n')
cat('Training n =', nrow(training), '\n')
cat('Testing n =', nrow(testing), '\n')

#*************************************************************
# Train Random Forest Model
#*************************************************************

# Prepare formula (all bands predicting pai)
formula <- as.formula(paste("rh98 ~", paste(bands, collapse = " + ")))

# Train Random Forest
# ntree = 250 trees, mtry = 20 variables tried at each split
classifier <- randomForest(
  formula = formula,
  data = training,
  ntree = 250,
  mtry = 20,
  importance = TRUE,
  na.action = na.omit
)

print(classifier)

#*************************************************************
# Test Model Accuracy
#*************************************************************

#*************************************************************
# Test Model Accuracy - For Continuous Variables
#*************************************************************

# Predict on testing data
predictions <- predict(classifier, testing)

# Calculate RMSE (Root Mean Square Error)
rmse <- sqrt(mean((predictions - testing$rh98)^2))
print(paste("RMSE:", round(rmse, 4)))

# Calculate R-squared
ss_res <- sum((testing$rh98 - predictions)^2)
ss_tot <- sum((testing$rh98 - mean(testing$rh98))^2)
r_squared <- 1 - (ss_res / ss_tot)
print(paste("R-squared:", round(r_squared, 4)))

# Calculate MAE (Mean Absolute Error)
mae <- mean(abs(predictions - testing$rh98))
print(paste("MAE:", round(mae, 4)))

# Plot observed vs predicted
plot(testing$rh98, predictions, 
     xlab = "Observed rh98", 
     ylab = "Predicted rh98",
     main = "Observed vs Predicted")
abline(0, 1, col = "red", lwd = 2)  # 1:1 line

#*************************************************************
# Classify the entire image
#*************************************************************

# Predict across the entire raster
classifiedrf <- predict(embed, classifier, type = "response")

#*************************************************************
# Visualize and Export Results
#*************************************************************

# Plot the classified raster
plot(classifiedrf, main = "rh98 Classification")

# Export the classification
writeRaster(classifiedrf, 
            filename = "/Users/abarenb1/Documents/PhD/Chapter3/Variables/RH98_classified_2023.tif", 
            overwrite = TRUE)

# Optional: View variable importance
varImpPlot(classifier, main = "Variable Importance")
importance(classifier)

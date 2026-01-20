library(terra)
library(tidyverse)
library(viridis)

# ============================================================
# STEP 1: Build Raster Stack of Covariates
# ============================================================
#Build Stack of Covariate Rasters
matching_tifs <- c("GEDI_biomass_filled","Dist_Wat_25m_Uganda_250mDist","Precipitation_WorldClim",
                   "WSCI","Bwindi_Heights_2mMaxarResample",
                   "Bwindi_Slope_2mMaxarResample","RH98_classified_2023","FHD_classified_2023",
                   "PAI_classified_2023","Cover_classified_2023","Population_Dens_Bwindi",
                   #"TMIN_WorldClim","TAVG_WorldClim","TMAX_WorldClim",
                   "Bwindi_tsri")
#,"Bwindi_tpi")

rasterList <- lapply(1:length(matching_tifs), function(x) {
  r <- rast(paste("/Users/abarenb1/Documents/PhD/Chapter3/Variables/",matching_tifs[[x]],".tif", sep=""))
  
  # For WSCI raster, select only the "Mean WSCI estimate" layer
  if(matching_tifs[[x]] == "WSCI") {
    # Option 1: If you know the layer name
    r <- r[["Mean WSCI estimate"]]
    
    # Option 2: If you know it's a specific band number (e.g., band 1)
    # r <- r[[1]]
  }
  
  return(r)
})

names(rasterList) <- c("biomass", "dist_water", "precip","wsci","height","slope","rh98","fhd","pai","cover","pop","tsri")

# Use the first raster as the template
template <- rasterList[[1]]

# Resample all others to match the template
rasterList_aligned <- lapply(rasterList, function(r) 
  project(r, template, method = "bilinear"))

# Now stack them
rasterStack <- do.call(c, rasterList_aligned)

#Bwindi BBOX
bwindi <- st_read("/Users/abarenb1/Documents/JGI/Bwindi Park Boundary/Bwindi_National_Park_Boundary.shp")

# ============================================================
# STEP 2: Standardize Rasters Using Training Data Parameters
# ============================================================

# IMPORTANT: Use the SAME standardization parameters from your training data
# Extract the means and SDs used during model fitting

standardization_params <- data.frame(
  covariate = c("biomass", "dist_water", "precip", "height", "slope",
                "rh98", "fhd", "pai", "cover", "pop",
                #"tmin", "tavg", "tmax", 
                "tsri"),
  mean = c(
    mean(site_covs$GEDI_biomass_filled, na.rm = TRUE),
    mean(site_covs$Dist_Wat_25m_Uganda_250mDist, na.rm = TRUE),
    mean(site_covs$Precipitation_WorldClim, na.rm = TRUE),
    mean(site_covs$Bwindi_Heights_2mMaxarResample, na.rm = TRUE),
    mean(site_covs$Bwindi_Slope_2mMaxarResample, na.rm = TRUE),
    mean(site_covs$RH98_classified_2023),
    mean(site_covs$FHD_classified_2023),
    mean(site_covs$PAI_classified_2023),
    mean(site_covs$Cover_classified_2023),
    mean(site_covs$Population_Dens_Bwindi),
   #mean(site_covs$TMIN_WorldClim)
   #mean(site_covs$TAVG_WorldClim)
   #mean(site_covs$TMAX_WorldClim)
    mean(site_covs$Bwindi_tsri)
  ),
  sd = c(
    sd(site_covs$GEDI_biomass_filled, na.rm = TRUE),
    sd(site_covs$Dist_Wat_25m_Uganda_250mDist, na.rm = TRUE),
    sd(site_covs$Precipitation_WorldClim, na.rm = TRUE),
    sd(site_covs$Bwindi_Heights_2mMaxarResample, na.rm = TRUE),
    sd(site_covs$Bwindi_Slope_2mMaxarResample, na.rm = TRUE),
    sd(site_covs$RH98_classified_2023),
    sd(site_covs$FHD_classified_2023),
    sd(site_covs$PAI_classified_2023),
    sd(site_covs$Cover_classified_2023),
    sd(site_covs$Population_Dens_Bwindi),
    #sd(site_cov$TMIN_WorldClim)
    #sd(site_cov$TAVG_WorldClim)
    #sd(site_cov$TMAX_WorldClim)
    sd(site_covs$Bwindi_tsri)
  )
)

# Save for future use
write.csv(standardization_params, "standardization_parameters.csv", row.names = FALSE)

# Standardize each layer in the stack
covariate_stack_std <- rast(rasterStack)

for(i in 1:nrow(standardization_params)) {
  cov_name <- standardization_params$covariate[i]
  cov_mean <- standardization_params$mean[i]
  cov_sd <- standardization_params$sd[i]
  
  covariate_stack_std[[cov_name]] <- (rasterStack[[cov_name]] - cov_mean) / cov_sd
}

# ============================================================
# STEP 3: Extract Model Coefficients
# ============================================================

# Community-level (average) predictions
community_intercept <- jags_out$mean$mu.lpsi

community_coefs <- data.frame(
  covariate = c("biomass", "dist_water", "precip", "height", "slope",
                "rh98", "fhd", "pai", "cover", "pop",
               # "tmin", "tavg", "tmax", 
               "tsri"),
  coefficient = c(
    jags_out$mean$mu.beta.biomass,
    jags_out$mean$mu.beta.dist_water,
    jags_out$mean$mu.beta.precip,
    jags_out$mean$mu.beta.height,
    jags_out$mean$mu.beta.slope,
    jags_out$mean$mu.beta.rh98,
    jags_out$mean$mu.beta.fhd,
    jags_out$mean$mu.beta.pai,
    jags_out$mean$mu.beta.cover,
    jags_out$mean$mu.beta.pop,
    #jags_out$mean$mu.beta.tmin,
    #jags_out$mean$mu.beta.tavg,
    #jags_out$mean$mu.beta.tmax,
    jags_out$mean$mu.beta.tsri
  )
)

print(community_coefs)

# ============================================================
# STEP 4: Predict Community-Level Occupancy (CORRECTED)
# ============================================================

# Method 1: Using app() function (most reliable for terra)
predict_community <- function(x, intercept, coefs_df) {
  # x is a matrix where each row is a pixel and columns are covariates
  
  # Start with intercept
  logit_psi <- intercept
  
  # Add each covariate contribution
  for(i in 1:nrow(coefs_df)) {
    cov_name <- coefs_df$covariate[i]
    cov_coef <- coefs_df$coefficient[i]
    
    # Find column index for this covariate
    col_idx <- which(names(covariate_stack_std) == cov_name)
    
    if(length(col_idx) > 0) {
      logit_psi <- logit_psi + cov_coef * x[, col_idx]
    }
  }
  
  # Convert to probability
  psi <- 1 / (1 + exp(-logit_psi))
  
  return(psi)
}

# Apply prediction function
community_occupancy <- app(covariate_stack_std, 
                           fun = function(x) {
                             predict_community(x, community_intercept, community_coefs)
                           })

names(community_occupancy) <- "community_occupancy_prob"

# Also calculate logit scale
community_logit <- app(covariate_stack_std,
                       fun = function(x) {
                         logit_psi <- community_intercept
                         for(i in 1:nrow(community_coefs)) {
                           cov_name <- community_coefs$covariate[i]
                           cov_coef <- community_coefs$coefficient[i]
                           col_idx <- which(names(covariate_stack_std) == cov_name)
                           if(length(col_idx) > 0) {
                             logit_psi <- logit_psi + cov_coef * x[, col_idx]
                           }
                         }
                         return(logit_psi)
                       })

names(community_logit) <- "community_occupancy_logit"

# Save predictions
writeRaster(community_occupancy, "community_occupancy_prediction.tif", overwrite = TRUE)
writeRaster(community_logit, "community_occupancy_logit.tif", overwrite = TRUE)

# Plot
plot(community_occupancy, 
     main = "Community-Level Occupancy Probability",
     col = viridis(100))

# ============================================================
# STEP 5: Extract Species-Specific Coefficients
# ============================================================
bird_data <- read.csv("/Users/abarenb1/Documents/PhD/Chapter3/Bwindi_Birds_Variables_Jan2026.csv")
sites <- sort(unique(bird_data$id))
species_list <- sort(unique(bird_data$Total))
species_names <- species_list

n_species <- length(species_names)

# Extract species-specific intercepts
species_intercepts <- jags_out$mean$lpsi

# Extract species-specific coefficients for each covariate
species_coefs <- data.frame(
  species = species_names,
  intercept = species_intercepts,
  biomass = jags_out$mean$beta.biomass,
  dist_water = jags_out$mean$beta.dist_water,
  precip = jags_out$mean$beta.precip,
  height = jags_out$mean$beta.height,
  slope = jags_out$mean$beta.slope,
  rh98 = jags_out$mean$beta.rh98,
  fhd = jags_out$mean$beta.fhd,
  pai = jags_out$mean$beta.pai,
  cover = jags_out$mean$beta.cover,
  pop = jags_out$mean$beta.pop,
  tsri = jags_out$mean$beta.tsri
)

# Save species coefficients
write.csv(species_coefs, "species_specific_coefficients.csv", row.names = FALSE)

print(species_coefs)

# ============================================================
# STEP 6: Predict Species-Specific Occupancy
# ============================================================

# Function to predict occupancy for a single species
predict_species <- function(x, species_row, covariate_stack_names) {
  # x is a matrix where each row is a pixel and columns are covariates
  
  # Start with species-specific intercept
  logit_psi <- species_row$intercept
  
  # Add each covariate contribution
  covariate_names <- c("biomass", "dist_water", "precip", "height", "slope",
                       "rh98", "fhd", "pai", "cover", "pop", "tsri")
  
  for(cov_name in covariate_names) {
    cov_coef <- species_row[[cov_name]]
    
    # Find column index for this covariate in the raster stack
    col_idx <- which(covariate_stack_names == cov_name)
    
    if(length(col_idx) > 0) {
      logit_psi <- logit_psi + cov_coef * x[, col_idx]
    }
  }
  
  # Convert to probability
  psi <- 1 / (1 + exp(-logit_psi))
  
  return(psi)
}

# Create a list to store species predictions
species_occupancy_list <- list()
species_logit_list <- list()

# Get covariate names from stack
stack_names <- names(covariate_stack_std)

# Loop through each species and create prediction maps
for(i in 1:n_species) {
  
  cat("Processing species", i, "of", n_species, ":", species_names[i], "\n")
  
  # Get coefficients for this species
  species_row <- species_coefs[i, ]
  
  # Predict occupancy probability
  species_occ <- app(covariate_stack_std,
                     fun = function(x) {
                       predict_species(x, species_row, stack_names)
                     })
  
  names(species_occ) <- paste0(species_names[i], "_occupancy_prob")
  
  # Predict logit scale
  species_log <- app(covariate_stack_std,
                     fun = function(x) {
                       logit_psi <- species_row$intercept
                       covariate_names <- c("biomass", "dist_water", "precip", "height", "slope",
                                            "rh98", "fhd", "pai", "cover", "pop", "tsri")
                       for(cov_name in covariate_names) {
                         cov_coef <- species_row[[cov_name]]
                         col_idx <- which(stack_names == cov_name)
                         if(length(col_idx) > 0) {
                           logit_psi <- logit_psi + cov_coef * x[, col_idx]
                         }
                       }
                       return(logit_psi)
                     })
  
  names(species_log) <- paste0(species_names[i], "_occupancy_logit")
  
  # Store in lists
  species_occupancy_list[[i]] <- species_occ
  species_logit_list[[i]] <- species_log
  
  # Save individual rasters
  writeRaster(species_occ, 
              paste0("/Users/abarenb1/Documents/PhD/Chapter3/species_predictions/", species_names[i], "_occupancy_prob.tif"),
              overwrite = TRUE)
  writeRaster(species_log,
              paste0("/Users/abarenb1/Documents/PhD/Chapter3/species_predictions/", species_names[i], "_occupancy_logit.tif"),
              overwrite = TRUE)
}

# Create multi-layer rasters
all_species_occupancy <- rast(species_occupancy_list)
all_species_logit <- rast(species_logit_list)

# Save combined rasters
writeRaster(all_species_occupancy, "/Users/abarenb1/Documents/PhD/Chapter3/species_predictions/all_species_occupancy_probabilities.tif", overwrite = TRUE)
writeRaster(all_species_logit, "/Users/abarenb1/Documents/PhD/Chapter3/species_predictions/all_species_occupancy_logit.tif", overwrite = TRUE)

# ============================================================
# STEP 7: Visualize Species Maps
# ============================================================

# Create output directory if it doesn't exist
dir.create("/Users/abarenb1/Documents/PhD/Chapter3/species_predictions/species_maps", showWarnings = FALSE)

# Plot individual species maps
pdf("species_maps/all_species_occupancy_maps.pdf", width = 11, height = 8.5)
for(i in 1:n_species) {
  plot(species_occupancy_list[[i]], 
       main = paste("Occupancy Probability:", species_names[i]),
       col = viridis(100),
       range = c(0, 1))
  plot(st_geometry(bwindi), add = TRUE, border = "black", lwd = 2)
}
dev.off()

# Create a multi-panel plot for all species
n_cols <- ceiling(sqrt(n_species))
n_rows <- ceiling(n_species / n_cols)

pdf("species_maps/all_species_grid.pdf", width = 14, height = 10)
par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 2, 2))
for(i in 1:n_species) {
  plot(species_occupancy_list[[i]], 
       main = species_names[i],
       col = viridis(100),
       range = c(0, 1),
       legend = FALSE)
}
dev.off()

# ============================================================
# STEP 8: Calculate Species Richness
# ============================================================

# Sum occupancy probabilities across species (expected richness)
species_richness <- app(all_species_occupancy, sum)
names(species_richness) <- "expected_species_richness"

# Save
writeRaster(species_richness, "species_richness_map.tif", overwrite = TRUE)

# Plot
pdf("species_maps/species_richness.pdf", width = 10, height = 8)
plot(species_richness,
     main = "Expected Species Richness",
     col = viridis(100))
plot(st_geometry(bwindi), add = TRUE, border = "black", lwd = 2)


# ============================================================
# STEP 9: Extract Uncertainty (Optional - if you saved SD/CI)
# ============================================================

# If you have credible intervals from JAGS output:
if("sd" %in% names(jags_out)) {
  
  # Extract species-specific SDs
  species_sd <- data.frame(
    species = species_names,
    intercept_sd = jags_out$sd$lpsi,
    biomass_sd = jags_out$sd$beta.biomass,
    dist_water_sd = jags_out$sd$beta.dist_water,
    precip_sd = jags_out$sd$beta.precip,
    height_sd = jags_out$sd$beta.height,
    slope_sd = jags_out$sd$beta.slope,
    rh98_sd = jags_out$sd$beta.rh98,
    fhd_sd = jags_out$sd$beta.fhd,
    pai_sd = jags_out$sd$beta.pai,
    cover_sd = jags_out$sd$beta.cover,
    pop_sd = jags_out$sd$beta.pop,
    tsri_sd = jags_out$sd$beta.tsri
  )
  
  write.csv(species_sd, "species_coefficients_uncertainty.csv", row.names = FALSE)
  
  # You could also map coefficient of variation (CV) for each species
  species_cv_list <- list()
  
  for(i in 1:n_species) {
    # Calculate approximate SE of prediction on logit scale
    # This is a simplified approach - full uncertainty propagation would be more complex
    
    cv <- species_sd$intercept_sd[i] / abs(species_coefs$intercept[i] + 0.001)
    
    species_cv_list[[i]] <- species_occupancy_list[[i]] * 0 + cv
    names(species_cv_list[[i]]) <- paste0(species_names[i], "_CV")
  }
  
  all_species_cv <- rast(species_cv_list)
  writeRaster(all_species_cv, "all_species_uncertainty.tif", overwrite = TRUE)
}

# ============================================================
# STEP 10: Summary Statistics
# ============================================================

# Extract summary statistics for each species
species_summary <- data.frame(
  species = species_names,
  mean_occupancy = sapply(species_occupancy_list, function(x) global(x, "mean", na.rm = TRUE)$mean),
  sd_occupancy = sapply(species_occupancy_list, function(x) global(x, "sd", na.rm = TRUE)$sd),
  min_occupancy = sapply(species_occupancy_list, function(x) global(x, "min", na.rm = TRUE)$min),
  max_occupancy = sapply(species_occupancy_list, function(x) global(x, "max", na.rm = TRUE)$max),
  median_occupancy = sapply(species_occupancy_list, function(x) global(x, "median", na.rm = TRUE)$median)
)

write.csv(species_summary, "species_occupancy_summary.csv", row.names = FALSE)

print(species_summary)

# Plot species comparison
pdf("species_maps/species_occupancy_comparison.pdf", width = 12, height = 6)
ggplot(species_summary, aes(x = reorder(species, mean_occupancy), y = mean_occupancy)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_occupancy - sd_occupancy, 
                    ymax = mean_occupancy + sd_occupancy),
                width = 0.2) +
  coord_flip() +
  labs(x = "Species", y = "Mean Occupancy Probability",
       title = "Mean Occupancy Probability by Species") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))
dev.off()

# ============================================================
# STEP 11: Create Interactive Comparison
# ============================================================

# Compare top 5 and bottom 5 species
species_summary <- species_summary %>%
  arrange(desc(mean_occupancy))

top5_species <- species_summary$species[1:min(5, n_species)]
bottom5_species <- species_summary$species[max(1, n_species-4):n_species]

pdf("species_maps/top_bottom_species.pdf", width = 14, height = 10)
par(mfrow = c(2, 5), mar = c(2, 2, 3, 2))

for(sp in top5_species) {
  idx <- which(species_names == sp)
  plot(species_occupancy_list[[idx]],
       main = paste("High:", sp),
       col = viridis(100),
       range = c(0, 1),
       legend = FALSE)
}

for(sp in bottom5_species) {
  idx <- which(species_names == sp)
  plot(species_occupancy_list[[idx]],
       main = paste("Low:", sp),
       col = viridis(100),
       range = c(0, 1),
       legend = FALSE)
}
dev.off()

cat("\n=== Species Mapping Complete ===\n")
cat("Total species processed:", n_species, "\n")
cat("Output files saved to:\n")
cat("  - species_predictions/ (individual GeoTIFFs)\n")
cat("  - species_maps/ (PDF visualizations)\n")
cat("  - all_species_occupancy_probabilities.tif\n")
cat("  - species_richness_map.tif\n")
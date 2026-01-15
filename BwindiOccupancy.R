#STEP 1: Process Environmental Variables with Bird Data
##############################################################################################
library(terra)
library(sf)

#Extract Environmental Variables Intersecting with Bird Sightings
birds<-st_read("/Users/abarenb1/Documents/PhD/Chapter3/Bwindi Endemic and threatened birds/Species_with_Coordinates_ESPG4326.shp")

matching_tifs <- c("GEDI_biomass_filled","Dist_Wat_25m_Uganda_250mDist","Precipitation_WorldClim",
                   "WSCI","Bwindi_Heights_2mMaxarResample","Bwindi_2mFilledDSM_MaxarResample",
                   "Bwindi_Slope_2mMaxarResample","RH98_classified_2023","FHD_classified_2023",
                   "PAI_classified_2023","Cover_classified_2023","Population_Dens_Bwindi",
                   "TMIN_WorldClim","TAVG_WorldClim","TMAX_WorldClim","Bwindi_tsri","Bwindi_tpi")
  
for (j in 1:length(matching_tifs)){
  print(matching_tifs[j])
  ras <- rast(paste("/Users/abarenb1/Documents/PhD/Chapter3/Variables/",matching_tifs[j],".tif", sep=""))
  ras_ex <- terra::extract(ras, birds, method="simple", factors=F)
  ras_ex[is.na(ras_ex)] <- 0
  nm <- names(ras)
  
  # Special handling for WSCI - specify which layer to use
  if(matching_tifs[j] == "WSCI"){
    birds[[matching_tifs[j]]] <- ras_ex[, "Mean WSCI estimate"]
  } else {
    birds[[matching_tifs[j]]] <- ras_ex[, nm] 
  }
}


st_drop_geometry(birds)


write.csv(birds,"/Users/abarenb1/Documents/PhD/Chapter3/Bwindi_Birds_Variables_Jan2026.csv")

#STEP 2: Examine Variable Correlation
##############################################################################################
# Load required packages
library(tidyverse)

# Read in your data
bird_data <- read.csv("/Users/abarenb1/Documents/PhD/Chapter3/Bwindi_Birds_Variables_Jan2026.csv")

#Examine Variable Correlation
library(corrplot)
data_subset <- bird_data[29:45]
# Remove columns with zero variance
data_subset <- data_subset[, sapply(data_subset, function(x) sd(x, na.rm = TRUE) > 0)]

# Check how many columns remain
cat("Columns after removing constant variables:", ncol(data_subset), "\n")

# Now calculate correlation
x <- cor(data_subset, use = "complete.obs")

colnames(x) <- substr(colnames(x), 1, 14)
rownames(x) <- substr(rownames(x), 1, 14)

corrplot(x, type="upper", order="hclust")

#STEP 3: Format Data in Correct Format for JAGS
##############################################################################################

# Prepare data for JAGS
# Extract response variable (bird presence/absence or counts)
sites <- sort(unique(bird_data$id))
species_list <- sort(unique(bird_data$Total))


# Clean the data
bird_data_clean <- bird_data %>%
  mutate(
    species = trimws(as.character(Total)),
    SiteKey.x = trimws(as.character(id))
  )

# Get unique species and sites FROM THE ACTUAL DATA
species_list <- sort(unique(bird_data_clean$species))
sites <- sort(unique(bird_data_clean$SiteKey.x))

n_species <- length(species_list)
n_sites <- length(sites)
n_visits <- 8

cat("Number of species:", n_species, "\n")
cat("Number of sites:", n_sites, "\n")
cat("Total rows in data:", nrow(bird_data_clean), "\n")
cat("Expected rows (species x sites):", n_species * n_sites, "\n")

# Check if we have all combinations
actual_combos <- bird_data_clean %>%
  distinct(species, SiteKey.x) %>%
  nrow()

cat("Actual species-site combinations:", actual_combos, "\n")

# If not all combinations exist, that's OK - we'll fill with zeros
if(actual_combos < n_species * n_sites) {
  cat("\nNote: Not all species are present at all sites (this is normal)\n")
}

# Create detection history array - initialize with 0
y <- array(0, dim = c(n_species, n_sites, n_visits))

# Fill array directly from bird_data_clean
visit_cols <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8")

for(i in 1:nrow(bird_data_clean)) {
  # Find indices
  sp <- bird_data_clean$species[i]
  st <- bird_data_clean$SiteKey.x[i]
  
  species_idx <- which(species_list == sp)
  site_idx <- which(sites == st)
  
  if(length(species_idx) == 1 && length(site_idx) == 1) {
    # Extract visits and convert to numeric
    visits <- as.numeric(bird_data_clean[i, visit_cols])
    
    # Replace NAs with 0 in visits
    visits[is.na(visits)] <- 0
    
    # Fill array
    y[species_idx, site_idx, ] <- visits
  } else {
    warning(paste("Could not match species:", sp, "site:", st))
  }
}

# Verify the array
cat("\n=== Detection History Array ===\n")
cat("Dimensions:", paste(dim(y), collapse = " x "), "\n")
cat("Number of NAs:", sum(is.na(y)), "\n")
cat("Number of 0s:", sum(y == 0), "\n")
cat("Number of 1s:", sum(y == 1), "\n")
cat("Proportion of 1s:", round(sum(y == 1) / length(y), 4), "\n")

# Show detection summary by species
cat("\n=== Detection Summary by Species ===\n")
for(i in 1:min(5, n_species)) {
  total_detections <- sum(y[i,,])
  sites_detected <- sum(apply(y[i,,], 1, sum) > 0)
  cat(sprintf("%-30s: %3d detections at %2d/%2d sites\n", 
              species_list[i], total_detections, sites_detected, n_sites))
}

# Extract site-level covariates (take first occurrence of each site)
site_covs <- bird_data_clean %>%
  group_by(SiteKey.x) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(SiteKey.x) %>%
  select(SiteKey.x, GEDI_biomass_filled, WSCI, Dist_Wat_25m_Uganda_250mDist,
         Precipitation_WorldClim, Bwindi_Heights_2mMaxarResample,
         Bwindi_Slope_2mMaxarResample, RH98_classified_2023,FHD_classified_2023,
         PAI_classified_2023,Cover_classified_2023,Population_Dens_Bwindi,
         TMIN_WorldClim,TAVG_WorldClim,TMAX_WorldClim,Bwindi_tsri)

# Check for NAs in covariates
cat("\n=== Covariate NAs ===\n")
cat("Biomass NAs:", sum(is.na(site_covs$GEDI_biomass_filled)), "\n")
cat("WSCI NAs:", sum(is.na(site_covs$WSCI)), "\n")
cat("Dist water NAs:", sum(is.na(site_covs$Dist_Wat_25m_Uganda_250mDist)), "\n")
cat("Precip NAs:", sum(is.na(site_covs$Precipitation_WorldClim)), "\n")
cat("Height NAs:", sum(is.na(site_covs$Bwindi_Heights_2mMaxarResample)), "\n")
cat("Slope NAs:", sum(is.na(site_covs$Bwindi_Slope_2mMaxarResample)), "\n")
cat("RH98 NAs:", sum(is.na(site_covs$RH98_classified_2023)), "\n")
cat("FHD NAs:", sum(is.na(site_covs$FHD_classified_2023)), "\n")
cat("PAI NAs:", sum(is.na(site_covs$PAI_classified_2023)), "\n")
cat("Cover NAs:", sum(is.na(site_covs$Cover_classified_2023)), "\n")
cat("Pop NAs:", sum(is.na(site_covs$Population_Dens_Bwindi)), "\n")
cat("Tmin NAs:", sum(is.na(site_covs$TMIN_WorldClim)), "\n")
cat("Tavg NAs:", sum(is.na(site_covs$TAVG_WorldClim)), "\n")
cat("Tmax NAs:", sum(is.na(site_covs$TMAX_WorldClim)), "\n")
cat("Tsri NAs:", sum(is.na(site_covs$Bwindi_tsri)), "\n")


# If there are any NAs in covariates, we need to handle them
# Option 1: Remove sites with NA covariates (not recommended)
# Option 2: Impute with mean (simple approach)
if(any(is.na(site_covs))) {
  cat("\nWarning: NAs found in covariates. Imputing with means.\n")
  for(col in names(site_covs)[-1]) {  # Skip SiteKey.x
    if(any(is.na(site_covs[[col]]))) {
      site_covs[[col]][is.na(site_covs[[col]])] <- mean(site_covs[[col]], na.rm = TRUE)
    }
  }
}

# Standardize covariates
biomass <- as.numeric(scale(site_covs$GEDI_biomass_filled))
wsci <- as.numeric(scale(site_covs$WSCI))
dist_water <- as.numeric(scale(site_covs$Dist_Wat_25m_Uganda_250mDist))
precip <- as.numeric(scale(site_covs$Precipitation_WorldClim))
height <- as.numeric(scale(site_covs$Bwindi_Heights_2mMaxarResample))
slope <- as.numeric(scale(site_covs$Bwindi_Slope_2mMaxarResample))
rh98<- as.numeric(scale(site_covs$RH98_classified_2023))
fhd<- as.numeric(scale(site_covs$FHD_classified_2023))
pai<- as.numeric(scale(site_covs$PAI_classified_2023))
cover<- as.numeric(scale(site_covs$Cover_classified_2023))
pop<- as.numeric(scale(site_covs$Population_Dens_Bwindi))
tmin<- as.numeric(scale(site_covs$TMIN_WorldClim))
tavg<- as.numeric(scale(site_covs$TAVG_WorldClim))
tmax<- as.numeric(scale(site_covs$TMAX_WorldClim))
tsri<- as.numeric(scale(site_covs$Bwindi_tsri))


# Verify covariates
cat("\n=== Standardized Covariates ===\n")
cat("Biomass - mean:", round(mean(biomass), 3), "sd:", round(sd(biomass), 3), "\n")
cat("WSCI - mean:", round(mean(wsci), 3), "sd:", round(sd(wsci), 3), "\n")
cat("Dist water - mean:", round(mean(dist_water), 3), "sd:", round(sd(dist_water), 3), "\n")
cat("Precip - mean:", round(mean(precip), 3), "sd:", round(sd(precip), 3), "\n")
cat("Height - mean:", round(mean(height), 3), "sd:", round(sd(height), 3), "\n")
cat("Slope - mean:", round(mean(slope), 3), "sd:", round(sd(slope), 3), "\n")
cat("RH98 - mean:", round(mean(rh98), 3), "sd:", round(sd(rh98), 3), "\n")
cat("FHD - mean:", round(mean(fhd), 3), "sd:", round(sd(fhd), 3), "\n")
cat("PAI - mean:", round(mean(pai), 3), "sd:", round(sd(pai), 3), "\n")
cat("Cover - mean:", round(mean(cover), 3), "sd:", round(sd(cover), 3), "\n")
cat("Pop - mean:", round(mean(pop), 3), "sd:", round(sd(pop), 3), "\n")
cat("Tmin - mean:", round(mean(tmin), 3), "sd:", round(sd(tmin), 3), "\n")
cat("Tavg - mean:", round(mean(tavg), 3), "sd:", round(sd(tavg), 3), "\n")
cat("Tmax - mean:", round(mean(tmax), 3), "sd:", round(sd(tmax), 3), "\n")
cat("Tsri - mean:", round(mean(tsri), 3), "sd:", round(sd(tsri), 3), "\n")


#########################################################################

# Create JAGS data list
jags_data <- list(
  y = y,
  n_species = n_species,
  n_sites = n_sites,
  n_visits = n_visits,
  biomass = biomass,
  wsci = wsci,
  dist_water = dist_water,
  precip = precip,
  height = height,
  slope = slope,
  rh98 = rh98,
  fhd = fhd,
  pai = pai,
  cover = cover,
  pop = pop,
  tmin = tmin,
  tavg = tavg,
  tmax = tmax,
  tsri = tsri
)

# Save species and site names for later reference
site_names <- sites
species_names <- species_list

write_rds(jags_data,"/Users/abarenb1/Documents/PhD/Chapter3/Bwindi_Birds_Variables_Jan2026.RDS")

#STEP 3: Load data and set up Occ Model
######################################################################################
library(dplyr) # for data organization
library(ggplot2) # for plotting

# Load array of community detection history 
# 261 x (sites) x 8 columns (surveys) x 22 species
# It is the same structure as a single species detection history, 
# but stacked with 17 additional species

community_detection_history <- read_rds("/Users/abarenb1/Documents/PhD/Chapter3/Bwindi_Birds_Variables_Jan2026.RDS")
# Structure
str(community_detection_history)

## ----msomCovariates--------------------------
library(jagsUI)

model_code <- "
model {
  
  # Priors for community-level hyperparameters
  mu.lpsi ~ dnorm(0, 0.1)
  tau.lpsi ~ dgamma(0.1, 0.1)
  sigma.lpsi <- 1/sqrt(tau.lpsi)
  
  mu.lp ~ dnorm(0, 0.1)
  tau.lp ~ dgamma(0.1, 0.1)
  sigma.lp <- 1/sqrt(tau.lp)
  
  # Priors for community-level covariate effects on occupancy
  mu.beta.biomass ~ dnorm(0, 0.1)
  tau.beta.biomass ~ dgamma(0.1, 0.1)
  
  mu.beta.wsci ~ dnorm(0, 0.1)
  tau.beta.wsci ~ dgamma(0.1, 0.1)
  
  mu.beta.precip ~ dnorm(0, 0.1)
  tau.beta.precip ~ dgamma(0.1, 0.1)
  
  mu.beta.slope ~ dnorm(0, 0.1)
  tau.beta.slope ~ dgamma(0.1, 0.1)
  
  mu.beta.dist_water ~ dnorm(0, 0.1)
  tau.beta.dist_water ~ dgamma(0.1, 0.1)
  
  mu.beta.height ~ dnorm(0, 0.1)
  tau.beta.height ~ dgamma(0.1, 0.1)
  
  # NEW COVARIATES
  mu.beta.rh98 ~ dnorm(0, 0.1)
  tau.beta.rh98 ~ dgamma(0.1, 0.1)
  
  mu.beta.fhd ~ dnorm(0, 0.1)
  tau.beta.fhd ~ dgamma(0.1, 0.1)
  
  mu.beta.pai ~ dnorm(0, 0.1)
  tau.beta.pai ~ dgamma(0.1, 0.1)
  
  mu.beta.cover ~ dnorm(0, 0.1)
  tau.beta.cover ~ dgamma(0.1, 0.1)
  
  mu.beta.pop ~ dnorm(0, 0.1)
  tau.beta.pop ~ dgamma(0.1, 0.1)
  
  mu.beta.tmin ~ dnorm(0, 0.1)
  tau.beta.tmin ~ dgamma(0.1, 0.1)
  
  mu.beta.tavg ~ dnorm(0, 0.1)
  tau.beta.tavg ~ dgamma(0.1, 0.1)
  
  mu.beta.tmax ~ dnorm(0, 0.1)
  tau.beta.tmax ~ dgamma(0.1, 0.1)
  
  mu.beta.tsri ~ dnorm(0, 0.1)
  tau.beta.tsri ~ dgamma(0.1, 0.1)
  
  # Species-specific parameters
  for(i in 1:n_species) {
    
    lpsi[i] ~ dnorm(mu.lpsi, tau.lpsi)
    lp[i] ~ dnorm(mu.lp, tau.lp)
    
    # Species-specific covariate effects
    beta.biomass[i] ~ dnorm(mu.beta.biomass, tau.beta.biomass)
    beta.wsci[i] ~ dnorm(mu.beta.wsci, tau.beta.wsci)
    beta.precip[i] ~ dnorm(mu.beta.precip, tau.beta.precip)
    beta.slope[i] ~ dnorm(mu.beta.slope, tau.beta.slope)
    beta.dist_water[i] ~ dnorm(mu.beta.dist_water, tau.beta.dist_water)
    beta.height[i] ~ dnorm(mu.beta.height, tau.beta.height)
    beta.rh98[i] ~ dnorm(mu.beta.rh98, tau.beta.rh98)
    beta.fhd[i] ~ dnorm(mu.beta.fhd, tau.beta.fhd)
    beta.pai[i] ~ dnorm(mu.beta.pai, tau.beta.pai)
    beta.cover[i] ~ dnorm(mu.beta.cover, tau.beta.cover)
    beta.pop[i] ~ dnorm(mu.beta.pop, tau.beta.pop)
    beta.tmin[i] ~ dnorm(mu.beta.tmin, tau.beta.tmin)
    beta.tavg[i] ~ dnorm(mu.beta.tavg, tau.beta.tavg)
    beta.tmax[i] ~ dnorm(mu.beta.tmax, tau.beta.tmax)
    beta.tsri[i] ~ dnorm(mu.beta.tsri, tau.beta.tsri)
    
    # Ecological model (occupancy)
    for(j in 1:n_sites) {
      
      # Occupancy probability as function of covariates
      logit(psi[i,j]) <- lpsi[i] + 
                         beta.biomass[i] * biomass[j] + 
                          beta.wsci[i] * wsci[j] + 
                         beta.precip[i] * precip[j] +
                         beta.slope[i] * slope[j] +
                         beta.dist_water[i] * dist_water[j] +
                         beta.height[i] * height[j] +
                         beta.rh98[i] * rh98[j] +
                         beta.fhd[i] * fhd[j] +
                         beta.pai[i] * pai[j] +
                         beta.cover[i] * cover[j] +
                         beta.pop[i] * pop[j] +
                         beta.tmin[i] * tmin[j] +
                         beta.tavg[i] * tavg[j] +
                         beta.tmax[i] * tmax[j] +
                         beta.tsri[i] * tsri[j]
      
      z[i,j] ~ dbern(psi[i,j])
      
      for(k in 1:n_visits) {
        logit(p[i,j,k]) <- lp[i]
        mu.y[i,j,k] <- z[i,j] * p[i,j,k]
        y[i,j,k] ~ dbern(mu.y[i,j,k])
      }
    }
    
    psi.fs[i] <- sum(z[i,])/n_sites
  }
  
  mean.psi <- mean(psi.fs[])
  
  for(j in 1:n_sites) {
    richness[j] <- sum(z[,j])
  }
  
  total.richness <- sum(psi.fs[])
}
"

writeLines(model_code, "multispecies_occupancy.txt")

# ============================================================
# STEP 3: Initial Values (UPDATED)
# ============================================================

inits <- function() {
  z_init <- apply(y, c(1,2), function(x) {
    if(sum(x) > 0) return(1) else return(0)
  })
  
  list(
    z = z_init,
    mu.lpsi = rnorm(1, 0, 0.5),
    mu.lp = rnorm(1, 0, 0.5),
    tau.lpsi = rgamma(1, 0.1, 0.1),
    tau.lp = rgamma(1, 0.1, 0.1),
    mu.beta.biomass = rnorm(1, 0, 0.5),
    mu.beta.wsci = rnorm(1, 0, 0.5),
    mu.beta.precip = rnorm(1, 0, 0.5),
    mu.beta.slope = rnorm(1, 0, 0.5),
    mu.beta.dist_water = rnorm(1, 0, 0.5),
    mu.beta.height = rnorm(1, 0, 0.5),
    mu.beta.rh98 = rnorm(1, 0, 0.5),
    mu.beta.fhd = rnorm(1, 0, 0.5),
    mu.beta.pai = rnorm(1, 0, 0.5),
    mu.beta.cover = rnorm(1, 0, 0.5),
    mu.beta.pop = rnorm(1, 0, 0.5),
    mu.beta.tmin = rnorm(1, 0, 0.5),
    mu.beta.tavg = rnorm(1, 0, 0.5),
    mu.beta.tmax = rnorm(1, 0, 0.5),
    mu.beta.tsri = rnorm(1, 0, 0.5),
    tau.beta.biomass = rgamma(1, 0.1, 0.1),
    tau.beta.wsci = rgamma(1, 0.1, 0.1),
    tau.beta.precip = rgamma(1, 0.1, 0.1),
    tau.beta.slope = rgamma(1, 0.1, 0.1),
    tau.beta.dist_water = rgamma(1, 0.1, 0.1),
    tau.beta.height = rgamma(1, 0.1, 0.1),
    tau.beta.rh98 = rgamma(1, 0.1, 0.1),
    tau.beta.fhd = rgamma(1, 0.1, 0.1),
    tau.beta.pai = rgamma(1, 0.1, 0.1),
    tau.beta.cover = rgamma(1, 0.1, 0.1),
    tau.beta.pop = rgamma(1, 0.1, 0.1),
    tau.beta.tmin = rgamma(1, 0.1, 0.1),
    tau.beta.tavg = rgamma(1, 0.1, 0.1),
    tau.beta.tmax = rgamma(1, 0.1, 0.1),
    tau.beta.tsri = rgamma(1, 0.1, 0.1)
  )
}

# ============================================================
# STEP 4: Parameters to Monitor (UPDATED)
# ============================================================

params <- c(
  # Community-level parameters
  "mu.lpsi", "sigma.lpsi",
  "mu.lp", "sigma.lp",
  "mu.beta.biomass","mu.beta.wsci", "mu.beta.precip", "mu.beta.slope", 
  "mu.beta.dist_water", "mu.beta.height",
  "mu.beta.rh98", "mu.beta.fhd", "mu.beta.pai", "mu.beta.cover",
  "mu.beta.pop", "mu.beta.tmin", "mu.beta.tavg", "mu.beta.tmax", "mu.beta.tsri",
  
  # Species-specific parameters
  "lpsi", "lp",
  "beta.biomass", "beta.wsci", "beta.precip", "beta.slope", "beta.dist_water", "beta.height",
  "beta.rh98", "beta.fhd", "beta.pai", "beta.cover", "beta.pop",
  "beta.tmin", "beta.tavg", "beta.tmax", "beta.tsri",
  
  # Derived parameters
  "psi.fs", "mean.psi", "richness", "total.richness", "z"
)

# ============================================================
# STEP 5: Run JAGS Model
# ============================================================

n.adapt <- 1000
n.iter <- 10000
n.burnin <- 5000
n.thin <- 5
n.chains <- 3

jags_out <- jagsUI::jags(
  data = jags_data,
  inits = inits,
  parameters.to.save = params,
  model.file = "multispecies_occupancy.txt",
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.iter = n.iter,
  n.burnin = n.burnin,
  n.thin = n.thin,
  parallel = TRUE
)

# ============================================================
# STEP 6: Examine Results
# ============================================================

print(jags_out)

max(jags_out$summary[, "Rhat"], na.rm = TRUE)

# Community-level results (ALL COVARIATES)
jags_out$summary[c("mu.lpsi", "sigma.lpsi", "mu.lp", "sigma.lp",
                   "mu.beta.biomass","mu.beta.wsci", "mu.beta.precip", "mu.beta.slope",
                   "mu.beta.dist_water", "mu.beta.height",
                   "mu.beta.rh98", "mu.beta.fhd", "mu.beta.pai", "mu.beta.cover",
                   "mu.beta.pop", "mu.beta.tmin", "mu.beta.tavg", "mu.beta.tmax", "mu.beta.tsri"), ]

# Species-specific occupancy (ALL COVARIATES)
species_results <- data.frame(
  species = species_names,
  occupancy_logit = jags_out$mean$lpsi,
  detection_logit = jags_out$mean$lp,
  beta_biomass = jags_out$mean$beta.biomass,
  beta_wsci = jags_out$mean$beta.wsci,
  beta_precip = jags_out$mean$beta.precip,
  beta_slope = jags_out$mean$beta.slope,
  beta_dist_water = jags_out$mean$beta.dist_water,
  beta_height = jags_out$mean$beta.height,
  beta_rh98 = jags_out$mean$beta.rh98,
  beta_fhd = jags_out$mean$beta.fhd,
  beta_pai = jags_out$mean$beta.pai,
  beta_cover = jags_out$mean$beta.cover,
  beta_pop = jags_out$mean$beta.pop,
  beta_tmin = jags_out$mean$beta.tmin,
  beta_tavg = jags_out$mean$beta.tavg,
  beta_tmax = jags_out$mean$beta.tmax,
  beta_tsri = jags_out$mean$beta.tsri,
  prop_occupied = jags_out$mean$psi.fs
)

print(species_results)

# ============================================================
# STEP 7: Visualize Results (ALL COVARIATES)
# ============================================================

library(ggplot2)

covariate_effects <- data.frame(
  covariate = rep(c("Biomass","WSCI", "Precipitation", "Slope", "DistWater", "Height",
                    "RH98", "FHD", "PAI", "Cover", "Pop", 
                    "Tmin", "Tavg", "Tmax", "TSRI"), 
                  each = length(jags_out$sims.list$mu.beta.biomass)),
  effect = c(jags_out$sims.list$mu.beta.biomass,
             jags_out$sims.list$mu.beta.wsci,
             jags_out$sims.list$mu.beta.precip,
             jags_out$sims.list$mu.beta.slope,
             jags_out$sims.list$mu.beta.dist_water,
             jags_out$sims.list$mu.beta.height,
             jags_out$sims.list$mu.beta.rh98,
             jags_out$sims.list$mu.beta.fhd,
             jags_out$sims.list$mu.beta.pai,
             jags_out$sims.list$mu.beta.cover,
             jags_out$sims.list$mu.beta.pop,
             jags_out$sims.list$mu.beta.tmin,
             jags_out$sims.list$mu.beta.tavg,
             jags_out$sims.list$mu.beta.tmax,
             jags_out$sims.list$mu.beta.tsri)
)

# Extract summary statistics for all covariate effects
covariate_names <- c("Biomass","WSCI", "Precipitation", "Slope", "Dist_Water", "Height",
                     "RH98", "FHD", "PAI", "Cover", "Pop", 
                     "Tmin", "Tavg", "Tmax", "TSRI")

covariate_params <- c("mu.beta.biomass","mu.beta.wsci", "mu.beta.precip", "mu.beta.slope", 
                      "mu.beta.dist_water", "mu.beta.height",
                      "mu.beta.rh98", "mu.beta.fhd", "mu.beta.pai", 
                      "mu.beta.cover", "mu.beta.pop",
                      "mu.beta.tmin", "mu.beta.tavg", "mu.beta.tmax", "mu.beta.tsri")

coef_summary <- data.frame(
  covariate = covariate_names,
  mean = jags_out$summary[covariate_params, "mean"],
  lower = jags_out$summary[covariate_params, "2.5%"],
  upper = jags_out$summary[covariate_params, "97.5%"],
  significant = ifelse(jags_out$summary[covariate_params, "2.5%"] > 0 | 
                         jags_out$summary[covariate_params, "97.5%"] < 0, "Yes", "No")
)

# Order by effect size
coef_summary <- coef_summary %>%
  arrange(mean) %>%
  mutate(covariate = factor(covariate, levels = covariate))

#Community effect

library(forestplot)

# Prepare data for forest plot
forest_data <- coef_summary %>%
  mutate(
    estimate = sprintf("%.3f", mean),
    ci = sprintf("(%.3f, %.3f)", lower, upper),
    text = paste0(covariate, "\n", estimate, " ", ci)
  )

# Create forest plot using ggplot
ggplot(forest_data, aes(y = reorder(covariate, mean))) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1, color = "gray40") +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = significant), 
                 height = 0.3, size = 1.2) +
  geom_point(aes(x = mean, color = significant), size = 4, shape = 18) +
  scale_color_manual(values = c("No" = "gray50", "Yes" = "#D55E00"),
                     labels = c("Not significant", "Significant")) +
  theme_classic(base_size = 14) +
  labs(title = "Community-level Covariate Effects",
       x = "Effect Size (95% Credible Interval)",
       y = NULL,
       color = NULL) +
  theme(
    legend.position = "top",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linetype = "dotted")
  )

ggsave("covariate_effects_forest_plot.png", width = 10, height = 6, dpi = 300)

saveRDS(jags_out, "multispecies_occupancy_results_full.rds")
write.csv(species_results, "species_occupancy_results_full.csv", row.names = FALSE)

# Plot species richness across sites
richness_summary <- data.frame(
  site = 1:n_sites,
  mean_richness = jags_out$mean$richness,
  sd_richness = jags_out$sd$richness
)

ggplot(richness_summary, aes(x = site, y = mean_richness)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_richness - sd_richness,
                    ymax = mean_richness + sd_richness),
                width = 0.2) +
  theme_minimal() +
  labs(title = "Estimated Species Richness by Site",
       x = "Site", y = "Species Richness")

# Plot species-specific responses to biomass
species_height <- data.frame(
  species = species_list,
  effect = jags_out$mean$beta.height,
  lower = jags_out$q2.5$beta.height,
  upper = jags_out$q97.5$beta.height
)

ggplot(species_height, aes(x = reorder(species, effect), y = effect)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_minimal() +
  labs(x = "Species", y = "Effect of Maxar Height on Occupancy")+
  theme(text = element_text(size = 20))      

# Save results
saveRDS(jags_out, "multispecies_occupancy_results.rds")
write.csv(species_results, "species_occupancy_results.csv", row.names = FALSE)

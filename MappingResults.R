library(terra)
library(viridis)
library(sf)
library(sp)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(patchwork)
library(tidyr)
library(corrplot)


#Resample Variables
#Build Stack of Covariate Rasters
matching_tifs <- c("GEDI_biomass_filled","Dist_Wat_25m_Uganda_250mDist","Precipitation_WorldClim",
                   "WSCI","Bwindi_Heights_2mMaxarResample","Bwindi_2mFilledDSM_MaxarResample",
                   "Bwindi_Slope_2mMaxarResample","FHD_classified_2023",
                   "ETH_GlobalCanopyHeight_2020_10m_v1","ETH_GlobalCanopyHeightSD_2020_10m_v1","Population_Dens_Bwindi",
                   "TMIN_WorldClim","TAVG_WorldClim","TMAX_WorldClim","Bwindi_tsri","Bwindi_tpi","DistanceToBorder",
                   "RH98_classified_2023_Mar12_2026","LANDSAT_08_SR_07-01-2022-01-30-2023_NDVI2")

# Load all rasters first
rasterList <- lapply(matching_tifs, function(x) 
  rast(paste0("/Users/abarenb1/Documents/PhD/Chapter3/Variables/", x, ".tif")))

# Set the first raster as template (or choose a specific one)
template <- rasterList[[1]]

# Align all rasters to template (same CRS, extent, resolution)
rasterList_aligned <- lapply(1:length(rasterList), function(i) {
  r <- rasterList[[i]]
  
  # Check if geometry matches
  if (!compareGeom(r, template, stopOnError = FALSE)) {
    print(paste("Aligning", matching_tifs[i]))
    
    # Determine method based on data type
    # Use "near" for classified layers, "bilinear" for continuous
    method <- ifelse(grepl("classified", matching_tifs[i]), "near", "bilinear")
    
    # Project and resample to template
    project(r, template, method = method)
  } else {
    r
  }
})

#Save new rasters
writeRaster(rasterList_aligned[[18]], "/Users/abarenb1/Documents/PhD/Chapter3/VariablesResample/RH98_classified_2023_Mar12_2026.tif", overwrite = TRUE)

##Check Variable Correlation
######################################################################
matching_tifs <- c("GEDI_biomass_filled","Dist_Wat_25m_Uganda_250mDist","Precipitation_WorldClim",
                   "WSCI","Bwindi_Heights_2mMaxarResample","Bwindi_2mFilledDSM_MaxarResample",
                   "Bwindi_Slope_2mMaxarResample","ETH_GlobalCanopyHeightSD_2020_10m_v1","Population_Dens_Bwindi",
                   "TMIN_WorldClim","TAVG_WorldClim","TMAX_WorldClim","Bwindi_tsri","DistanceToBorder",
                   "RH98_classified_2023_Mar12_2026","LANDSAT_08_SR_07-01-2022-01-30-2023_NDVI")

# Load all rasters first
rasterList <- lapply(matching_tifs, function(x) 
  rast(paste0("/Users/abarenb1/Documents/PhD/Chapter3/VariablesResample/", x, ".tif")))

birds<-st_read("/Users/abarenb1/Documents/PhD/Chapter3/Bwindi Endemic and threatened birds/Species_with_Coordinates_ESPG4326.shp")

for (j in 1:length(matching_tifs)){
  print(paste("Extracting from:", matching_tifs[j]))
  
  ras <- rasterList[[j]]
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


birds<-st_drop_geometry(birds)

data_subset <- birds[c(27:34,36:42)]
# Remove columns with zero variance
data_subset <- data_subset[, sapply(data_subset, function(x) sd(x, na.rm = TRUE) > 0)]

colnames(data_subset)<-c("Biomass","Distance Water","Precipitation","WSCI","Canopy Height", "DSM",
                         "Slope","SD Height","Temp Min","Temp Avg","Temp Max", "TSRI", "Distance Border",
                         "RH98","NDVI")

# Check how many columns remain
cat("Columns after removing constant variables:", ncol(data_subset), "\n")

# Now calculate correlation
x <- cor(data_subset, use = "complete.obs")

colnames(x) <- substr(colnames(x), 1, 15)
rownames(x) <- substr(rownames(x), 1, 15)

corrplot(x, type="upper", order="hclust")


##Mapping species results
########################################################################
# Define models and their paths
model_info <- list(
  LAND = "/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-03-26.SSDM/Species",
  STRU  = "/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Structure-03-13-26.SSDM/Species",
  LATO  = "/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-Topography-03-26.SSDM/Species",
  LATOE = "/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-Topography-Elevation-03-13-26.SSDM/Species",
  LTSR  = "/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-Topography-Structure-03-15-26.SSDM/Species",
  LTSRE = "/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-Topography-Structure-Elevation.SSDM/Species"
)

# Function to load and name rasters for a given model
load_model_rasters <- function(path, model_name) {
  
  rastlist      <- list.files(path = path, pattern = '.tif$', all.files = TRUE, full.names = TRUE)
  rastlistshort <- list.files(path = path, pattern = '.tif$', all.files = TRUE, full.names = FALSE)
  
  allrasters <- lapply(rastlist, rast)
  names(allrasters) <- tools::file_path_sans_ext(rastlistshort)
  
  # Remove model-specific probability prefix — adjust pattern per model if needed
  names(allrasters) <- gsub("ProbabilityGEDIModel", "", names(allrasters))
  
  return(allrasters)
}

# Load all models
allrastersLAND  <- load_model_rasters(model_info$LAND,  "LAND")  
allrastersSTRU  <- load_model_rasters(model_info$STRU,  "STRU")    # STRU uses "ProbabilityGEDIModel"
allrastersLATO  <- load_model_rasters(model_info$LATO,  "LATO")
allrastersLATOE <- load_model_rasters(model_info$LATOE, "LATOE")
allrastersLTSR  <- load_model_rasters(model_info$LTSR,  "LTSR")
allrastersLTSRE <- load_model_rasters(model_info$LTSRE, "LTSRE")

redFaced <- list(
  #STRU  = allrastersSTRU[[2]],
  #LAND  = mask(allrastersLAND[[10]], allrastersLATO[[13]]),
  LATO  = allrastersLATO[[13]],
  LATOE = allrastersLATOE[[13]],
  LTSR  = allrastersLTSR[[13]],
  LTSRE = allrastersLTSRE[[13]]
)

Invert <-allrasters[c(1,2,5,7,8,10,11,13,14,16:20)]
Nectar <-allrasters[c(3,12,15)]
Grani <-allrasters[c(4,9)]
Omni <-allrasters[c(6)]

#Ordered <- allrastersS[c(1,2,4,5,7,8,10,11,13,14,15,16,17,9,12,6,3)] #Land
#Ordered <- allrastersS[c(1,2,3,5,7,8,11,12,13,14,15,9,10,6,4)] #Struct
Ordered <- allrastersS[c(1,2,5,7,8,10,11,13,14,16:20,3,12,15,4,9,6)] #LandTopoStruct

#FOR SPECIFIC LIST
plot_rasters_facet <- function(raster_list) {
  
  df_combined <- bind_rows(lapply(seq_along(raster_list), function(i) {
    r <- raster_list[[i]]
    
    # If multi-layer, collapse to single layer by selecting the first layer
    if (nlyr(r) > 1) {
      r <- r[[1]]
    }
    
    df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
    
    # Rename value column
    colnames(df)[ncol(df)] <- "value"
    
    # Assign name label
    df$name <- if (!is.null(names(raster_list)[i]) && names(raster_list)[i] != "") {
      gsub(" ", "\n", names(raster_list)[i])
    } else {
      paste0("Raster_", i)
    }
    
    return(df)
  }))
  
  # Remove NAs
  df_combined <- df_combined %>%
    filter(!is.na(value)) %>%
    drop_na()
  
  ggplot(df_combined, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Presence \n Probability") +
    facet_wrap(~name, ncol = 2) +
    coord_equal() +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    )
}

plot <- plot_rasters_facet(redFaced)
plot

#FOR FULL LIST
plot_rasters_facet <- function(raster_list) {
  # Combine all rasters into one data frame
  df_combined <- bind_rows(lapply(seq_along(raster_list), function(i) {
    df <- as.data.frame(raster_list[[i]], xy = TRUE, na.rm = FALSE)
    
    # Rename the last column to "value"
    colnames(df)[ncol(df)] <- "value"
    
    df$name <- if(!is.null(names(raster_list)[i]) && names(raster_list)[i] != "") {
      # Replace spaces with newlines for wrapping at every word
      gsub(" ", "\n", names(raster_list)[i])
    } else {
      paste0("Raster_", i)
    }
    
    return(df)
  }))
  
  # Remove NAs
  df_combined <- df_combined %>% 
    filter(!is.na(value)) %>%
    drop_na()
  
  # Get the order of names as they appear
  name_order <- unique(df_combined$name)
  
  # Manually insert dummy placeholders to force line breaks
  # Break after position 14, 17, and 19
  #ALL
  final_order <- c(
    name_order[1:14],           # First 14 rasters (positions 1-14)
    "dummy1",                    # Force raster 15 to new line
    name_order[15:17],           # Rasters 15-17
    "dummy3", "dummy4",  # Force raster 18 to new line
    name_order[18:19],           # Rasters 18-19
    "dummy6","dummy7","dummy8", 
    name_order[20],"dummy10","dummy11","dummy12","dummy13"# Force raster 20 to new line
  )
  
  #final_order <- c(
  #  name_order[1:11],           # First 14 rasters (positions 1-14)
  #  "dummy1","dummy2","dummy3", "dummy4",                    # Force raster 15 to new line
  #  name_order[12:13],           # Rasters 15-17
  #  "dummy5","dummy6","dummy7",  # Force raster 18 to new line
  #  name_order[14],           # Rasters 18-19
  #  "dummy8", "dummy9","dummy10","dummy11",
  #  name_order[15],"dummy12","dummy13","dummy14","dummy15"# Force raster 20 to new line
  #)
  
  #LAND
  #final_order <- c(
  #  name_order[1:13],           # First 14 rasters (positions 1-14)
  #  "dummy1","dummy2",                    # Force raster 15 to new line
  #  name_order[14:15],           # Rasters 15-17
  #  "dummy3", "dummy4","dummy5",  # Force raster 18 to new line
  #  name_order[16],           # Rasters 18-19
  #  "dummy6","dummy7","dummy8", "dummy9",
  #  name_order[17],"dummy10","dummy11","dummy12","dummy13"# Force raster 20 to new line
  #)
  
  # Add remaining rasters if they exist
  if(length(name_order) >= 20) {
    final_order <- c(final_order, name_order[20:length(name_order)])
  }
  final_order <-final_order[1:30]
  
  # Convert name to factor
  df_combined$name <- factor(df_combined$name, levels = final_order)
  
  # Create the plot
  p <- ggplot(df_combined, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(name="Presence \n Probability") +
    facet_wrap(~name, ncol = 5, drop = FALSE) +
    coord_equal() +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    )
  
  return(p)
}

plot <- plot_rasters_facet(Ordered)
plot

ggsave("/Users/abarenb1/Documents/PhD/Chapter3/PaperFigures/SDMsMapped_Landscape-Topography-Structure-Elevation-03-16-26.png", 
       width = 30, height = 50, units = "cm")

#Charting Species Detections
############################################

data <- read.csv("/Users/abarenb1/Documents/PhD/Chapter3/Sightings.csv")

ordered_plot <- ggplot(data, aes(x = reorder(Common.Name, -Sites), y = Sites)) +
  geom_bar(stat = "identity", fill = "purple") +
  xlab("Species") +
  ylab("Sites Detected At (n=261)") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text( 
    angle = 60, 
    size=12,vjust = 1, hjust=1)
  )+
theme(axis.text.y = element_text( 
  size=12)
)+theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))+
  geom_text(aes(label = Label), nudge_y = 3)+
  annotate("text", x=19, y=110, label= "* = Near Threatened")+
  annotate("text", x=18.5, y=100, label= "** = Vulnerable")
 # geom_vline(xintercept=c(15.5), linetype='dashed', color=c('red'))
ordered_plot

#Charting Variable Importance
############################################
dataLAND <- read.csv("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-03-26.SSDM/Stack/Tables/VarImp.csv")
dataSTRU <- read.csv("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Structure-03-13-26.SSDM/Stack/Tables/VarImp.csv")
dataLATO <- read.csv("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-Topography-03-26.SSDM/Stack/Tables/VarImp.csv")
dataLATOE <- read.csv("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-Topography-Elevation-03-13-26.SSDM/Stack/Tables/VarImp.csv")
dataLTSR <- read.csv("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-Topography-Structure-03-15-26.SSDM/Stack/Tables/VarImp.csv")
dataLTSRE <- read.csv("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-Topography-Structure-Elevation.SSDM/Stack/Tables/VarImp.csv")


dataInvert <- read.csv("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_FEB/Invertivore-03-02-26.SSDM/Stack/Tables/VarImp.csv")
dataNectar <- read.csv("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_FEB/Nectar.SSDM/Stack/Tables/VarImp.csv")
dataGrani <- read.csv("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_FEB/Granivore-03-02-26.SSDM/Stack/Tables/VarImp.csv")
dataOmni <-read.csv("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_FEB/Omnivore-03-03-26.Ensemble.SDM/Tables/VarImp.csv")

colnames(dataLAND)<-c("X","Distance Water","NDVI","Precipitation")
colnames(dataSTRU)<-c("X","Canopy Height","SD Height","Biomass","RH98","WSCI")
colnames(dataLATO)<-c("X","Slope","TSRI","Distance Water","NDVI","Precipitation")
colnames(dataLATOE)<-c("X","DSM","Slope","TSRI","Distance Water","NDVI","Precipitation")
colnames(dataLTSR)<-c("X","Canopy Height","Slope","TSRI","Distance Water","SD Height","Biomass","NDVI","Precipitation","RH98","WSCI")
colnames(dataLTSRE)<-c("X","DSM","Canopy Height","Slope","TSRI","Distance Water","SD Height","Biomass","NDVI","Precipitation","RH98","WSCI")


dataLAND <- dataLAND %>% pivot_longer(cols=c("Distance Water","NDVI","Precipitation"),
                            names_to='Variable',
                            values_to=c('Value'))

dataLAND_wide <- dataLAND %>%
  pivot_wider(names_from = X, values_from = Value)

dataSTRU <- dataSTRU %>% pivot_longer(cols=c("Canopy Height","SD Height","Biomass","RH98","WSCI"),
                                      names_to='Variable',
                                      values_to=c('Value'))

dataSTRU_wide <- dataSTRU %>%
  pivot_wider(names_from = X, values_from = Value)

dataLATO <- dataLATO %>% pivot_longer(cols=c("Slope","TSRI","Distance Water","NDVI","Precipitation"),
                                      names_to='Variable',
                                      values_to=c('Value'))

dataLATO_wide <- dataLATO %>%
  pivot_wider(names_from = X, values_from = Value)

dataLATOE <- dataLATOE %>% pivot_longer(cols=c("DSM","Slope","TSRI","Distance Water","NDVI","Precipitation"),
                                      names_to='Variable',
                                      values_to=c('Value'))

dataLATOE_wide <- dataLATOE %>%
  pivot_wider(names_from = X, values_from = Value)

dataLTSR <- dataLTSR %>% pivot_longer(cols=c("Canopy Height","Slope","TSRI","Distance Water","SD Height","Biomass","NDVI","Precipitation","RH98","WSCI"),
                                        names_to='Variable',
                                        values_to=c('Value'))

dataLTSR_wide <- dataLTSR %>%
  pivot_wider(names_from = X, values_from = Value)

dataLTSRE <- dataLTSRE %>% pivot_longer(cols=c("DSM","Canopy Height","Slope","TSRI","Distance Water","SD Height","Biomass","NDVI","Precipitation","RH98","WSCI"),
                                      names_to='Variable',
                                      values_to=c('Value'))

dataLTSRE_wide <- dataLTSRE %>%
  pivot_wider(names_from = X, values_from = Value)

##START PLOTTING
#####################################################################
  
ggplot(dataLAND_wide , aes(x=reorder(Variable, Mean), y=Mean)) + 
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.3) +
    geom_point(size=2)+
    labs(y = "Variable Importance ",
         x = NULL,
         color = NULL)+theme(text = element_text(size = 18))+coord_flip()
  

ggplot(dataSTRU_wide , aes(x=reorder(Variable, Mean), y=Mean)) + 
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.3) +
    geom_point(size=2)+
  labs(y = "Variable Importance ",
       x = NULL,
       color = NULL)+theme(text = element_text(size = 18))+coord_flip()

ggplot(dataLATO_wide , aes(x=reorder(Variable, Mean), y=Mean)) + 
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.3) +
  geom_point(size=2)+
  labs(y = "Variable Importance ",
       x = NULL,
       color = NULL)+theme(text = element_text(size = 18))+coord_flip()

ggplot(dataLATOE_wide , aes(x=reorder(Variable, Mean), y=Mean)) + 
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.3) +
  geom_point(size=2)+
  labs(y = "Variable Importance ",
       x = NULL,
       color = NULL)+theme(text = element_text(size = 15))+coord_flip()

ggplot(dataLTSR_wide , aes(x=reorder(Variable, Mean), y=Mean)) + 
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.3) +
  geom_point(size=2)+
  labs(y = "Variable Importance ",
       x = NULL,
       color = NULL)+theme(text = element_text(size = 15))+coord_flip()

ggplot(dataLTSRE_wide , aes(x=reorder(Variable, Mean), y=Mean)) + 
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.3) +
  geom_point(size=2)+
  labs(y = "Variable Importance ",
       x = NULL,
       color = NULL)+theme(text = element_text(size = 15))+coord_flip()
  
ggplot(dataLATOE_wide , aes(x=reorder(Variable, Axes.evaluation), y=Axes.evaluation)) + 
    geom_point(size=2)+
    labs(y = "Variable Importance ",
         x = NULL,
         color = NULL)+theme(text = element_text(size = 20))+coord_flip()+
  scale_x_discrete(label = labelsO)

##Variable Importance for Each Species
################################################################
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(reshape2)

# Set working directory to where your CSV files are located
setwd("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-Topography-Structure-Elevation.SSDM/Species")

# Get all CSV files matching your pattern
csv_files <- list.files(pattern = "GEDIModel.*\\.csv", full.names = TRUE)

# Read all CSV files and combine into one dataframe
all_species_data <- list()

# Get all CSV files matching your pattern
csv_files <- list.files(pattern = "GEDIModel.*\\.csv", full.names = TRUE)

# Read all CSV files and combine into one dataframe
all_species_data <- list()

for (file in csv_files) {
  # Extract species name from filename
  species_name <- gsub("GEDIModel", "", basename(file))
  species_name <- gsub("\\.csv", "", species_name)
  
  # Read the CSV file
  df <- read.csv(file, row.names = 1)
  
  # Extract the first row (variable importance values)
  # Transpose to get variables as rows
  importance_values <- as.numeric(df[1, ])
  names(importance_values) <- colnames(df)
  
  # Store in list with species name
  all_species_data[[species_name]] <- importance_values
}

# Combine into a dataframe
combined_df <- do.call(cbind, all_species_data)
combined_df <- as.data.frame(combined_df)

# Display the comparison table
print("Variable Importance Comparison Across Species:")
print(combined_df)

combined_long <- combined_df %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Species", values_to = "Importance")

combined_long <- combined_long %>%
  mutate(Species_Wrapped = str_replace_all(Species, " ", "\n"))

variable_labels <- c(
  "WSCI" = "WSCI",
  "ETH_GlobalCanopyHeightSD_2020_10m_v1" = "SD Height",
  "RH98_classified_2023_Mar12_2026" = "RH98", 
  "GEDI_biomass_filled" = "Biomass",
  "Bwindi_Heights_2mMaxarResample" = "Height",
  "Bwindi_2mFilledDSM_MaxarResample" = "DSM",
  "Bwindi_Slope_2mMaxarResample" = "Slope",
  "Bwindi_tsri" = "TSRI",
  "Dist_Wat_25m_Uganda_250mDist" = "Water",
  "Precipitation_WorldClim" = "PPTN",
  "LANDSAT_08_SR_07.01.2022.01.30.2023_NDVI" = "NDVI"
)

# Add short labels to the long dataframe
combined_long <- combined_long %>%
  mutate(Variable_Label = variable_labels[Variable])

# Identify top variable for each species
top_vars <- combined_df %>%
  summarise(across(everything(), ~ rownames(combined_df)[which.max(.)]))

# Create a dataframe marking which observations are top variables
combined_long <- combined_long %>%
  group_by(Species) %>%
  mutate(Is_Top = Variable == Variable[which.max(Importance)]) %>%
  ungroup()



group3_vars <- c("GEDI_biomass_filled", "Bwindi_Heights_2mMaxarResample", 
                 "RH98_classified_2023_Mar12_2026", 
                 "ETH_GlobalCanopyHeightSD_2020_10m_v1", "WSCI")

group2_vars <- c("Bwindi_2mFilledDSM_MaxarResample", 
                 "Bwindi_Slope_2mMaxarResample","Bwindi_tsri")

group1_vars <- c("Precipitation_WorldClim","Dist_Wat_25m_Uganda_250mDist","LANDSAT_08_SR_07.01.2022.01.30.2023_NDVI")

group3_label <- "Structure"
group2_label <- "Topography"
group1_label <- "Landscape"

# Create grouped data
# Create grouped data
combined_grouped <- combined_long %>%
  mutate(Variable_Group = ifelse(Variable %in% group1_vars, 
                                 group1_label, 
                                 ifelse(Variable %in% group2_vars,
                                        group2_label,
                                        group3_label))) %>%
  group_by(Species, Species_Wrapped, Variable_Group) 

# Set consistent color palette for all variables
n_vars <- length(unique(combined_long$Variable))
color_palette <- setNames(
  rainbow(n_vars, s = 0.7, v = 0.8),
  unique(combined_long$Variable)
)

# Calculate global y-axis limits
y_max <- max(combined_long$Importance) * 1.05  # Add 5% padding


p3 <- p3 <- ggplot(combined_grouped, aes(fill=Variable_Group,
                                         x = reorder(Variable_Label, as.numeric(factor(Variable_Group))), 
                                         y = Importance,
                                         group=factor(Variable_Group,ordered=TRUE))) +
  geom_bar(stat = "identity", 
           aes(alpha = Is_Top),
           color = "black", 
           size = 0.3) +
  geom_text(data = filter(combined_grouped, Is_Top == TRUE),
            aes(label = "*", y = Importance), 
            vjust = .3, size = 7, color = "red", fontface = "bold") +
  scale_fill_manual(values = c("Structure" = "#d90368",
                               "Topography" = "#820263","Landscape" = "#fb8b24"))+
  scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1.0), guide = "none") +
  facet_wrap(~Species_Wrapped, ncol = 5) +
  coord_cartesian(ylim = c(0, 75)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, size = 0.25)
  ) +
  labs(y = "Importance",
       x = "Variable")

p3

ggsave("/Users/abarenb1/Documents/PhD/Chapter3/PaperFigures/VariableImportanceSpecies.png", 
       width = 30, height = 30, units = "cm")

# Filter for only 4 specific species
selected_species <- c("Dwarf Honeyguide","Handsome Spurfowl","Oberlander's Ground-thrush", "Purple-breasted Sunbird", 
                      "Strange Weaver", "Neumann's Warbler")

combined_long_filtered <- combined_grouped %>%
  filter(Species %in% selected_species) %>%
  mutate(Species_Wrapped = str_replace_all(Species, " ", "\n"))


p4 <- ggplot(combined_long_filtered, aes(fill=Variable_Group,
                                         x = reorder(Variable_Label, as.numeric(factor(Variable_Group))), 
                                         y = Importance,
                                         group=factor(Variable_Group,ordered=TRUE))) +
  geom_bar(stat = "identity", 
           aes(alpha = Is_Top),
           color = "black", 
           size = 0.3) +
  geom_text(data = filter(combined_long_filtered, Is_Top == TRUE),
            aes(label = "*", y = Importance), 
            vjust = .3, size = 7, color = "red", fontface = "bold") +
  scale_fill_manual(values = c("Structure" = "#d90368",
                               "Topography" = "#820263","Landscape" = "#fb8b24"))+
  scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1.0), guide = "none") +
  facet_wrap(~Species_Wrapped, ncol = 3) +
  coord_cartesian(ylim = c(0, 40)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, size = 0.25)
  ) +
  labs(y = "Importance",
       x = "Variable")

p4

##Graph Variable Importance Compared by Grouping

group3_vars <- c("GEDI_biomass_filled", "Bwindi_Heights_2mMaxarResample", 
                 "RH98_classified_2023_Feb23_2026", 
                 "ETH_GlobalCanopyHeightSD_2020_10m_v1", "WSCI")

group2_vars <- c("Bwindi_2mFilledDSM_MaxarResample", 
                 "Bwindi_Slope_2mMaxarResample","Bwindi_tsri")

group1_vars <- c("DistanceToBorder","Precipitation_WorldClim","Dist_Wat_25m_Uganda_250mDist","LANDSAT_08_SR_07.01.2022.01.30.2023_NDVI")

group3_label <- "Structure"
group2_label <- "Topography"
group1_label <- "Landscape"

# Create grouped data
# Create grouped data
combined_grouped <- combined_long %>%
  mutate(Variable_Group = ifelse(Variable %in% group1_vars, 
                                 group1_label, 
                                 ifelse(Variable %in% group2_vars,
                                        group2_label,
                                        group3_label))) %>%
  group_by(Species, Species_Wrapped, Variable_Group) %>%
  summarise(Combined_Importance = sum(Importance), .groups = "drop")

# Create side-by-side bar chart
p_grouped <- ggplot(combined_grouped, 
                    aes(x = Species, 
                        y = Combined_Importance, 
                        fill = Variable_Group)) +
  geom_bar(stat = "identity", 
           position = "stack",
           color = "black", 
           size = 0.3,
           width = 0.7) +
  #geom_text(aes(label = round(Combined_Importance, 1)), 
  #          position = position_dodge(width = 0.8),
  #          vjust = -0.5, 
  #          size = 3.5,
  #          fontface = "bold") +
  scale_fill_manual(values = c("Structure" = "#d90368",
                               "Topography" = "#820263","Landscape" = "#fb8b24")) +
  scale_x_discrete(limits = c("Archer's Robin-chat","Black-faced Apalis","Dwarf Honeyguide","Grauer's Swamp-Warbler",
                              "Grauer's Warbler","Neumann's Warbler","Oberlander's Ground-thrush","Red-faced Woodland-Warbler",
                              "Red-throated Alethe","Rwenzori Apalis","Rwenzori Batis","Strange Weaver","Stripe-breasted Tit",
                              "Yellow-eyed Black-Flycatcher","Regal Sunbird","Blue-headed Sunbird",
                              "Purple-breasted Sunbird","Dusky Crimsonwing","Handsome Spurfowl","Grauer's Broadbill"))+
  theme_minimal() +
  geom_vline(xintercept = c(14.5,17.5,19.5), colour="blue", linetype = "longdash",lwd=1)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "italic"),
    axis.text.y = element_text(size = 12),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, size = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(y = "Importance",
       x = "Species") + theme(ylim=100)

p_grouped

##Compare AUC by species
setwd("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-03-26.SSDM/Species")
# Get all folders in the directory (each folder is a species)
all_folders <- list.dirs(path = getwd(), full.names = FALSE, recursive = FALSE)

# Remove any non-species folders if needed (e.g., hidden folders or system folders)
all_folders <- all_folders[!grepl("^\\.", all_folders)]  # Remove hidden folders starting with .

all_species_auc <- list()

# Read AUC data from all species folders
for (species in all_folders) {
  species_folder <- file.path(getwd(), species,"Tables")
  
  # Find the AlgoEval CSV file
  csv_files <- list.files(species_folder, 
                          pattern = "AlgoEval\\.csv", 
                          full.names = TRUE,
                          ignore.case = TRUE)
  
  if (length(csv_files) > 0) {
    csv_file <- csv_files[1]
    
    # Read the CSV file
    df <- read.csv(csv_file, row.names = 1, check.names = FALSE)
    
    # Extract AUC column
    if ("AUC" %in% colnames(df)) {
      auc_values <- as.numeric(df[, "AUC"])
      avg_auc <- mean(auc_values, na.rm = TRUE)
      all_species_auc[[species]] <- list(
        values = auc_values,
        mean = avg_auc,
        sd = sd(auc_values, na.rm = TRUE),
        min = min(auc_values, na.rm = TRUE),
        max = max(auc_values, na.rm = TRUE)
      )
      cat(paste("✓ Loaded AUC for:", species, "- Mean AUC:", round(avg_auc, 3), "\n"))
    } else {
      cat(paste("✗ No AUC column found for species:", species, "\n"))
      cat("Available columns:", paste(colnames(df), collapse = ", "), "\n")
    }
  } else {
    cat(paste("✗ No AlgoEval.csv file found for species:", species, "\n"))
  }
}

# Create summary dataframe
auc_summary <- data.frame(
  Species = names(all_species_auc),
  Mean_AUC = sapply(all_species_auc, function(x) x$mean)
) %>% arrange(desc(Mean_AUC))

# Create visualization - Bar plot
auc_summary <- auc_summary %>%
  mutate(Species_Wrapped = str_replace_all(Species, " ", "\n"))

### Repeat for LTSR
setwd("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_MAR/Landscape-Topography-Structure-Elevation.SSDM/Species")
# Get all folders in the directory (each folder is a species)
all_folders <- list.dirs(path = getwd(), full.names = FALSE, recursive = FALSE)

# Remove any non-species folders if needed (e.g., hidden folders or system folders)
all_folders <- all_folders[!grepl("^\\.", all_folders)]  # Remove hidden folders starting with .

all_species_auc <- list()

# Read AUC data from all species folders
for (species in all_folders) {
  species_folder <- file.path(getwd(), species,"Tables")
  
  # Find the AlgoEval CSV file
  csv_files <- list.files(species_folder, 
                          pattern = "AlgoEval\\.csv", 
                          full.names = TRUE,
                          ignore.case = TRUE)
  
  if (length(csv_files) > 0) {
    csv_file <- csv_files[1]
    
    # Read the CSV file
    df <- read.csv(csv_file, row.names = 1, check.names = FALSE)
    
    # Extract AUC column
    if ("AUC" %in% colnames(df)) {
      auc_values <- as.numeric(df[, "AUC"])
      avg_auc <- mean(auc_values, na.rm = TRUE)
      all_species_auc[[species]] <- list(
        values = auc_values,
        mean = avg_auc,
        sd = sd(auc_values, na.rm = TRUE),
        min = min(auc_values, na.rm = TRUE),
        max = max(auc_values, na.rm = TRUE)
      )
      cat(paste("✓ Loaded AUC for:", species, "- Mean AUC:", round(avg_auc, 3), "\n"))
    } else {
      cat(paste("✗ No AUC column found for species:", species, "\n"))
      cat("Available columns:", paste(colnames(df), collapse = ", "), "\n")
    }
  } else {
    cat(paste("✗ No AlgoEval.csv file found for species:", species, "\n"))
  }
}

# Create summary dataframe
auc_summaryLTSRE <- data.frame(
  Species = names(all_species_auc),
  Mean_AUC = sapply(all_species_auc, function(x) x$mean)
) %>% arrange(desc(Mean_AUC))

# Create visualization - Bar plot
auc_summaryLTSRE <- auc_summaryLTSRE %>%
  mutate(Species_Wrapped = str_replace_all(Species, " ", "\n"))

##Get difference in Average AUC

auc_summaryBOTH<-within(merge(auc_summaryLTSRE,auc_summary,by="Species_Wrapped"), {
  Mean_AUC <- Mean_AUC.x - Mean_AUC.y
  })[,c("Mean_AUC","Species_Wrapped")]


ggplot(auc_summaryBOTH, aes(x = reorder(Species_Wrapped, Mean_AUC), y = Mean_AUC,, fill = Mean_AUC > 0)) +
  geom_bar(stat = "identity", ) +
  #geom_text(aes(label = round(Mean_AUC, 3)), 
  #          vjust = -0.5, hjust = -0.2, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "#555555"), guide = "none") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  labs(title = "LAND vs LTSR-E",
       x = "Species",
       y = "Difference Mean AUC") +
  ylim(min(auc_summaryBOTH$Mean_AUC) * 1.1, max(auc_summaryBOTH$Mean_AUC) * 1.1)

#Mapping Uncertainty
##################################################################
rastlist <- list.files(path = "/Users/abarenb1/Documents/PhD/Chapter3/SSDM/GEDISSDM_5StackModel.SSDM/Species", pattern='Uncertainty', 
                       all.files=TRUE, full.names=TRUE)
rastlistshort <- list.files(path = "/Users/abarenb1/Documents/PhD/Chapter3/SSDM/GEDISSDM_5StackModel.SSDM/Species", pattern='Uncertainty', 
                            all.files=TRUE, full.names=FALSE)

rastlistNON <- list.files(path = "/Users/abarenb1/Documents/PhD/Chapter3/SSDM/NonGEDISSDM_5StackModel.SSDM/Species", pattern='Uncertainty', 
                       all.files=TRUE, full.names=TRUE)
rastlistshortNON <- list.files(path = "/Users/abarenb1/Documents/PhD/Chapter3/SSDM/NonGEDISSDM_5StackModel.SSDM/Species", pattern='Uncertainty', 
                            all.files=TRUE, full.names=FALSE)
#
##import all raster files in folder using lapply
allrasters <- lapply(rastlist, rast)
allrastersNON <- lapply(rastlistNON, rast)

names(allrasters) <- tools::file_path_sans_ext(rastlistshort)
names(allrastersNON) <- tools::file_path_sans_ext(rastlistshortNON)

names(allrasters) <- gsub("UncertaintyGEDIModel", "", names(allrasters))
names(allrastersNON) <- gsub("UncertaintyGEDIModel", "", names(allrastersNON))


plot_rasters_facet <- function(raster_list) {
  # Combine all rasters into one data frame
  df_combined <- bind_rows(lapply(seq_along(raster_list), function(i) {
    df <- as.data.frame(raster_list[[i]], xy = TRUE, na.rm = FALSE)
    
    # Rename the last column to "value"
    colnames(df)[ncol(df)] <- "value"
    
    df$name <- if(!is.null(names(raster_list)[i]) && names(raster_list)[i] != "") {
      # Replace spaces with newlines for wrapping at every word
      gsub(" ", "\n", names(raster_list)[i])
    } else {
      paste0("Raster_", i)
    }
    
    return(df)
  }))
  
  # Remove NAs
  df_combined <- df_combined %>% 
    filter(!is.na(value)) %>%
    drop_na()
  
  ggplot(df_combined, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "cividis",name="Uncertainty") +
    facet_wrap(~name, ncol = 5) +
    coord_equal() +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    )
}


plot <- plot_rasters_facet(allrasters)
plot

ggsave(plot = plot, width = 5, height = 5, dpi = 300, filename = "/Users/abarenb1/Documents/PhD/Chapter3/GEDIUncertaintySpec.png")


##Uncertainty Correlation
#################################################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)

matching_tifs <- c("GEDI_biomass_filled","Dist_Wat_25m_Uganda_250mDist","Precipitation_WorldClim",
                   "WSCI","Bwindi_Heights_2mMaxarResample","Bwindi_2mFilledDSM_MaxarResample",
                   "Bwindi_Slope_2mMaxarResample","ETH_GlobalCanopyHeightSD_2020_10m_v1","Bwindi_tsri",
                   "RH98_classified_2023_Mar12_2026","LANDSAT_08_SR_07-01-2022-01-30-2023_NDVI")

# Human-readable labels aligned to matching_tifs order
raster_labels <- c("Biomass", "DistanceWater", "Precipitation", "WSCI", "CanopyHeight",
                   "DSM", "Slope", "SDHeight", "TSRI", "RH98", "NDVI")

# Load all rasters, subsetting WSCI to "mean" layer only
rasterList <- lapply(seq_along(matching_tifs), function(i) {
  r <- rast(paste0("/Users/abarenb1/Documents/PhD/Chapter3/VariablesResample/", matching_tifs[i], ".tif"))
  
  # Pull only the "mean" layer from WSCI
  if (matching_tifs[i] == "WSCI") {
    if ("Mean WSCI estimate" %in% names(r)) {
      r <- r[["Mean WSCI estimate"]]
    } else {
      stop("No layer named 'mean' found in WSCI raster. Available layers: ", paste(names(r), collapse = ", "))
    }
  }
  
  return(r)
})

# Use the Maxar to mask
mask_raster <- rasterList[[5]]

# Mask all rasters by Maxar extent
rasterList <- lapply(seq_along(rasterList), function(i) {
  mask(rasterList[[i]], mask_raster)
})

# Extract each raster as a named column, joined on x and y
df_list <- lapply(seq_along(rasterList), function(i) {
  df <- as.data.frame(rasterList[[i]], xy = TRUE, na.rm = FALSE)
  colnames(df) <- c("x", "y", raster_labels[i])  # name the value column directly
  return(df)
})

# Join all by x and y coordinates
df_combined <- Reduce(function(a, b) full_join(a, b, by = c("x", "y")), df_list)

# Remove rows where ALL raster columns are NA
df_combined <- df_combined %>%
  filter(if_any(all_of(raster_labels), ~ !is.na(.)))

uncertainGEDI <-rast("/Users/abarenb1/Documents/PhD/Chapter3/SSDM/GEDISSDM_5StackModel.SSDM/Species/UncertaintyGEDIModelDwarf Honeyguide.tif")
uncertainNON<-rast("/Users/abarenb1/Documents/PhD/Chapter3/SSDM/NonGEDISSDM_5StackModel.SSDM/Species/UncertaintyGEDIModelDwarf Honeyguide.tif")

#uncert_values <- values(uncertainGEDI)
#uncert_valuesNON <- values(uncertainNON)
#
#df <- data.frame(
#  elevation = elev_values,
#  uncertaintyNON = uncert_valuesNON,
#  uncertaintyGEDI = uncert_values
#)
#df <- na.omit(df)
#
#colnames(df)<-c('elevation','NON','GEDI')
#head(df)

#Plot.
ggplot(df, aes(x = elevation, y = GEDI)) +
  geom_point(alpha = 0.15, size = 0.3, color = "steelblue") +
  geom_smooth(method = "gam", color = "red", fill = "red", alpha = 0.2) +
  labs(title = "GEDI",
       x = "Elevation",
       y = "Uncertainty") +
  theme_minimal()+ylim(0,0.25)

ggplot(df, aes(x = elevation, y = NON)) +
  geom_point(alpha = 0.15, size = 0.3, color = "steelblue") +
  geom_smooth(method = "gam", color = "red", fill = "red", alpha = 0.2) +
  labs(title = "Non GEDI",
       x = "Elevation",
       y = "Uncertainty") +
  theme_minimal()+ylim(0,0.25)


presenceGEDI <-rast("/Users/abarenb1/Documents/PhD/Chapter3/SSDM_FEB/Landscape-Topography-Structure-02-27-26.SSDM/Stack/Rasters/Diversity.tif")
presenceNON<-rast("/Users/abarenb1/Documents/PhD/Chapter3/SSDM/NonGEDISSDM_5StackModel.SSDM/Species/NonGEDIModelDwarf Honeyguide.tif")

pres_values <- values(presenceGEDI)
pres_valuesNON <- values(presenceNON)

dfpres <- data.frame(
  elevation = elev_values,
  fhd= fhd_values,
  wsci = wsci_values,
  presGEDI = pres_values
)
dfpres <- na.omit(dfpres)

colnames(dfpres)<-c('elevation','fhd','wsci','GEDI')
head(dfpres)

a<-ggplot(dfpres, aes(x = elevation, y = GEDI)) +
  geom_point(alpha = 0.15, size = 0.3, color = "orange") +
  geom_smooth(method = "gam", color = "red", fill = "red", alpha = 0.2) +
  labs(x = NULL, y = "Species Observed") +
  theme_minimal() +
  theme(text = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

b<-ggplot(dfpres, aes(x = fhd, y = GEDI)) +
  geom_point(alpha = 0.15, size = 0.3, color = "purple") +
  geom_smooth(method = "gam", color = "red", fill = "red", alpha = 0.2) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

c<-ggplot(dfpres, aes(x = wsci, y = GEDI)) +
  geom_point(alpha = 0.15, size = 0.3, color = "#1EB3A5") +
  geom_smooth(method = "gam", color = "red", fill = "red", alpha = 0.2) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

###DensityPlots
#############################################################################################

d<-ggplot(df_combined, aes(x = DSM))+geom_density(alpha=0.3,color="orange", fill='orange')+
  labs(x = "DSM", y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = 'none')  

d

e<-ggplot(df_combined, aes(x = RH98))+geom_density(alpha=0.3,color="purple", fill='purple')+
  labs(x = "RH98(m)", y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
e  

f<- ggplot(df_combined, aes(x = WSCI))+geom_density(alpha=0.3,color="#1EB3A5", fill='#1EB3A5')+
  labs(x = "WSCI", y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.ticks.y = element_blank())

f

g<- ggplot(df_combined, aes(x = Biomass))+geom_density(alpha=0.3,color="#F527DA", fill='#F527DA')+
  labs(x = "Biomass(Mg/ha)", y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.ticks.y = element_blank())

g

h<- ggplot(df_combined, aes(x = DistanceWater))+geom_density(alpha=0.3,color="#259BC2", fill='#259BC2')+
  labs(x = "Distance Water(m)", y = "Relative Frequency") +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.ticks.y = element_blank())

h

i<- ggplot(df_combined, aes(x = Precipitation))+geom_density(alpha=0.3,color="#C2BA25", fill='#C2BA25')+
  labs(x = "Precipitation(mm)", y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.ticks.y = element_blank())

i

j<- ggplot(df_combined, aes(x = CanopyHeight))+geom_density(alpha=0.3,color="#C25425", fill='#C25425')+
  labs(x = "Canopy Height(m)", y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.ticks.y = element_blank())

j

k<- ggplot(df_combined, aes(x = Slope))+geom_density(alpha=0.3,color="#C22554", fill='#C22554')+
  labs(x = "Slope", y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.ticks.y = element_blank())

k

l<- ggplot(df_combined, aes(x = SDHeight))+geom_density(alpha=0.3,color="#4F2569", fill='#4F2569')+
  labs(x = "SD Height", y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.ticks.y = element_blank())

l

m<- ggplot(df_combined, aes(x = TSRI))+geom_density(alpha=0.3,color="#516925", fill='#516925')+
  labs(x = "TSRI", y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.ticks.y = element_blank())

m

n<- ggplot(df_combined, aes(x = NDVI))+geom_density(alpha=0.3,color="#46D454", fill='#46D454')+
  labs(x = "NDVI", y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 15),
        legend.position = 'none',
        axis.ticks.y = element_blank())

n




#(d | k | m | n)/(h | i | j | e)/(l | f | g |n)
(d | k | m | n | i)/(h | j | e| l |f)/( g |n|n|n|n)

ggsave("/Users/abarenb1/Documents/PhD/Chapter3/PaperFigures/RelativeFrequenciesAll.png", 
       width = 20, height = 20, units = "cm")

ggplot(dfpres, aes(x = elevation, y = NON)) +
  geom_point(alpha = 0.15, size = 0.3, color = "orange") +
  geom_smooth(method = "gam", color = "red", fill = "red", alpha = 0.2) +
  labs(title = "Non GEDI",
       x = "Elevation",
       y = "Presence") +
  theme_minimal()+ylim(0,0.7)

a<-ggplot(df, aes(x = Bwindi_2mFilledDSM_MaxarResample, y = UncertaintyGEDIModelStrange.Weaver)) +
  geom_point(alpha = 0.3) +
  #geom_smooth(method = "lm", color = "red") +
  labs(title = "GEDI",
       x = "Elevation",
       y = "Uncertainty") +
  theme_minimal()
a





b<-ggplot(dfNON, aes(x = Bwindi_2mFilledDSM_MaxarResample, y = UncertaintyGEDIModelStrange.Weaver)) +
  geom_point(alpha = 0.3) +
  #geom_smooth(method = "lm", color = "red") +
  labs(title = "Non GEDI",
       x = "Elevation",
       y = "Uncertainty") +
  theme_minimal()
b

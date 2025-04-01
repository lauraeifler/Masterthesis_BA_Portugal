# remove all variables
rm(list=ls())

# load packages
library(raster)
library(caret)
library(sf)
#library(mapview)
library(dplyr)
library(CAST)
library(randomForest)
library(lme4)
library(performance)
library(ggspatial) 
library(ggplot2)
library(quantreg)
library(mgcv)
library(cluster)
library(rnaturalearth)
library(rnaturalearthdata)
library(permimp)
library(vegan)
library(ineq)

setwd("~/Masterthesis/Data")
path_fig <- "~/Masterthesis/Figures/"
################################################################################
data <- read.csv('Data_BA_forest.csv')

table(data['lccs'])
table(data['PA_status'])

# Create the interaction table
interaction_table <- table(data$lccs, data$PA_status)

# Print the table
print(interaction_table)

##### BALANCE #####
# Calculate the minimum number of samples
min_samples <- min(interaction_table)
print(min_samples)


# Stratified sampling
balanced_data <- data %>%
  group_by(lccs, PA_status) %>%
  slice_sample(n = min_samples) %>%
  ungroup()

# Check the new balanced data
table(balanced_data$lccs, balanced_data$PA_status)

write.csv(balanced_data, 'Data_BA_forest_balanced.csv')

################################################################################
data <- read.csv('Data_BA_forest_balanced.csv')


## clean, because first three columns some indexes
data <- data %>%
  select(-(1:3)) %>%
  mutate(PA_status = as.factor(PA_status),
         lccs = as.factor(lccs))

## normalization of Lat and Lon
data$Latitude_norm <- (data$Latitude - mean(data$Latitude, na.rm = TRUE)) / 
  sd(data$Latitude, na.rm = TRUE)

data$Longitude_norm <- (data$Longitude - mean(data$Longitude, na.rm = TRUE)) / 
  sd(data$Longitude, na.rm = TRUE)



## PA shape 
portugal_sf <- ne_countries(scale = "medium", country = "Portugal", returnclass = "sf")
# Define a bounding box for mainland Portugal (approximate)
mainland_bbox <- st_bbox(c(xmin = -10, xmax = -6, ymin = 36.8, ymax = 42), crs = st_crs(portugal_sf))

# Crop Portugal to the bounding box (removes islands)
pt_sf <- st_crop(portugal_sf, mainland_bbox)



## CLUSTER
## create sf dataset for later accounting for spatial autocorrelation
set.seed(84)  # For reproducibility
clusters <- kmeans(data[,  c("Longitude", "Latitude")], centers =10)$cluster  # Adjust centers as needed

# Add clusters to the dataset
data$cluster <- as.factor(clusters)


##### CV #####
##### Random folds, balanced
set.seed(42)
# data$fold <- sample(1:10, size = nrow(data), replace = TRUE)
randomfold <- createFolds(data$PA_status, k = 10, list = TRUE, returnTrain = FALSE)
data$randomfold <- NA

# Assign fold numbers to each row
for (i in seq_along(randomfold)) {
  data$randomfold[randomfold[[i]]] <- i  # Assign fold number i
}


## balanced folds based on PA status AND LC
set.seed(42)
# Create a new factor combining PA_status and lccs
data$PA_lccs <- interaction(data$PA_status, data$lccs, drop = TRUE)

randomfold_b <- createFolds(data$PA_lccs, k = 10, list = TRUE, returnTrain = FALSE)
data$randomfold_b <- NA

# Assign fold numbers to each row
for (i in seq_along(randomfold_b)) {
  data$randomfold_b[randomfold_b[[i]]] <- i  # Assign fold number i
}


table(data$randomfold, data$PA_status)
table(data$randomfold_b, data$PA_status)


#### CLUSTERFOLDS 
## balances clusterfold
data$clusterfold <- NA

for (fold in unique(data$cluster)) {
  fold_indices <- which(data$cluster == fold)
  fold_data <- data[fold_indices, ]  # Subset fold data
  
  # Find the minimum count between PA_status classes
  min_count <- min(table(fold_data$PA_status))
  
  # Downsample both classes to match the smaller count
  # if balanced PA and LC change teh data$PA_lccs
  balanced_indices <- unlist(lapply(split(fold_indices, data$PA_status[fold_indices]), function(idx) {
    sample(idx, min_count)
  }))
  
  # Assign fold values only for the balanced subset
  data$clusterfold[balanced_indices] <- fold
}


###### balanced for PA AND LC
data$clusterfold_b <- NA

for (fold in unique(data$cluster)) {
  fold_indices <- which(data$cluster == fold)
  fold_data <- data[fold_indices, ]  # Subset fold data
  
  ## balanced PA and LC
  # Create combined factor for stratification (PA_status + lccs)
  # fold_data$PA_lccs <- interaction(fold_data$PA_status, fold_data$lccs, drop = TRUE)
  # Find the minimum count of any PA_status + lccs combination
  min_count <- min(table(fold_data$PA_lccs))
  
  # Downsample both classes to match the smaller count
  # if balanced PA and LC change teh data$PA_lccs
  balanced_indices <- unlist(lapply(split(fold_indices, data$PA_lccs[fold_indices]), function(idx) {
    sample(idx, min_count)
  }))
  
  # Assign fold values only for the balanced subset
  data$clusterfold_b[balanced_indices] <- fold
}


# Convert to factor and remove NA (only keep balanced rows)
# dataTrain <- dataTrain[!is.na(dataTrain$clusterfold_balanced), ]
data$clusterfold_b <- as.factor(data$clusterfold_b)

###### tables 
# Check balance per fold
table(data$clusterfold_b, data$PA_status)

# Check the distribution of folds
table(data$clusterfold) #, data$PA_status)
table(data$clusterfold_b) # , data$lccs

table(data$fold, data$lccs)
table(data$fold_balance, data$lccs)


#### test homogeneity per folds

class_counts <- table(data$fold, data$PA_status)

# Calculate homogeneity per fold
homogeneity_metrics <- apply(class_counts, 1, function(x) {
  mean_x <- mean(x)       # Mean of the counts per fold
  sd_x <- sd(x)           # Standard deviation
  min_max_ratio <- min(x) / max(x)  # Min/Max ratio
  
  homogeneity_score <- 1 - (sd_x / mean_x)  # Higher means more homogeneous
  
  return(c(Homogeneity = homogeneity_score, MinMaxRatio = min_max_ratio))
})

# Convert to data frame for readability
homogeneity_df <- as.data.frame(t(homogeneity_metrics))
homogeneity_df$Fold <- rownames(homogeneity_df)
homogeneity_df

mean(homogeneity_df$Homogeneity)


#############
## separate to training and test data
set.seed(84)
trainingids <- createDataPartition(data$PA_lccs, list=FALSE, p=0.8) # (80% training data can also be incraesed)
dataTrain <- data[trainingids,]
dataTest <- data[-trainingids,]

table(dataTrain$PA_status)

# Check the distribution of folds
table(dataTest$clusterfold) #, dataTrain$lccs)
table(dataTrain$clusterfold_b) #, dataTrain$lccs)

table(dataTrain$randomfold, dataTrain$lccs)
table(dataTrain$randomfold_b, dataTrain$lccs)

table(data$PA_lccs)

###################
## data as spatial feature
data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs=4326)
# data_sf
data_sf$clusterfold <- factor(data_sf$clusterfold, levels = 1:10)
data_sf$clusterfold_b <- factor(data_sf$clusterfold_b, levels = 1:10)

#######################################
library(spdep)
### testing for spatial autocorrelation in the data with Moran's I test

# Assume data_sf is your spatial points data as an sf object
# Create a spatial weights matrix (neighbors)
coords <- st_coordinates(data_sf)
nb <- knn2nb(knearneigh(coords, k = 10))  # Adjust k for neighbors, e.g., 4 nearest neighbors

# Compute spatial weights matrix
lw <- nb2listw(nb, style = "W")

# Calculate Moran's I
moran_test <- moran.test(data_sf$Burned_Area, lw)
print(moran_test)

########################################
## Visualize cluster
clusters_sf <- data_sf %>% group_by(cluster)
# For scale bar and north arrow

# Create convex hulls for each cluster
cluster_polygons <- clusters_sf %>%
  summarise(geometry = st_convex_hull(st_union(geometry)))

ggplot() +
  geom_sf(data = pt_sf, fill = NA, color = "black", size = 0.8) +
  # Plot cluster polygons
  geom_sf(data = cluster_polygons, aes(fill = cluster), alpha = 0.4, color = "black") +
  
  # Plot individual points
  geom_sf(data = data_sf, aes(color = cluster), size = 1) +
  
  # Customize appearance
  scale_fill_viridis_d(name = "Clusters") +
  scale_color_viridis_d(name = "Clusters") +
  
  # Add map elements
  #annotation_scale(location = "br") +
  annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering) +
  
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))+
  # Explicitly set longitude breaks for 9°W, 8°W, and 7°W
  scale_x_continuous(breaks = c(-9, -8, -7), labels = c("9°W", "8°W", "7°W"))+
  
  coord_sf(expand = FALSE)

ggsave(filename = paste0(path_fig, "Map_PT_cluster_dist.png"), plot = last_plot(),
       width = 5, height = 7, dpi = 300, limitsize = FALSE)


## visualize random
#### coloring folds random sampling, no account for spatial autocorrelation
fold_colors <- c("navy", "maroon", "yellowgreen", "purple", "orange", "yellow", "darkcyan", 
                 "deeppink", "grey", "darkgreen")
ggplot() +
  geom_sf(data = pt_sf, fill = NA, color = "black", size = 0.8) +
  # Plot individual points
  geom_sf(data = data_sf, aes(color = as.factor(randomfold)), size = 1) +
  
  # Add map elements
  #annotation_scale(location = "br") +
  annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering) +
  
  scale_color_manual(values = fold_colors, name = "Group") +
  
  # Explicitly set longitude breaks for 9°W, 8°W, and 7°W
  scale_x_continuous(breaks = c(-9, -8, -7), labels = c("9°W", "8°W", "7°W")) +
  
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))

ggsave(filename = paste0(path_fig, "Map_PT_random_sampling_dist.png"), plot = last_plot(),
  width = 5, height = 8, dpi = 300)
  

## cluster folds
ggplot() +
  geom_sf(data = pt_sf, fill = NA, color = "black", size = 0.8) +
  
  # Plot individual points
  geom_sf(data = data_sf %>% filter(!is.na(clusterfold)), aes(color = as.factor(clusterfold)), size = 1) +
  
  # Add map elements
  #annotation_scale(location = "br") +
  annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering) +
  
  scale_color_manual(values = fold_colors, name = "Cluster") +
  
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude") +
  
  # Explicitly set longitude breaks for 9°W, 8°W, and 7°W
  scale_x_continuous(breaks = c(-9, -8, -7), labels = c("9°W", "8°W", "7°W")) +
  
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))

ggsave(filename = paste0(path_fig, "Map_PT_spatial_sampling_dist.png"), plot = last_plot(),
       width = 5, height = 8, dpi = 300)


#### balanced clusterfolds ####
ggplot() +
  geom_sf(data = pt_sf, fill = NA, color = "black", size = 0.8) +
  
  # Plot individual points
  geom_sf(data = na.omit(data_sf), aes(color = as.factor(clusterfold_b)), size = 1) +
  
  # Add map elements
  #annotation_scale(location = "br") +
  annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering) +
  
  scale_color_manual(values = fold_colors, name = "Cluster") +
  
  # Explicitly set longitude breaks for 9°W, 8°W, and 7°W
  scale_x_continuous(breaks = c(-9, -8, -7), labels = c("9°W", "8°W", "7°W")) +
  
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))

ggsave(filename = paste0(path_fig, "Map_PT_spatial_balanced_sampling_dist.png"), plot = last_plot(),
       width = 5, height = 8, dpi = 300)



### all data points
forest_colors <- c("yellowgreen", "darkolivegreen4", "darkgreen")
ggplot() +
  geom_sf(data = pt_sf, fill = NA, color = "black", size = 0.8) +
  # Plot individual points
  geom_sf(data = data_sf, aes(color = as.factor(lccs)), size = 1) +
  
  # Add map elements
  #annotation_scale(location = "br") +
  annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering) +
  
  scale_color_manual(values = forest_colors, name = "LC class") +
  
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")

table(dataTrain$lccs, dataTrain$PA_status)



#### BA and PA histogram and Boxplot ####
# Plot the stacked histogram
ggplot(dataTrain, aes(x = Burned_Area, fill = as.factor(PA_status))) +
  geom_histogram(binwidth = 5, position = "stack", color = "black") +
  scale_fill_manual(values = c("0" = "lightblue", "1" = "darkcyan"),
                    name = "Protection Status",
                    labels = c("Non-protected", "Protected")) +
  labs(x = "Burned Area [%]", y = "Count") +
  theme_minimal() +
  theme(legend.position = c(0, 1),        # Position the legend in the top left corner
        legend.justification = c(0, 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))
ggsave(filename = paste0(path_fig, "Histogram_BA_PA_traindata.png"), plot = last_plot(),
       width = 8, height = 4, dpi = 300)



ggplot(dataTrain, aes(x = as.factor(PA_status), y = Burned_Area, fill = as.factor(PA_status))) +
  geom_boxplot(width = 0.3) +  # Narrower boxes
  scale_fill_manual(
    values = c("0" = "lightblue", "1" = "darkcyan")) +
  labs(x = "Protection Status", y = "Burned Area [%]") +
  scale_x_discrete(labels = c("0" = "Non-Protected", "1" = "Protected")) +  # Change x-axis labels
  theme_minimal() +
  # Add median line in orange and text next to it
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 1), group = PA_status),
               color = "purple4", size = 6, fontface = "bold", vjust = -0.5, hjust = 0.5) +  
  theme(
    axis.text.x = element_text(size = 14, vjust = 0.5),  # X-axis labels larger
    axis.text.y = element_text(size = 14),  # Y-axis tick labels larger
    axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
    axis.title.y = element_text(size = 16, face = "bold"),  # Y-axis title larger
    axis.ticks.length = unit(0.3, "cm"),  # Extend tick marks
    legend.position = "none"  # Remove legend
  )

ggsave(filename = paste0(path_fig, "Boxplot_BA_PA_traindata.png"), plot = last_plot(),
       width = 8, height = 4, dpi = 300)


###############################################

################################################################################
## Quantile regression


# Define the models
models <- list()
# models[["M1"]] <- "Burned_Area ~ fwi_norm"
# models[["M2"]] <- "Burned_Area ~ fwi_norm * PA_status"
# models[["M3"]] <- "Burned_Area ~ fwi_norm * Latitude * Longitude"
# models[["M4"]] <- "Burned_Area ~ fwi_norm + Latitude * PA_status"
models[["M5"]] <- "Burned_Area ~ fwi_norm + PA_status * lccs"
# models[["M6"]] <- "Burned_Area ~ PA_status"
# models[["M7"]] <- "Burned_Area ~ PA_status * lccs"
# models[["M8"]] <- "Burned_Area ~ PA_status + Latitude_norm + Longitude_norm"


print(cor(dataTrain[, c("PA_status", "lccs")], use = "complete.obs"))

#### FINAL approach
set.seed(84)
results <- list()

for (mod_name in names(models)) {
  for (quant in c(0.25, 0.5, 0.75)) {
    
    coef_list <- list()  # Stores coefficients for each fold
    sd_list <- list()     # Stores coefficient SDs per fold
    r2_train_list <- c()  # Stores training R² per fold
    r2_test_list <- c()   # Stores test R² per fold
    rmse_list <- c()      # Stores RMSE per fold
    
    for (fold in 1:10) {
      # Split into train and test
      train_data <- dataTrain %>% filter(!is.na(clusterfold) & clusterfold != fold)
      test_data <- dataTrain %>% filter(!is.na(clusterfold) & clusterfold == fold)
      
      # Skip fold if there are not enough data points
      # if (nrow(train_data) < 2 | nrow(test_data) < 2) {
        # next  
      # }
      
      # Fit quantile regression model
      mod_formula <- as.formula(models[[mod_name]])
      model <- rq(mod_formula, data = train_data, tau = quant)
      
      # Store coefficients & SD per fold
      coef_val <- coef(model)
      coef_list[[fold]] <- coef_val
      sd_list[[fold]] <- summary(model)$coefficients[, 2]  # Extract SDs
      
      # Predict on training data (for R² train)
      train_predictions <- predict(model, newdata = train_data)
      valid_train_idx <- !is.na(train_predictions) & !is.na(train_data$Burned_Area)
      
      if (sum(valid_train_idx) > 1) {
        ss_total_train <- sum((train_data$Burned_Area[valid_train_idx] - mean(train_data$Burned_Area[valid_train_idx]))^2)
        ss_residual_train <- sum((train_data$Burned_Area[valid_train_idx] - train_predictions[valid_train_idx])^2)
        r2_train <- ifelse(ss_total_train > 0, 1 - (ss_residual_train / ss_total_train), NA)
      } else {
        r2_train <- NA
      }
      r2_train_list <- c(r2_train_list, r2_train)
      
      # Predict on test data (for R² test & RMSE)
      predictions <- predict(model, newdata = test_data)
      valid_test_idx <- !is.na(predictions) & !is.na(test_data$Burned_Area)
      
      if (sum(valid_test_idx) > 1) {
        ss_total_test <- sum((test_data$Burned_Area[valid_test_idx] - mean(test_data$Burned_Area[valid_test_idx]))^2)
        ss_residual_test <- sum((test_data$Burned_Area[valid_test_idx] - predictions[valid_test_idx])^2)
        r2_test <- ifelse(ss_total_test > 0, 1 - (ss_residual_test / ss_total_test), NA)
        rmse <- sqrt(mean((test_data$Burned_Area[valid_test_idx] - predictions[valid_test_idx])^2, na.rm = TRUE))
      } else {
        r2_test <- NA
        rmse <- NA
      }
      
      r2_test_list <- c(r2_test_list, r2_test)
      rmse_list <- c(rmse_list, rmse)
    }
    
    # Compute summary statistics
    coef_df <- do.call(rbind, coef_list)  # Combine coefficients from each fold
    sd_df <- do.call(rbind, sd_list)      # Combine SDs from each fold
    
    coef_summary <- data.frame(
      Variable = names(coef_val),
      Mean_Coef = apply(coef_df, 2, mean, na.rm = TRUE),
      Mean_SD = apply(sd_df, 2, mean, na.rm = TRUE)
    )
    
    # Store all results
    results[[as.character(quant)]] <- list(
      coef_summary = coef_summary,
      Mean_R2_Train = mean(r2_train_list, na.rm = TRUE),
      Mean_R2_Test = mean(r2_test_list, na.rm = TRUE),
      Mean_RMSE = mean(rmse_list, na.rm = TRUE)
    )
  }
}

for (quant in c(0.25, 0.5, 0.75)) {
  cat(sprintf("\n--- Results for Quantile: Q%2.0f ---\n", 100 * quant))
  
  # Extract the coefficient summary and mean R² values for this quantile
  coef_summary <- results[[as.character(quant)]]$coef_summary
  mean_r2_train <- results[[as.character(quant)]]$Mean_R2_Train
  mean_r2_test <- results[[as.character(quant)]]$Mean_R2_Test
  mean_rmse <- results[[as.character(quant)]]$Mean_RMSE
  
  # Print coefficient summary
  cat("\nCoefficient Summary:\n")
  print(coef_summary)
  
  # Print performance metrics
  cat(sprintf("\nMean R² (Training Data): %.3f\n", mean_r2_train))
  cat(sprintf("Mean R² (Test Data): %.3f\n", mean_r2_test))
  cat(sprintf("Mean RMSE (Test Data): %.3f\n", mean_rmse))
}






################################################################################
## GAM

### FINAL APPROACH GAM
# Initialize an empty list to store model results
gam_results <- list()
valid_folds <- unique(na.omit(dataTrain$clusterfold)) 

# Lists to store coefficients and R-squared values
coefficients_list <- list()
adj_r2_values <- c()
rmse_list <- c()  
r2_list <- c()

set.seed(84)
# Loop through the 5 folds
for (fold in valid_folds) {
  
  cat(sprintf("\n--- Fold %s ---\n", fold))
  
  # Split the data: Use all data except the current fold for training
  train_data <- dataTrain %>% filter(!is.na(clusterfold) & clusterfold != fold)
  test_data <- dataTrain %>% filter(!is.na(clusterfold) & clusterfold == fold)
  
  # Fit the GAM model
  mod_gam <- gam(Burned_Area ~ s(Latitude_norm, bs = "cr") + s(Longitude_norm, bs = "cr") + 
                   PA_status, data = train_data, method = "REML") #*lccs + fwi_norm
  # mod_gam <- gam(Burned_Area ~ PA_status, data = train_data)
  
  # Store the model results
  gam_results[[paste0("Fold_", fold)]] <- mod_gam
  
  # Extract coefficients and adjusted R²
  coefficients_list[[paste0("Fold_", fold)]] <- coef(mod_gam)  # Store coefficients
  # adj_r2_values <- c(adj_r2_values, summary(mod_gam)$r.sq)  # Store adjusted R²
  if (!is.na(summary(mod_gam)$r.sq)) {
    adj_r2_values <- c(adj_r2_values, summary(mod_gam)$r.sq)
  }
  
  # Print summary for each fold
  print(summary(mod_gam))
  
  
  # Make predictions on the test set
  predictions <- predict(mod_gam, newdata = test_data)
  
  # Compute RMSE
  rmse_val <- sqrt(mean((test_data$Burned_Area - predictions)^2, na.rm = TRUE))
  
  # Compute R²
  ss_total <- sum((test_data$Burned_Area - mean(test_data$Burned_Area, na.rm = TRUE))^2, na.rm = TRUE)
  ss_residual <- sum((test_data$Burned_Area - predictions)^2, na.rm = TRUE)
  r2_val <- 1 - (ss_residual / ss_total)
  
  # Store RMSE and R² values
  rmse_list <- c(rmse_list, rmse_val)
  r2_list <- c(r2_list, r2_val)
  
  # Print summary for each fold
  print(summary(mod_gam))
  
  # Print RMSE and R² for the fold
  cat(sprintf("Fold %s: RMSE = %.3f, R² = %.3f\n", fold, rmse_val, r2_val))
}


#### Coefficients of model
# Convert coefficients_list into a data frame for easier computation
coefficients_df <- do.call(rbind, coefficients_list)  # Convert list to data frame

# Compute mean and standard deviation for each coefficient
mean_coefficients <- colMeans(coefficients_df, na.rm = TRUE)
sd_coefficients <- apply(coefficients_df, 2, sd, na.rm = TRUE)

# Compute mean adjusted R²
mean_adj_r2 <- mean(adj_r2_values, na.rm = TRUE)

# Print the results
cat("\n*** Summary of GAM Model Over All Folds ***\n")
cat("\nMean Coefficients:\n")
print(mean_coefficients)

cat("\nStandard Deviation of Coefficients:\n")
print(sd_coefficients)

cat(sprintf("\nMean Adjusted R²: %.4f\n", mean_adj_r2))

####### Prediction
# Compute mean RMSE and mean R² across folds
mean_rmse <- mean(rmse_list, na.rm = TRUE)
mean_r2 <- mean(r2_list, na.rm = TRUE)

# Print final results
cat("\n*** GAM Model Performance Across Folds ***\n")
cat(sprintf("Mean RMSE: %.3f\n", mean_rmse))
cat(sprintf("Mean R²: %.3f\n", mean_r2))


# Initialize an empty list to store coefficients
coef_list <- list()

# Loop through folds to extract coefficients
for (fold in names(gam_results)) {
  mod_gam <- gam_results[[fold]]  # Get the model
  
  # Extract parametric coefficients (not smooth terms)
  coef_vals <- coef(mod_gam)
  
  # Convert to a data frame
  df <- data.frame(
    Variable = names(coef_vals),
    Coefficient = coef_vals,
    Fold = fold
  )
  
  # Store in list
  coef_list[[fold]] <- df
}

# Combine all folds into one data frame
coef_df <- bind_rows(coef_list)

# Filter out intercept and smooth terms if needed
coef_df <- coef_df %>% filter(!grepl("s\\(", Variable))  # Remove smooth terms

# Extract numeric part from Fold column (e.g., "Fold_3" -> 3)
coef_df$Fold <- gsub("Fold_", "", coef_df$Fold)  # Remove "Fold_" part
coef_df$Fold <- as.numeric(coef_df$Fold)  # Convert to numeric

# Convert Fold to a factor with correct order (1 to 10)
coef_df$Fold <- factor(coef_df$Fold, levels = 1:10)

library(viridis)
library(RColorBrewer)

# Plot stacked histogram of coefficients
ggplot(coef_df, aes(x = Variable, y = Coefficient, fill = Fold)) +
  geom_bar(stat = "identity", position = "dodge") +  # Dodged bars per fold
  # scale_fill_viridis_d(option = "cividis") +  # Apply magma color palette
  scale_fill_manual(values = colorRampPalette(c("lightgoldenrod2", "deeppink4"))(10)) +
  # scale_fill_brewer(palette = "Blues") +
  scale_x_discrete(labels = c("(Intercept)" = "nPA : LC 1",  # Manually change fwi_norm to FWI
                              "lccs2" = "LC 2",  # lccs1 to LCCS1
                              "lccs3" = "LC 3",                              
                              "PA_status1:lccs2" = "PA : LC 2", 
                              "PA_status1:lccs3" = "PA : LC 3",
                              "PA_status1" = "PA",
                              "fwi_norm"="FWI")) +
  labs(x = "Variable", 
       y = "Coefficient Value",
       fill = "Fold") +
  theme_minimal()
ggsave(filename = paste0(path_fig, "GAM_Coef_all.png"), plot = last_plot(),
       width = 10, height = 6, dpi = 300)



# Assuming coef_df contains the original data
coef_mean_df <- coef_df %>%
  group_by(Variable) %>%
  summarise(Mean_Coefficient = mean(Coefficient, na.rm = TRUE))  # Calculate the mean coefficient per variable


ggplot(coef_mean_df, aes(x = Variable, y = Mean_Coefficient, fill = Variable)) +
  geom_bar(stat = "identity") +  # Plot bars with the mean coefficient value
  scale_fill_manual(values = colorRampPalette(c("lightgoldenrod2", "deeppink4"))(10)) +  # Apply color palette
  scale_x_discrete(labels = c("(Intercept)" = "nPA : LC 1",  # Manually change labels
                              "lccs2" = "LC 2", 
                              "lccs3" = "LC 3",                              
                              "PA_status1:lccs2" = "PA : LC 2", 
                              "PA_status1:lccs3" = "PA : LC 3",
                              "PA_status1" = "PA",
                              "fwi_norm" = "FWI")) +
  labs(x = "Variable", 
       y = "Mean Coefficient Value",
       fill = "Variable") +
  theme_minimal()+
  theme(legend.position = "None")
ggsave(filename = paste0(path_fig, "GAM_Coef_all_mean.png"), plot = last_plot(),
       width = 8, height = 5, dpi = 300)

################################################################################
## SIMPLE Random Forest

predictors <- c("lccs", "PA_status", "Latitude_norm", "Longitude_norm", "fwi_norm") #"Month", "Year", 
response <- "Burned_Area"


## split test and training so that they have the same number of LC classes, and that one 
# class is not under-represented


boxplot(dataTrain$Burned_Area ~ dataTrain$lccs, las=2)
boxplot(dataTrain$Burned_Area ~ dataTrain$PA_status, las=2)
boxplot(dataTrain$Burned_Area ~ dataTrain$Month, las=2)
boxplot(dataTrain$Burned_Area ~ dataTrain$Year, las=2)


data$interaction <- interaction(data$PA_status, data$lccs)


############## TUNE RF MODEL ###############

### optimal mtry value
fit1 <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude + Longitude, 
                     data = dataTrain, mtry = 1)
fit2 <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude + Longitude, 
                     data = dataTrain, mtry = 2)
fit3 <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude + Longitude, 
                     data = dataTrain, mtry = 3)

print(fit1)
print(fit2)
print(fit3)

#### mtry = 2


# Try different values of ntree (e.g., 100, 200, 500)
fit1 <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude + Longitude, 
                     data = dataTrain, ntree = 100)
fit2 <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude + Longitude, 
                     data = dataTrain, ntree = 200)
fit3 <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude + Longitude, 
                     data = dataTrain, ntree = 250)

# Compare performance (e.g., OOB error rate) of different models
print(fit1)
print(fit2)
print(fit3)

#### better performance with ntree = 200


# Try different values of nodesize (e.g., 1, 5, 10)
fit1 <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude + Longitude, 
                     data = dataTrain, nodesize = 1)
fit2 <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude + Longitude, 
                     data = dataTrain, nodesize = 5)
fit3 <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude + Longitude, 
                     data = dataTrain, nodesize = 10)

# Compare performance of different nodesize values
print(fit1)
print(fit2)
print(fit3)

#### best wiith nodesize = 5


mRF <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude + Longitude, 
                    data = dataTrain, nodesize = 5, ntree = 200, mtry=2)

plot(mRF)  # Plot OOB error for the final model



#####################
## RF per fold
importance_list <- list()
set.seed(84)
# Loop over each fold
for (fold in 1:10) {
  # Define train-test split for this fold
  train_data <- dataTrain %>% filter(!is.na(clusterfold) & clusterfold != fold) ## fold random doesn't make sense
  test_data <- dataTrain %>% filter(!is.na(clusterfold) & clusterfold == fold)
  
  # Train Random Forest on the training set
  modelRF <- randomForest(x = train_data[, predictors], 
                          y = train_data$Burned_Area, 
                          importance = TRUE, 
                          mtry=2,
                          ntree = 200,
                          nodesize=5)
  
  # Extract variable importance for this fold and store it
  fold_importance <- as.data.frame(importance(modelRF))
  # Add fold number as a new column
  fold_importance$Fold <- fold
  # store in list
  importance_list[[fold]] <- fold_importance
}

# Combine all fold importances into a single data frame
all_importance <- do.call(rbind, importance_list)




### Plotting the importances
all_importance$Variable <- rep(rownames(importance_list[[1]]), times = length(importance_list))
all_importance$Fold <- rep(1:10, each = nrow(importance_list[[1]]))

## Numbers
mean_incMSE <- all_importance %>%
  group_by(Variable) %>%
  summarise(mean_IncMSE = mean(`%IncMSE`, na.rm = TRUE),
            max_IncMSE = max(`%IncMSE`, na.rm = TRUE))

# Print the result
print(mean_incMSE)
  


# Define custom colors for each variable
var_cols <- c('fwi_norm' = 'deeppink4', 
              'Latitude_norm' = 'lightgoldenrod3', 
              'Longitude_norm' = 'goldenrod3', 
              'lccs' = 'olivedrab4', 
              'PA_status' = 'steelblue3')

# Convert variable column to a factor for proper ordering
all_importance$Variable <- factor(all_importance$Variable, levels = names(var_cols))

var_labels <- c('fwi_norm' = 'FWI', 
                'Latitude_norm' = 'Latitude', 
                'Longitude_norm' = 'Longitude', 
                'lccs' = 'LCC', 
                'PA_status' = 'PA Status')

# Create stacked histogram
ggplot(all_importance, aes(x = factor(Fold), y = `%IncMSE`, fill = Variable)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars
  scale_fill_manual(values = var_cols, labels = var_labels) +  # Apply custom colors
  theme_minimal() +
  labs(x = "Fold",
       y = "% Increase in MSE",
       fill = "Variable") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))
ggsave(filename = paste0(path_fig, "RF_VarImp_folds.png"), plot = last_plot(),
       width = 8, height = 6, dpi = 300)



########## MODEL COMMPARISON ##############

# Initialize an empty data frame to store performance results
performance_results <- data.frame(Fold = integer(),
                                  Model = character(),
                                  MSE = numeric(),
                                  R2 = numeric(),
                                  stringsAsFactors = FALSE)


# Loop over folds (1:10)
#  + fwi_norm + Latitude + Longitude
for (fold in 1:10) {
  
  # Subset train-test data
  train_data <- dataTrain %>% filter(clusterfold != fold)
  test_data <- dataTrain %>% filter(clusterfold == fold)
  
  # --- 1. Quantile Regression ---
  modelQR <- rq(Burned_Area ~ PA_status * lccs + fwi_norm, data = train_data, tau = 0.5)
  pred_QR <- predict(modelQR, test_data)
  mse_QR <- mean((test_data$Burned_Area - pred_QR)^2)
  r2_QR <- 1 - sum((test_data$Burned_Area - pred_QR)^2) / sum((test_data$Burned_Area - mean(test_data$Burned_Area))^2)
  
  # --- 2. Generalized Additive Model (GAM) ---
  modelGAM <- gam(Burned_Area ~ PA_status * lccs + fwi_norm + s(Latitude, bs='cr') + s(Longitude, bs='cr'), data = train_data, method="REML")
  pred_GAM <- predict(modelGAM, test_data)
  mse_GAM <- mean((test_data$Burned_Area - pred_GAM)^2)
  r2_GAM <- 1 - sum((test_data$Burned_Area - pred_GAM)^2) / sum((test_data$Burned_Area - mean(test_data$Burned_Area))^2)
  
  # --- 3. Random Forest ---
  modelRF <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude + Longitude, 
                          data = train_data, ntree = 200, mtry=2, nodesize=5)
  pred_RF <- predict(modelRF, test_data)
  mse_RF <- mean((test_data$Burned_Area - pred_RF)^2)
  r2_RF <- 1 - sum((test_data$Burned_Area - pred_RF)^2) / sum((test_data$Burned_Area - mean(test_data$Burned_Area))^2)
  
  # Store results (MSE and R2)
  performance_results <- rbind(performance_results,
                               data.frame(Fold = fold, Model = "Quantile Regression", MSE = mse_QR, R2 = r2_QR),
                               data.frame(Fold = fold, Model = "GAM", MSE = mse_GAM, R2 = r2_GAM),
                               data.frame(Fold = fold, Model = "Random Forest", MSE = mse_RF, R2 = r2_RF))
}


######
# Initialize a dataframe to store performance metrics
performance_results <- data.frame(Model = character(),
                                  MSE = numeric(),
                                  R2 = numeric(),
                                  stringsAsFactors = FALSE)

# --- 1. Quantile Regression (QR) ---
modelQR <- rq(Burned_Area ~ PA_status * lccs + fwi_norm, data = dataTrain, tau = 0.5)
pred_QR <- predict(modelQR, dataTest)
mse_QR <- mean((dataTest$Burned_Area - pred_QR)^2)
r2_QR <- 1 - sum((dataTest$Burned_Area - pred_QR)^2) / sum((dataTest$Burned_Area - mean(dataTest$Burned_Area))^2)

# --- 2. Generalized Additive Model (GAM) ---
modelGAM <- gam(Burned_Area ~ PA_status * lccs + fwi_norm + s(Latitude, bs='cr') + s(Longitude, bs='cr'), 
                data = dataTrain, method="REML")
pred_GAM <- predict(modelGAM, dataTest)
mse_GAM <- mean((dataTest$Burned_Area - pred_GAM)^2)
r2_GAM <- 1 - sum((dataTest$Burned_Area - pred_GAM)^2) / sum((dataTest$Burned_Area - mean(dataTest$Burned_Area))^2)

# --- 3. Random Forest (RF) ---
modelRF <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude + Longitude, 
                        data = dataTrain, ntree = 200, mtry=2, nodesize=5)
pred_RF <- predict(modelRF, dataTest)
rmse_RF <- sqrt(mean((dataTest$Burned_Area - pred_RF)^2))
r2_RF <- 1 - sum((dataTest$Burned_Area - pred_RF)^2) / sum((dataTest$Burned_Area - mean(dataTest$Burned_Area))^2)

# Store results
performance_results <- rbind(performance_results,
                             data.frame(Model = "Quantile Regression", MSE = mse_QR, R2 = r2_QR),
                             data.frame(Model = "GAM", MSE = mse_GAM, R2 = r2_GAM),
                             data.frame(Model = "Random Forest", MSE = mse_RF, R2 = r2_RF))

# Print performance results
print(performance_results)



performance_results$RMSE <- sqrt(performance_results$MSE)

anova_results <- aov(RMSE ~ Model, data = performance_results)
summary(anova_results)


kruskal.test(RMSE ~ Model, data = performance_results)


# Compute the mean R² for each model
mean_r2_results <- performance_results %>%
  group_by(Model) %>%
  summarize(mean_R2 = min(RMSE))

# Print the mean R² vamin()# Print the mean R² values
print(mean_r2_results)

# Define custom colors for each variable
model_cols <- c('Quantile Regression' = 'hotpink4', 
                'GAM' = 'darkslategray4', 
                'Random Forest' = 'goldenrod2')

# Convert variable column to a factor for proper ordering
performance_results$Model <- factor(performance_results$Model, levels = names(model_cols))



## with median value written above line
ggplot(performance_results, aes(x = Model, y = RMSE, fill = Model)) +
  geom_boxplot() +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = round(..y.., 2)), 
    vjust = 1.5, 
    color = "black", 
    size = 6) +
  scale_fill_manual(values = model_cols) +
  scale_x_discrete(labels = c("Quantile Regression" = "QR", "GAM" = "GAM", "Random Forest" = "RF")) +
  
  theme_minimal() +
  labs(y = "RMSE", 
       x = "Model") + 
  theme(legend.position = "none",
        axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))

ggsave(filename = paste0(path_fig, "Anova_allPred_boxplot_RMSE.png"), plot = last_plot(),
       width = 8, height = 6, dpi = 300)


ggplot(performance_results, aes(x = Model, y = R2, fill = Model)) +
  geom_boxplot() +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = round(..y.., 3)), 
    vjust = -0.45, 
    color = "black", 
    size = 6) +
  scale_fill_manual(values = model_cols) +
  scale_x_discrete(labels = c("Quantile Regression" = "QR", "GAM" = "GAM", "Random Forest" = "RF")) +
  theme_minimal() +
  labs(y = "R2", 
       x = "Model") + 
  theme(legend.position = "none",
        axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))
ggsave(filename = paste0(path_fig, "Anova_allPred_boxplot_R2.png"), plot = last_plot(),
       width = 8, height = 6, dpi = 300)




mean(performance_results$RMSE)

mean_rmse_per_model <- performance_results %>%
  group_by(Model) %>%
  summarise(Mean_RMSE = min(R2, na.rm = TRUE))

print(mean_rmse_per_model)


#### PACKAGE compare performance
performance_list <- list()
for (fold in 1:10) {
  
  # Subset train-test data
  train_data <- dataTrain %>% filter(clusterfold != fold)
  test_data <- dataTrain %>% filter(clusterfold == fold)
  
  # QR models
  qr0 <- rq(Burned_Area ~ PA_status, data = train_data, tau = 0.5)
  qr1 <- rq(Burned_Area ~ PA_status * lccs, data = train_data, tau = 0.5)
  qr2 <- rq(Burned_Area ~ PA_status * lccs + fwi_norm, data = train_data, tau = 0.5)
  # qr3 <- rq(Burned_Area ~ PA_status * lccs + fwi_norm + Latitude + Longitude, data = train_data, tau = 0.5)
  # qr4 <- rq(Burned_Area ~ PA_status * lccs + fwi_norm, data = train_data, tau = 0.25)
  # qr5 <- rq(Burned_Area ~ PA_status * lccs + fwi_norm, data = train_data, tau = 0.75)
  
  # GAM models
  m0 <- gam(Burned_Area ~ PA_status, data = train_data)
  m1 <- gam(Burned_Area ~ PA_status + s(Latitude_norm, bs="cr") + s(Longitude_norm, bs="cr"), data = train_data)
  m2 <- gam(Burned_Area ~ PA_status * lccs + s(Latitude_norm, bs="cr") + s(Longitude_norm, bs="cr"), data = train_data, method="REML")
  m3 <- gam(Burned_Area ~ PA_status * lccs + fwi_norm + s(Latitude_norm, bs="cr") + s(Longitude_norm, bs="cr"), data = train_data, method="REML")
  
  
  # Compare model performance
  perform_res <- compare_performance(qr0, qr1, qr2, m0, m1, m2, m3)
  
  # Store results
  performance_list[[fold]] <- as.data.frame(perform_res) %>% mutate(Fold = fold)
}


# Combine all folds into a single dataframe
performance_df <- bind_rows(performance_list)

# Compute the mean performance metrics for each model
performance_summary <- performance_df %>%
  group_by(Name) %>%
  summarise(
    AIC = mean(AIC, na.rm = TRUE),
    BIC = mean(BIC, na.rm = TRUE),
    R2 = mean(R2, na.rm = TRUE),
    RMSE = mean(RMSE, na.rm = TRUE),
    .groups = "drop"
  )

# Print the mean performance summary
print(performance_summary)



### on dataTrain
# QR
qr0 <- rq(Burned_Area ~ PA_status, data = dataTrain, tau = 0.5)
qr1 <- rq(Burned_Area ~ PA_status * lccs, data = dataTrain, tau = 0.5)
qr2 <- rq(Burned_Area ~ PA_status * lccs + fwi_norm, data = dataTrain, tau = 0.5)
# qr3 <- rq(Burned_Area ~ PA_status * lccs + fwi_norm + Latitude + Longitude, data = dataTrain, tau = 0.5)

compare_performance(qr0, qr1, qr2, qr3)

# GAM
m0 <- gam(Burned_Area ~ PA_status, data = dataTrain)
m1 <- gam(Burned_Area ~ PA_status + s(Latitude, bs="cr") + s(Longitude, bs="cr"), data = dataTrain)
m2 <- gam(Burned_Area ~ PA_status * lccs + s(Latitude, bs="cr") + s(Longitude, bs="cr"), data = dataTrain, method="REML")
m3 <- gam(Burned_Area ~ PA_status * lccs + fwi_norm + s(Latitude, bs="cr") + s(Longitude, bs="cr"), data = dataTrain, method="REML")

compare_performance(m0, m1, m2, m3)

summary(m3)

## compare all QR and GAM with dataTrain
compare_performance(qr0, qr1, qr2, m0, m1, m2, m3)

### on dataTest
# QR
qr0 <- rq(Burned_Area ~ PA_status, data = dataTest, tau = 0.5)
qr1 <- rq(Burned_Area ~ PA_status * lccs, data = dataTest, tau = 0.5)
qr2 <- rq(Burned_Area ~ PA_status * lccs + fwi_norm, data = dataTest, tau = 0.5)
# qr3 <- rq(Burned_Area ~ PA_status * lccs + fwi_norm + Latitude + Longitude, data = dataTrain, tau = 0.5)

compare_performance(qr0, qr1, qr2, qr3)

# GAM
m0 <- gam(Burned_Area ~ PA_status, data = dataTest)
m1 <- gam(Burned_Area ~ PA_status + s(Latitude, bs="cr") + s(Longitude, bs="cr"), data = dataTest)
m2 <- gam(Burned_Area ~ PA_status * lccs + s(Latitude, bs="cr") + s(Longitude, bs="cr"), data = dataTest, method="REML")
m3 <- gam(Burned_Area ~ PA_status * lccs + fwi_norm + s(Latitude_norm, bs="cr") + s(Longitude_norm, bs="cr"), data = dataTest, method="REML")

compare_performance(m0, m1, m2, m3)

summary(m3)

compare_performance(qr0, qr1, qr2, m0, m1, m2, m3)


rf1 <- randomForest(Burned_Area ~ PA_status + lccs + fwi_norm + Latitude_norm + Longitude_norm, 
                        data = dataTest, ntree = 200, mtry=2, nodesize=5)

compare_performance(qr2, m3, rf1)

############# NORTH SOUTH ###############

train_south <- dataTrain %>%
  filter(clusterfold %in% c(2, 4, 5))

train_north <- dataTrain %>%
  filter(!(clusterfold %in% c(2, 4, 5)))

data_combined <- bind_rows(
  train_south %>% mutate(Region = "South"),
  train_north %>% mutate(Region = "North")
)

median(train_south$Burned_Area[train_south$PA_status == 1])
median(train_south$Burned_Area[train_south$PA_status == 0])

median(train_north$Burned_Area[train_north$PA_status == 1])
median(train_north$Burned_Area[train_north$PA_status == 0])

# Create the boxplot
ggplot(data_combined, aes(x = Region, y = Burned_Area, fill = PA_status)) +
  geom_boxplot(width = 0.4) +  # Narrower boxes
  scale_fill_manual(
    values = c("0" = "lightblue", "1" = "darkcyan"),
    name = "Protection Status",
    labels = c("Non-Protected", "Protected")) +
  labs(x = "Region", y = "Burned Area [%]") +
  theme_minimal() +
  # Add median line in orange and text next to it
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 1), group = PA_status),
               color = "purple4", size = 5, fontface = "bold", vjust = -0.5, hjust = 0.5) +  # Position text near median
  theme(axis.text.x = element_text(size = 14, hjust = 0.5))+
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))

ggsave(filename = paste0(path_fig, "Boxplot_BA_PA_NorthSouth.png"), plot = last_plot(),
       width = 8, height = 4, dpi = 300)

# different version
ggplot(data_combined, aes(x = as.factor(PA_status), y = Burned_Area, fill = Region)) +
  geom_boxplot(width = 0.3) +  # Narrower boxes
  scale_fill_manual(
    values = c("South" = "plum", "North" = "steelblue"),
    name = "Region",
    labels = c("North", "South")) +
  scale_x_discrete(labels = c("0" = "Non-Protected", "1" = "Protected")) +  # Label x-axis for PA_status
  labs(x = "Protection Status", y = "Burned Area [%]") +
  theme_minimal() +
  # Add median line in orange and text next to it
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 1), group = PA_status),
               color = "maroon4", size = 5, fontface = "bold", vjust = -0.5, hjust = 0.5) +  # Position text near median
  theme(axis.text.x = element_text(size = 14, hjust = 1))


model_south <- gam(Burned_Area ~ PA_status * lccs + fwi_norm + s(Latitude_norm, bs="cr") + s(Longitude_norm, bs="cr"), data = train_south, method="REML")
summary(model_south)


model_north <- gam(Burned_Area ~ PA_status * lccs + fwi_norm + s(Latitude_norm, bs="cr") + s(Longitude_norm, bs="cr"), data = train_north, method="REML")
summary(model_north)


### Predict
test_south <- dataTest %>%
  filter(clusterfold %in% c(2, 4, 5))

test_north <- dataTest %>%
  filter(!(clusterfold %in% c(2, 4, 5)))

predictSouth <- predict(model_south, test_south)
rmse_south <- sqrt(mean((test_south$Burned_Area - predictSouth)^2))
r2_south <- 1 - sum((test_south$Burned_Area - predictSouth)^2) / sum((test_south$Burned_Area - mean(test_south$Burned_Area))^2)

predictNorth <- predict(model_north, test_north)
rmse_north <- sqrt(mean((test_north$Burned_Area - predictNorth)^2))
r2_north <- 1 - sum((test_north$Burned_Area - predictNorth)^2) / sum((test_north$Burned_Area - mean(test_north$Burned_Area))^2)


########################## PA LAG BA #############################

dataPA <- dataTrain %>%
  filter(PA_status == 1) %>%
  mutate(Lag =  Year - PA_year)

dataPA <- dataTrain %>%
  filter(PA_status == 1) %>%  # Filter for PA_status == 1
  mutate(Lag = Year - PA_year) %>%  # Compute Lag
  group_by(Latitude, Longitude) %>%  # Group by unique location
  slice_min(Year, with_ties = FALSE) %>%  # Keep only the row with the lowest Year per location
  ungroup() 

summary(dataPA$Lag)
quantile(dataPA$Lag, probs = seq(0, 1, 0.25))

# dataPA <- dataPA %>%
  # mutate(Lag_bin = cut(Lag, 
    #                    breaks = c(-Inf, 5, 13, 31, Inf), 
      #                  labels = c("0-5", "6-13", "14-31", ">31")))

dataPA <- dataPA %>%
  mutate(Lag_bin = cut(Lag, 
                       breaks = c(-Inf, 3, 11, 29, Inf), 
                       labels = c("0-3", "4-11", "12-29", ">29")))

# breaks = c(-Inf, 4, 12, 30, Inf), 
# labels = c("0-4", "5-12", "13-30", ">30")))

# Check the distribution of 'Lag_binned'
table(dataPA$Lag_bin)

### LAG AND SPATIAL
mLag <- gam(Burned_Area ~ Lag_bin + s(Latitude_norm, bs="cr") + s(Longitude_norm, bs="cr"), 
            data = dataPA, method="REML")

# Median outcome:
# mLag <- gam(Burned_Area ~ Lag_bin, data = dataPA)

summary(mLag)

### ALL PREDICTORS
mLag <- gam(Burned_Area ~ Lag_bin * lccs + fwi_norm + s(Latitude, bs="cr") + s(Longitude, bs="cr"), 
            data = dataPA, method="REML")

summary(mLag)


aggregate(Burned_Area ~ Lag_bin, data = dataPA, FUN = median)
table(dataPA$Lag_bin)




############### PLOTS

ggplot(dataPA, aes(x = Lag)) +
  geom_histogram(binwidth = 5, color = "black", fill="cadetblue4") +
  labs(x = "Lag [Years]", y = "Count") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, max(dataPA$Lag), by = 10)) +
  theme(legend.position = c(0, 1),        # Position the legend in the top left corner
        legend.justification = c(0, 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))

ggsave(filename = paste0(path_fig, "Histogram_PA_Lag.png"), plot = last_plot(),
       width = 6, height = 4, dpi = 300)


lag_cols <- c("darkcyan", "hotpink2", "darkgoldenrod2", "yellowgreen")

ggplot(dataPA, aes(x = Lag_bin, y = Burned_Area)) +
  geom_boxplot(width = 0.3, fill=lag_cols) +  # Narrower boxes
  labs(x = "Lag (years)", y = "Burned Area [%]") +
  theme_minimal() +
  # Add median line in orange and text next to it
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 1)),
               color = "black", size = 5, fontface = "bold", vjust = -0.5, hjust = 0.5) +  # Position text near median
  theme(axis.text.x = element_text(size = 14, hjust = 1))+
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))
ggsave(filename = paste0(path_fig, "Boxplot_BA_PA_LagBin.png"), plot = last_plot(),
       width = 6, height = 4, dpi = 300)



dataPA_sf <- st_as_sf(dataPA, coords = c("Longitude", "Latitude"), crs=4326)
ggplot() +
  geom_sf(data = pt_sf, fill = NA, color = "black", size = 0.8) +
  # Plot individual points
  geom_sf(data = dataPA_sf, aes(color = as.factor(Lag_bin)), size = 1) +
  
  # Add map elements
  #annotation_scale(location = "br") +
  annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering) +
  
  scale_color_manual(values = lag_cols, name = "Lag Group") +
  
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))+
  # Explicitly set longitude breaks for 9°W, 8°W, and 7°W
  scale_x_continuous(breaks = c(-9, -8, -7), labels = c("9°W", "8°W", "7°W"))

ggsave(filename = paste0(path_fig, "Map_PT_Lag.png"), plot = last_plot(),
       width = 5, height = 8, dpi = 300)


########### GINI #############

data_gini <- dataTrain %>%
  filter(Gini > 0.75)

data_gini <- data_gini %>%
  group_by(PA_lccs) %>% 
  slice_sample(n = min(table(data_gini$PA_lccs))) %>%  # Ensures equal samples per PA_lccs
  ungroup()

table(data_gini$lccs)

mGini <- gam(Burned_Area ~ PA_status * lccs + fwi_norm + s(Latitude_norm, bs="cr") + s(Longitude_norm, bs="cr"), 
             data = data_gini, method="REML")

summary(mGini)


########### OTHER STATS ##############

# Create a summary of Burned_Area by Region
region_summary <- dataPA %>%
  group_by(cluster) %>%
  summarise(median_BA = median(Burned_Area), mean_BA = mean(Burned_Area), .groups = 'drop')


### BA across Clusters
ggplot(data, aes(x = factor(cluster), y = Burned_Area)) + 
  geom_boxplot(fill = fold_colors, color = "black") +  # Customize the boxplot color
  labs(x = "Cluster", y = "Burned Area [%]") + 
  theme_minimal() + 
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 1)), color = "cyan2", size = 4.5, vjust = -0.5) +
  theme(axis.text.x = element_text(hjust = 1))+
  theme(axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))
ggsave(filename = paste0(path_fig, "Boxplot_BA_per_Cluster.png"), plot = last_plot(),
       width = 8, height = 6, dpi = 300)


### BA across Clusters
ggplot(data, aes(x = factor(Month), y = Burned_Area)) + 
  geom_boxplot(fill = "lightblue", color = "black") +  # Customize the boxplot color
  labs(x = "Month", y = "Burned Area [%]") + 
  theme_minimal() + 
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 1)), color = "lightblue1", size = 4, vjust = -0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(path_fig, "Boxplot_BA_per_Cluster.png"), plot = last_plot(),
       width = 8, height = 6, dpi = 300)



lc_cols <- c("orange", "maroon", "darkcyan")
ggplot(data_gini, aes(x = as.factor(lccs), y = Burned_Area, fill = as.factor(lccs), alpha = as.factor(PA_status))) +
  geom_boxplot(position = position_dodge(width = 0.75)) +  # Dodge for separation
  scale_fill_manual(values = forest_colors) + 
  scale_alpha_manual(values = c("0" = 0.6, "1" = 1)) +  # Set transparency for PA_status
  labs(x = "Land Cover Class (LCC)", y = "Burned Area [%]", fill = "LCCS", alpha = "PA Status") +
  theme_minimal()+
  theme(legend.position = "None") +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = round(..y.., 1), group = interaction(lccs, PA_status)),  # Separate medians by PA_status
    position = position_dodge(width = 0.75),  # Align within boxes
    color = "black", 
    size = 6, 
    vjust = -0.3  # Adjust height to stay within the box
  ) +
  theme(axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title larger
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),  # X-axis labels larger
        axis.text.y = element_text(size = 14))

ggsave(filename = paste0(path_fig, "Boxplot_BA_per_LCC.png"), plot = last_plot(),
       width = 8, height = 6, dpi = 300)




### BA map 
ggplot() +
  # Plot the spatial polygons with black borders
  geom_sf(data = pt_sf, fill = NA, color = "black", size = 0.8) +
  
  # Plot individual points with color mapped to Burned_Area
  geom_sf(data = data_sf, aes(color = Burned_Area), size = 0.8) +
  
  # Apply a continuous color scale to Burned_Area
  scale_color_gradient(low = "gold", high = "darkred", name = "BA [%]") +
  
  # Add a north arrow to the map
  annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering) +
  
  # Add labels for the axes
  labs(x = "Longitude", y = "Latitude") +
  
  # Apply minimal theme
  theme_minimal()


### Fire frequency
library(viridis)
grid_size <- 0.05  # Adjust resolution (smaller = finer grid)
bbox <- st_bbox(data_sf)

# Create grid over bounding box
grid <- st_make_grid(st_as_sfc(bbox), cellsize = grid_size, what = "polygons")
grid_sf <- st_sf(geometry = grid)

# Count fire occurrences per grid cell
fire_counts <- st_join(grid_sf, data_sf) %>%
  group_by(geometry) %>%
  summarise(Fire_Count = n(), .groups = "drop")  # Count occurrences

# Set missing fire count values to 0
fire_counts$Fire_Count[is.na(fire_counts$Fire_Count)] <- 0

# Plot the fire frequency map
ggplot() +
  # Add fire frequency grid
  geom_sf(data = fire_counts, aes(fill = Fire_Count), color = NA) +
  # Set color scale: 0 is white, higher values use "inferno" colors
  # scale_fill_gradientn(colors = c("white", viridis(5, option = "magma")), name = "Fire Frequency") +
  scale_fill_gradientn(colors = c("white", "lightgreen", "orange", "red", "darkred"), name = "Fire Frequency") +
  # Add Portugal shape outline
  geom_sf(data = pt_sf, fill = NA, color = "black", size = 0.8) +
  labs(title = "Fire Frequency Map of Portugal", x = "Longitude", y = "Latitude") +
  theme_minimal()

hist(log(fire_counts$Fire_Count))

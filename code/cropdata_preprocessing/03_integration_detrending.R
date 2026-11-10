library(ncdf4)        # For reading NetCDF files
library(raster)       # For working with raster data (e.g., NetCDF via brick/stack)
library(tidyverse)    # For data manipulation (select, mutate, filter, etc.)

# Construct the paths
repo_path <- here::here()
tgrid_path   <- file.path(repo_path, "data", "processed", "other", "tgrid_ISIMIP_ctrs.RData")
lonlat_path  <- file.path(repo_path, "data", "processed", "other", "lonlat_ISIMIP_ctrs.RData")
LU_list_path <- file.path(repo_path, "data", "processed", "other", "LU_list.RData")

# Load the files
load(tgrid_path)
load(lonlat_path)
load(LU_list_path)

# Post-load operation
tgrid_ctr <- sort(unique(tgrid$ctr))

detrend_yield <- function(yield, year){
  dat <- tibble("log_yield" = log(yield), "year" = year)
  fitted <- predict(lm(log_yield ~ log(year), dat))
  # prevent non-positive trend values
  fitted <- pmax(exp(fitted), 0.001)
  return(fitted)
}

crops <- c("mai", "ri1", "ri2", "soy", "swh", "wwh")
models <- c("acea", "ldndc", "lpjml", "pdssat", "pepic", "promet", "simplace-lintul5", "isam", "epic-iiasa", "cygma1p74", "crover", "dssat-pythia", "lpj-guess")

for (crop in crops) {  # Loop through each crop
  
  # Create another name based on the crop for accessing the LU_list later 
  if (crop == "ri1" | crop == "ri2") {
    LU_crop <- "rice"
  } else if (crop == "wwh" | crop == "swh") {
    LU_crop <- "temperate_cereals"
  } else if (crop == "soy") {
    LU_crop <- "oil_crops_soybean"
  } else {
    LU_crop <- "maize"
  }
  
  irrigated_LU <- LU_list[[paste("LU_", LU_crop, "_irrigated", sep = "")]] %>% 
    mutate(irr_area =  `2015`) %>% 
    dplyr::select(irr_area, lon, lat)
  rainfed_LU <- LU_list[[paste("LU_", LU_crop, "_rainfed", sep = "")]]  %>% 
    mutate(rain_area =  `2015`) %>% 
    dplyr::select(lon, lat, rain_area)
  
  # Specify the directory where your NetCDF files are located
  directory <- file.path(
    repo_path,
    "data", "raw", "benchmark_yields",
    crop
  )
  
  # Use list.files() to find files matching the specified pattern
  file_paths <- list.files(directory, pattern = "yield_\\d{4}\\.nc4", full.names = TRUE)
  
  # Extract years from file names
  years <- unique(as.numeric(sub(".*/yield_(\\d{4})\\.nc4", "\\1", file_paths)))
  
  # Read NetCDF files into a list
  nc_list <- lapply(file_paths, raster::brick)
  
  # Stack RasterBricks into a single RasterStack
  stacked_data <- stack(nc_list)
  
  # Extract lon and lat coordinates
  lonlat <- raster::xyFromCell(stacked_data, 1:ncell(stacked_data))
  lon <- lonlat[, 1]
  lat <- lonlat[, 2]
  
  # Convert to Dataframe
  gdhy <- data.frame(lon = lon, lat = lat, values(stacked_data))
  
  # Rename the 'layer' variables to correspond to the years
  names(gdhy)[-c(1:2)] <- years
  
  # Convert to long format
  gdhy <- pivot_longer(gdhy, cols = -c(lon, lat), names_to = "year", values_to = "yield") %>% 
    filter(!is.na(yield))
  
  # Change longitudes larger than 180 to negative
  gdhy$lon[gdhy$lon > 180] <- gdhy$lon[gdhy$lon > 180] - 360
  
  # Create an empty list to store datasets for each model
  crop_datasets <- list()
  
  for (model in models) {  # Loop through each model available for the crop
    
    # Get all runs 
    
    tryCatch({
      # Construct the path to the .RData file
      load_path <- file.path(
        repo_path,
        "data", "processed", "GGCMI_dataframes",
        crop,
        paste0(model, "_", crop, ".RData")
      )
      
      # Load the file
      load(load_path)
      
      
      years <- 1901:2016
      
      # First irrigated
      irr <- merge(yld_long_irr, tgrid, by = c("lon", "lat"))  
      
      
      # Rename the 'v' variables to correspond to the years
      names(irr)[grepl("^V", names(irr))] <- years
      
      # Convert to long format
      irr <- pivot_longer(irr, cols = -c(ctr, cellarea, lon, lat), names_to = "year", values_to = "irr") %>% 
        dplyr::select(ctr, lon, lat, year, irr)
      
      # Then rainfed 
      rain <- merge(yld_long_rain, tgrid, by = c("lon", "lat")) 
      
      # Rename the 'v' variables to correspond to the years
      names(rain)[grepl("^V", names(rain))] <- years
      
      # Convert to long format
      rain <- pivot_longer(rain, cols = -c(ctr, cellarea, lon, lat), names_to = "year", values_to = "rain") %>% 
        dplyr::select(cellarea, lon, lat, year, rain)
      
      
      # Join simulation (irrigated and rainfed) with GDHY benchmark yield data
      dat <- merge(gdhy, rain, by = c("lon", "lat", "year")) %>%
        merge(irr, by = c("lon", "lat", "year")) %>%
        mutate(crop = crop,
               model = model)
      
      # Store the dataset in the list
      crop_datasets[[model]] <- dat
      
    }, error = function(e) {
      # Print a message if an error occurs (e.g., file not found)
      cat("No combination:", model, "and", crop, "\n")
      # Print the error message
      cat("Error message:", conditionMessage(e), "\n")
      
    })  # end of the trycatch-error
    
  } #end of the model loop
  
  # Combine datasets for all models into one dataset for this crop
  combined_dat <- do.call(rbind, crop_datasets) %>% 
    merge(irrigated_LU, by = c("lon", "lat")) %>% 
    merge(rainfed_LU, by = c("lon", "lat"))
  
  combined_obs <- combined_dat %>%
    dplyr::select(yield, lon, lat, year, ctr, rain_area, irr_area) %>% 
    distinct() %>% 
    mutate(year = as.numeric(year)) %>%
    filter(yield != 0 ) %>% 
    group_by(lon, lat) %>%
    drop_na %>% 
    filter(n() > 4) %>%   # Filter out groups with fewer than 5 non-zero observations
    mutate(
      fitted = detrend_yield(yield, year),
      divtrend_obs = yield / fitted, # Relative detrending
      difftrend_obs = yield - fitted  # Absolute detrending
    )
  
  combined_sim <- combined_dat %>%
    dplyr::select(irr, rain, lon, lat, year, model, ctr, rain_area, irr_area) %>% 
    mutate(year = as.numeric(year),
           tot_yld = (irr * irr_area + rain * rain_area) / (irr_area + rain_area)) %>% #Integrate rainfed and irrigated simulations using the corresponding crop areas from the land-use data
    filter(tot_yld != 0 ) %>% 
    group_by(lon, lat, model) %>%
    drop_na() %>% 
    filter(n() > 4) %>%   # Filter out groups with fewer than 5 non-zero observations
    mutate(
      fitted = detrend_yield(tot_yld, year),
      divtrend_sim = tot_yld / fitted, # Relative detrending
      difftrend_sim = tot_yld - fitted # Absolute detrending
    )
  
  #Join benchmark and simulated data
  combined <- combined_sim %>%
    left_join(combined_obs, by = c("year", "lon", "lat", "ctr", "irr_area", "rain_area"))
  
  # Dynamically assign the correct crop name to the combined dataset 
  assign(paste0("dat_", crop), combined)
}

# Construct the path to the output file
save_path <- file.path(
  repo_path,
  "data", "processed", "integrated_cropdata",
  "crop_specific_data.RData"
)

# Save the data
save(dat_mai, dat_ri1, dat_ri2, dat_swh, dat_wwh, dat_soy, file = save_path)

crops <- c("mai", "ri1", "ri2", "soy", "swh", "wwh")
models <- c("acea", "ldndc", "lpjml", "pdssat", "pepic", "promet", "simplace-lintul5", "isam", "epic-iiasa", "cygma1p74", "crover", "dssat-pythia", "lpj-guess")

# Initialize empty lists to store results for each crop
crop_productions <- list()

for (crop in crops) {  # Loop through each crop
  
  # Create another name based on the crop for accessing the LU_list later 
  if (crop == "ri1" | crop == "ri2") {
    LU_crop <- "rice"
  } else if (crop == "wwh" | crop == "swh") {
    LU_crop <- "temperate_cereals"
  } else if (crop == "soy") {
    LU_crop <- "oil_crops_soybean"
  } else {
    LU_crop <- "maize"
  }
  
  irrigated_LU <- LU_list[[paste("LU_", LU_crop, "_irrigated", sep = "")]] %>% 
    mutate(irr_area =  `2015`) %>% 
    dplyr::select(irr_area, lon, lat)
  rainfed_LU <- LU_list[[paste("LU_", LU_crop, "_rainfed", sep = "")]]  %>% 
    mutate(rain_area =  `2015`) %>% 
    dplyr::select(lon, lat, rain_area)
  
  # Specify the directory where your NetCDF files are located
  directory <- file.path(
    repo_path,
    "data", "raw", "benchmark_yields",
    crop
  )
  
  # Use list.files() to find files matching the specified pattern
  file_paths <- list.files(directory, pattern = "yield_\\d{4}\\.nc4", full.names = TRUE)
  
  # Read NetCDF files into a list
  nc_list <- lapply(file_paths, raster::brick)
  
  # Stack RasterBricks into a single RasterStack
  stacked_data <- stack(nc_list)
  
  # Extract lon and lat coordinates
  lonlat <- raster::xyFromCell(stacked_data, 1:ncell(stacked_data))
  lon <- lonlat[, 1]
  lat <- lonlat[, 2]
  
  # Convert to Data Frame
  gdhy <- data.frame(lon = lon, lat = lat, values(stacked_data))
  
  # Rename the 'layer' variables to correspond to the years
  years <- unique(as.numeric(sub(".*/yield_(\\d{4})\\.nc4", "\\1", file_paths)))
  names(gdhy)[-c(1:2)] <- years
  
  # Convert to long format
  gdhy <- pivot_longer(gdhy, cols = -c(lon, lat), names_to = "year", values_to = "yield") %>% 
    filter(!is.na(yield))
  
  # Change longitudes larger than 180 to negative
  gdhy$lon[gdhy$lon > 180] <- gdhy$lon[gdhy$lon > 180] - 360
  
  # Create an empty list to store datasets for each model
  crop_datasets <- list()
  
  for (model in models) {  # Loop through each model available for the crop
    tryCatch({
      # Construct the path to the .RData file
      load_path <- file.path(
        repo_path,
        "data", "processed", "GGCMI_dataframes",
        crop,
        paste0(model, "_", crop, ".RData")
      )
      
      # Load the file
      load(load_path)
      
      years <- 1901:2016
      
      # First irrigated
      irr <- merge(yld_long_irr, tgrid, by = c("lon", "lat"))  
      names(irr)[grepl("^V", names(irr))] <- years
      irr <- pivot_longer(irr, cols = -c(ctr, cellarea, lon, lat), names_to = "year", values_to = "irr") %>% 
        dplyr::select(ctr, lon, lat, year, irr)
      
      # then rainfed
      rain <- merge(yld_long_rain, tgrid, by = c("lon", "lat")) 
      names(rain)[grepl("^V", names(rain))] <- years
      rain <- pivot_longer(rain, cols = -c(ctr, cellarea, lon, lat), names_to = "year", values_to = "rain") %>% 
        dplyr::select(cellarea, lon, lat, year, rain)
      
      # merge the datasets
      dat <- merge(gdhy, rain, by = c("lon", "lat", "year")) %>%
        merge(irr, by = c("lon", "lat", "year")) %>%
        mutate(crop = crop, model = model)
      
      crop_datasets[[model]] <- dat
      
    }, error = function(e) {
      cat("No combination:", model, "and", crop, "\n")
      cat("Error message:", conditionMessage(e), "\n")
    })
  } # end of the model loop
  
  # Combine datasets for all models into one dataset for this crop
  combined_dat <- do.call(rbind, crop_datasets) %>% 
    merge(irrigated_LU, by = c("lon", "lat")) %>% 
    merge(rainfed_LU, by = c("lon", "lat"))
  
  # Calculate production for this crop for each model
  total_crop_production <- combined_dat %>%
    mutate(production = (irr * irr_area + rain * rain_area),
           area = irr_area + rain_area) 
  
  # Store the production for each crop
  crop_productions[[crop]] <- total_crop_production
}

# Combine production for all crops
final_total_yields <- do.call(rbind, crop_productions) %>%
  group_by(lon, lat, year, model) %>%
  summarise(total_production_sim = sum(production, na.rm = TRUE),
            total_area = sum(area, na.rm = TRUE),
            total_production_obs = sum(yield * area),
            ctr = first(ctr)) %>%
  mutate(aggregated_yield_sim = total_production_sim / total_area,
         aggregated_yield_obs = total_production_obs / total_area)

# Detrend observed yields
detr_obs <- final_total_yields %>%
  dplyr::select(aggregated_yield_obs, lon, lat, year, ctr, total_area) %>% 
  distinct() %>% 
  mutate(year = as.numeric(year)) %>%
  filter(aggregated_yield_obs != 0 ) %>% 
  group_by(lon, lat) %>%
  drop_na() %>% 
  filter(n() > 4) %>%   # Filter out groups with fewer than 5 non-zero observations
  mutate(
    fitted = detrend_yield(aggregated_yield_obs, year),
    divtrend_obs = aggregated_yield_obs / fitted, # Relative detrending
    difftrend_obs = aggregated_yield_obs - fitted, # absolute detrending
  )

# Detrend simulated yields
detr_sim <- final_total_yields %>%
  dplyr::select(aggregated_yield_sim, lon, lat, year, model, ctr, total_area) %>% 
  mutate(year = as.numeric(year)) %>%
  filter(aggregated_yield_sim != 0 ) %>% 
  group_by(lon, lat, model) %>%
  drop_na() %>% 
  filter(n() > 4) %>%   # Filter out groups with fewer than 5 non-zero observations
  mutate(
    fitted = detrend_yield(aggregated_yield_sim, year),
    divtrend_sim = aggregated_yield_sim / fitted, # Relative detrending
    difftrend_sim = aggregated_yield_sim - fitted, # absolute detrending
  )

# Save the results
# Construct output paths
detr_sim_path <- file.path(
  repo_path,
  "data", "processed", "integrated_cropdata",
  "aggr_sim.RData"
)
detr_obs_path <- file.path(
  repo_path,
  "data", "processed", "integrated_cropdata",
  "aggr_bench.RData"
)

# Save the objects
save(detr_sim, file = detr_sim_path)
save(detr_obs, file = detr_obs_path)

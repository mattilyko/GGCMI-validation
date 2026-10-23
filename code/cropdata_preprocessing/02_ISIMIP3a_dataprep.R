library(ncdf4)
library(stringr)

# Base path to the repository
repo_path <- here::here()

# Set working directory to the repo path (optional)
# setwd(repo_path)

# Construct the full path to the NetCDF file
country_mask_path <- file.path(repo_path, "data", "raw", "other", "countrymasks.nc")

# Open the NetCDF file
ctr.bin <- ncdf4::nc_open(filename = country_mask_path)

# Optional: check contents
# print(ctr.bin)  226 variables, for each country one variable

var.names=names(ctr.bin$var)
var.iso=gsub(x=var.names, pattern="m_", replacement="", fixed=TRUE)

# Create the matrix
lon=ncvar_get(ctr.bin,"lon")
lat=ncvar_get(ctr.bin,"lat")
lonlat <- as.matrix(expand.grid(lon,lat))

# Add the countries 
lonlat=cbind(lonlat, matrix(NA, nrow=nrow(lonlat), ncol=length(var.names)))
colnames(lonlat)=c("lon", "lat", var.iso)

for(i in 1:length(var.names)){
  xdata=ncvar_get(ctr.bin, var.names[i])
  xdata_long=matrix(xdata, prod(dim(xdata)[1:2]))
  lonlat[, 2+i]=xdata_long
}

lonlat=as.data.frame(lonlat)

country_cols <- names(lonlat)[3:216]
# Add country ID column
lonlat$ctr <- apply(lonlat[, country_cols], 1, function(row) {
  num_countries <- sum(row == 1)
  if (num_countries == 0) {
    return("0") # No country assignment 
  } else if (num_countries == 1) {
    return(country_cols[which(row == 1)]) # 1 unique country assignment
  } else {
    return(NA_character_) # cell is part of two countries 
  }
})

lonlat$ctr[68044]="AND"		
lonlat$ctr[85481]=0				# between China + India. We omit this cell. 
lonlat$ctr[93808]="BHS"		
lonlat$ctr[94525]="BHS"
lonlat$ctr[96692]="BHS"
lonlat$ctr[103144]="BLZ"
lonlat$ctr[103863]="BLZ"
lonlat$ctr[103864]="BLZ"
lonlat$ctr[103911]="VIR"
lonlat$ctr[104583]="BLZ"
lonlat$ctr[105303]="BLZ"
lonlat$ctr[106023]="BLZ"
lonlat$ctr[110398]="VCT"
lonlat$ctr[126766]="KIR"
lonlat$ctr[146607]="COM"
lonlat$ctr[158876]="MUS"

lonlat=lonlat[, c("lon", "lat", "ctr")] 
lonlat=lonlat[order(lonlat$lon, lonlat$lat), ]

cellwidth=6371000.785*pi/180 			# earthradius [m] times Pi
cellarea=(cellwidth*0.5)*(cellwidth*0.5)*cos(lonlat$lat/180*pi) # replace grid[, 2] with lat/long those ranging from -90 to 90
lonlat$cellarea=cellarea

tgrid=lonlat[lonlat$ctr!=0, ]		
tgrid=tgrid[order(tgrid$lon, tgrid$lat), ]
tgrid_ctr=sort(unique(tgrid$ctr))

# Construct paths relative to the repository
output_dir <- file.path(repo_path, "data", "processed", "other")

# Save the lonlat and tgrid objects
save(lonlat, file = file.path(output_dir, "lonlat_ISIMIP_ctrs.RData"))
save(tgrid,  file = file.path(output_dir, "tgrid_ISIMIP_ctrs.RData"))

# Construct the path to the LU (land use) NetCDF file
LU_path <- file.path(
  repo_path,
  "data", "raw", "other",
  "landuse-15crops_2015soc_annual_1901_2021_landonly.nc"
)

# Open the NetCDF file
LU <- ncdf4::nc_open(filename = LU_path)

# extract coordinates
lon.LU=ncvar_get(LU,"lon")
lat.LU=ncvar_get(LU,"lat")
lonlat.LU <- as.matrix(expand.grid(lon.LU, lat.LU))

vars = c("maize_irrigated", "maize_rainfed", "oil_crops_soybean_irrigated", "oil_crops_soybean_rainfed", "rice_irrigated", "rice_rainfed", "temperate_cereals_irrigated", "temperate_cereals_rainfed")

LU_list = c()

start_LU = "1901"
end_LU = "2021"


for (var in vars){
  # get the actual data
  target=var
  xdata=ncvar_get(LU, target)								
  fillvalue=ncatt_get(LU, target, "_FillValue") 
  unitvalue=ncatt_get(LU, target, "units")		# in %
  xdata[xdata==fillvalue$value] <- NA				# replace fill values with NAs
  # restructure the 3d array of dimensions (x,y,z) to 2d array of (x*y, z)
  xdata_long=matrix(xdata, prod(dim(xdata)[1:2]), dim(xdata)[3])		
  xdata_long=as.data.frame(cbind(lonlat.LU, xdata_long))  
  colnames(xdata_long)=c("lon", "lat", start_LU:end_LU)
  xdata_long=xdata_long[order(xdata_long$lon, xdata_long$lat),] 
  # subset cells to tgrid countries 
  xdata_long=xdata_long[lonlat$ctr%in%tgrid_ctr, ] 
  # Subset cells to only year 2015 for 2015_soc experiment 
  xdata_long = xdata_long[, c("lon", "lat", "2015")]
  # Determine crop area
  xdata_long_crop=xdata_long
  xdata_long_crop[, 3]= xdata_long[, 3]*tgrid$cellarea	 
  
  # Create a name for the data frame based on the variable
  data_frame_name <- paste("LU_", target, sep = "")
  
  # Assign the data frame to the list
  LU_list[[data_frame_name]] <- xdata_long_crop  
  
}

# Save the data in the repository
save(LU_list, file = file.path(output_dir, "LU_list.RData"))

isimip_base <- file.path(
  repo_path,
  "data", "processed", "GGCMI_calendar_adjusted"
)

# List NetCDF files for each crop type
outputs_mai  <- list.files(file.path(isimip_base, "mai"), pattern = "nc")
outputs_ri1  <- list.files(file.path(isimip_base, "ri1"), pattern = "nc")
outputs_ri2  <- list.files(file.path(isimip_base, "ri2"), pattern = "nc")
outputs_soy  <- list.files(file.path(isimip_base, "soy"), pattern = "nc")
outputs_swh  <- list.files(file.path(isimip_base, "swh"), pattern = "nc")
outputs_wwh  <- list.files(file.path(isimip_base, "wwh"), pattern = "nc")

crops <- c("mai", "ri1", "ri2", "soy", "swh", "wwh")

info_list <- list()

for (crop in crops) {
  outputs <- get(paste("outputs_", crop, sep = ""))
  
  info_df <- as.data.frame(matrix(unlist(strsplit(outputs, split="_", fixed=TRUE)), nrow=length(outputs), byrow=TRUE))
  info_df[,10] <- gsub(x=info_df[,10], pattern=".nc", replacement="", fixed=TRUE)
  colnames(info_df) <- c("model", "climate", "obsclim", "soc", "default", "variable", "global", "annual", "startyear", "endyear", "calendar")  # Change if filenames of the model runs are modified
  info_df$startyear <- as.numeric(as.character(info_df$startyear))
  info_df$endyear <- as.numeric(info_df$endyear)
  
  info_list[[paste("info_", crop, sep = "")]] <- info_df
}

models_mai=sort(as.character(unique(info_list$info_mai$model)))
models_ri1=sort(as.character(unique(info_list$info_ri1$model)))
models_ri2=sort(as.character(unique(info_list$info_ri2$model)))
models_soy=sort(as.character(unique(info_list$info_soy$model)))
models_swh=sort(as.character(unique(info_list$info_swh$model)))
models_wwh=sort(as.character(unique(info_list$info_wwh$model)))

crops_mai=sort(as.character(unique(str_extract(info_list$info_mai$variable, "(?<=-)\\w+(?=-)"))))

for (crop in crops) {  # Loop through each crop
  
  # Construct the output list name based on the current crop
  current_list_name <- paste("models_", crop, sep = "")
  current_list <- get(current_list_name)
  current_outputs_name <- paste("outputs_", crop, sep = "")
  current_outputs <- get(current_outputs_name)
  current_info <- info_list[[paste("info_", crop, sep = "")]]
  
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
  
  for (model in current_list) {  # Loop through each model available for the crop
    
    x <- which(current_info$model == model)
    
    # Extract yield + calculate production for irrigated
    nc_file_path_irr <- file.path(
      repo_path,
      "data", "processed", "GGCMI_calendar_adjusted",
      crop,
      current_outputs[x][1]   # Get the irrigated file
    )
    
    # Open the NetCDF file
    yld <- ncdf4::nc_open(filename = nc_file_path_irr)
    
    target <- unique(current_info$variable)[1]
    
    # Extract lon and lat
    xlon <- ncvar_get(yld, "lon")
    xlat <- ncvar_get(yld, "lat")
    xlonlat <- as.matrix(expand.grid(xlon, xlat))
    
    # Get yield data
    xdata <- ncvar_get(yld, target)
    fillvalue <- ncatt_get(yld, target, "_FillValue")
    xdata[xdata == fillvalue$value] <- NA
    yld_long_irr <- matrix(xdata, prod(dim(xdata)[1:2]), dim(xdata)[3])
    
    # Restructure yield data
    yld_long_irr <- as.data.frame(cbind(xlonlat, yld_long_irr))
    colnames(yld_long_irr)[1:2] <- c("lon", "lat")  
    yld_long_irr <- yld_long_irr[order(yld_long_irr$lon, yld_long_irr$lat), ]
    yld_long_irr <- yld_long_irr[lonlat$ctr %in% tgrid_ctr, ]  # Subset to tgrid
    
    # Multiply yield by irrigated crop area for 2015      
    irrigated_LU_list <- LU_list[[paste("LU_", LU_crop, "_irrigated", sep = "")]]
    prod_irr <- yld_long_irr[, 3:ncol(yld_long_irr)] * irrigated_LU_list$`2015` 
    
    # Extract yield + calculate production for rainfed
    nc_file_path_noirr <- file.path(
      repo_path,
      "data", "processed", "GGCMI_calendar_adjusted",
      crop,
      current_outputs[x][2]   # Get the rainfed file
    )
    
    # Open the NetCDF file
    yld <- ncdf4::nc_open(filename = nc_file_path_noirr)
    
    target <- unique(current_info$variable)[2]
    
    # Extract lon and lat
    xlon <- ncvar_get(yld, "lon")
    xlat <- ncvar_get(yld, "lat")
    xlonlat <- as.matrix(expand.grid(xlon, xlat))
    
    # Get yield data
    xdata <- ncvar_get(yld, target)
    xdata[xdata == fillvalue$value] <- NA
    yld_long_rain <- matrix(xdata, prod(dim(xdata)[1:2]), dim(xdata)[3])
    
    # Restructure yield data
    yld_long_rain <- as.data.frame(cbind(xlonlat, yld_long_rain))
    colnames(yld_long_rain)[1:2] <- c("lon", "lat")  
    yld_long_rain <- yld_long_rain[order(yld_long_rain$lon, yld_long_rain$lat), ]
    yld_long_rain <- yld_long_rain[lonlat$ctr %in% tgrid_ctr, ]
    
    # Multiply yield by rainfed crop area for 2015
    rainfed_LU_list <-  LU_list[[paste("LU_", LU_crop, "_rainfed", sep = "")]]
    prod_rain <- yld_long_rain[, 3:ncol(yld_long_rain)] * rainfed_LU_list$`2015`
    nc_close(yld)
    
    # Aggregate production by country
    cntr_prod_irr <- matrix(0, length(tgrid_ctr), ncol(prod_irr))
    cntr_prod_rain <- matrix(0, length(tgrid_ctr), ncol(prod_rain))
    area_irr <- matrix(0, length(tgrid_ctr), 1)
    area_rain <- matrix(0, length(tgrid_ctr), 1)
    
    for (l in 1:ncol(prod_irr)) {  # For each year
      for (k in 1:length(tgrid_ctr)) {  # For each country ID
        cntr_prod_irr[k, l] <- sum(prod_irr[tgrid$ctr == tgrid_ctr[k], l], na.rm = TRUE)
        cntr_prod_rain[k, l] <- sum(prod_rain[tgrid$ctr == tgrid_ctr[k], l], na.rm = TRUE)
        area_irr[k, ] <- sum(irrigated_LU_list[tgrid$ctr == tgrid_ctr[k], ], na.rm = TRUE)
        area_rain[k, ] <- sum(rainfed_LU_list[tgrid$ctr == tgrid_ctr[k], ], na.rm = TRUE)
      } # close year loop
    } # close country loop
    
    # Sum up production for irrigated and rainfed
    cntr_prod <- as.data.frame(cntr_prod_irr + cntr_prod_rain)
    area <- as.data.frame(area_rain + area_irr)
    
    # Calculate yield
    yield <- data.frame(matrix(NA, nrow = nrow(cntr_prod), ncol = ncol(cntr_prod)))
    for (i in 1:ncol(cntr_prod)) {
      yield[, i] <- cntr_prod[, i] / area[, 1]
    }
    
    end_year <- 1900 + ncol(yield)
    years <- 1901:end_year
    colnames(yield) <- c(as.character(years))
    yield$country <- tgrid_ctr
    
    # Construct the full path to the output file
    save_path <- file.path(
      repo_path,
      "data", "processed", "GGCMI_dataframes",
      crop,
      paste0(model, "_", crop, ".RData")
    )
    
    # Save the results
    save(
      cntr_prod_irr, cntr_prod_rain, area_irr, area_rain, area, cntr_prod, yield,
      yld_long_irr, yld_long_rain,
      file = save_path
    )
    
  } # close model loop
} # close crop loop

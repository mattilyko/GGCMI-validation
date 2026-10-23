# Script to adjust the time in the GGCMI yield data from growing seasons to calendar years
# Note that the script handles ALL GGCMI simulation runs/data. For this analysis only yield data from ISIMIP3a, 2015soc, default 
# are used. 
# We downloaded the script from its original Github repository: https://github.com/AgMIP-GGCMI/phase3_processing_pipeline/blob/main/time_adjustment.R

# load libraries
library(ncdf4)
library(testthat)

# settings

base_path <- here::here("data", "raw", "GGCMI_yields")
modelnames <- c(
  "ACEA", "CROVER", "CYGMA1p74", "Daycent", "EPIC-IIASA", "DSSAT-Pythia", 
  "ISAM", "LDNDC", "LPJ-GUESS", "LPJmL", "PEPIC", "pDSSAT", 
  "PROMET", "SIMPLACE-LINTUL5"
)
crops <- c("mai", "ri1", "ri2", "soy", "swh", "mil", "sor", "wwh", "bea") #"sgb","pot","rap","bar")
irrigs <- c("firr", "noirr")
rcsps <- c("picontrol", "historical", "ssp126", "ssp585", "ssp370")
socs <- c("histsoc", "2015soc")
sens <- c("default", "2015co2", "1850co2")
gcms <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
vars <- c("yield", "biom", "cnyield", "plantday", "plantyear", "harvyear", "anthday",
      "matyday", "pirnreq", "aet", "soilmoist1m", "soilmoistmat",
      "transp", "evap", "soilevap", "runoff", "rootm", "tnrup", "tnrin", "tnrloss",
      "n2oemis", "n2emis", "nleach", "tcemis", "ch4emis", "maturitystatus",
      "maturityindex")[1]
calenadar_path <- here::here("data", "raw", "other")
phases <- c("phase3a", "phase3b")

soc_scenario1 <- soc_scenario <- socs[2]
sens_scenario <- sens[1]
bias_adjustment <- "w5e5"
phase <- phases[1]
region <- "landonly"
timestep <- "annual-gs"

if (phase == phases[1]) {
  gcms <- c("gswp3-w5e5")
  rcsps <- "obsclim"
} 


# Function to adjust temporal vec
adjust_temporal_vec <- function(vec, ntimesteps, start_year) {
  # return if there are no valid values in the vec or
  # if first value is start_year (nothing to correct)
  # Check if all values in 'vec' are either NA or 0, or if the first value is equal to 'start_year'
  # The logical operator '|' is used here instead of '||' to perform element-wise comparison
  # '|' returns a vector of the same length as the input vectors, while '||' returns a single logical value
  # This is important because we want to check each element of 'vec' individually, rather than treating it as a single logical value
  
  if (all(is.na(vec) | vec == 0)) {
    return(list(to=c(1:ntimesteps), from=c(1:ntimesteps)))
  }
  if (!is.na(vec[1]) && vec[1] == start_year) {
    return(list(to=c(1:ntimesteps), from=c(1:ntimesteps)))
  }
  
  # calculate the difference between the first valid value and start_year
  # cases with no valid entries are handled above
  first_valid_id <- which(!is.na(vec))[1]
  last_valid_id <- which(!is.na(vec))[length(which(!is.na(vec)))]

  # calculate vec index for the first valid value in vec c(start_year:(start_year + ntimesteps))
  sindex <- which(start_year:(start_year + ntimesteps -1) == vec[first_valid_id])
  if(length(sindex) == 0){
    return(list(to=c(1:ntimesteps), from=c(1:ntimesteps)))
  }

  eindex <-c(last_valid_id - first_valid_id + sindex)
  # check while eindex is larger than ntimesteps, reduce eindex and last_valid_id by 1
  while (eindex > ntimesteps) {
    eindex <- eindex - 1
    last_valid_id <- last_valid_id - 1
  }
  adjusted_vec <- sindex:eindex	
  # fine to return vectors that are shorter than original ntimesteps
  return(list(to=adjusted_vec, from=c(first_valid_id:last_valid_id)))

}


# Main script
main <- function(modelname, climate_forcing, bias_adjustment, climate_scenario,
  soc_scenario, sens_scenario, variable, crop, irrigation, region,
  timestep, start_year, end_year, use_hyear = TRUE) {

  # Concatenate the file name
  file_name <- tolower(paste0(modelname, "_", climate_forcing, "_",
    if(phase == phases[2]) {paste0(bias_adjustment, "_")} else {""},
    climate_scenario, "_", soc_scenario, "_", sens_scenario, "_", variable, "-", crop,
    "-", irrigation, "_", region, "_", timestep, "_", start_year, "_", end_year, ".nc"))
  
  # check if the file exists
  if(
    !file.exists(paste(base_path, crop, file_name, sep = "/"))
  ){
    print(
      paste(
        "File does not exist:", 
        paste(base_path, crop, file_name, sep = "/")
      )
    )  
    return()
  }
  if(use_hyear){
    # concatenate the file name for variable `harvyear`
    file_name_hy <- tolower(paste0(modelname, "_", climate_forcing, "_",
      if(phase == phases[2]) {paste0(bias_adjustment, "_")} else {""},
      climate_scenario, "_", soc_scenario, "_", sens_scenario, "_harvyear-", crop,
      "-", irrigation, "_", region, "_", timestep, "_", start_year, "_", end_year, ".nc"))
    # open the netcdf file for variable `harvyear` 
    nc_hy <- nc_open(
      paste(base_path, crop, file_name_hy, sep = "/")
    )
    # get the data from the netcdf file for variable `harvyear`
    harvyear <-ncvar_get(nc_hy, varid = paste("harvyear", crop, irrigation, sep = "-"))
    # close the netcdf files  
    nc_close(nc_hy)
  } else {
    print(paste("Using crop calendar for", modelname, "for", crop, "and", irrigation))
    file_name_calendar <- file.path(
      calenadar_path, 
      paste(
        "ggcmi-crop-calendar-phase3",
        soc_scenario, crop, irrigation, paste0(region, ".nc"),
        sep = "_"
      )
    )
    nc_calendar <- nc_open(file_name_calendar)
    planting_day <- ncvar_get(nc_calendar, varid = "planting_day")
    maturity_day <- ncvar_get(nc_calendar, varid = "maturity_day")
    nc_close(nc_calendar)
  }
  # open the netcdf file and read metadata on time and space dimensions
  nc <- nc_open(paste(base_path, crop, file_name, sep = "/"))
  ntimesteps <- nc$dim$time$len
  # extract reference year from time dimension unit string "growing seasons since 1661-01-01, 00:00:00"
  ref_year <- strsplit(nc$dim$time$units, "since ")[[1]][2]
  ref_year <- as.numeric(strsplit(ref_year, "-")[[1]][1])
  # compute calendar year of first time step
  # sy <- nc$dim$time$vals[1] + ref_year -1 # this is handled differently per model group, so we always use the start_year
  sy <- start_year

  
  # create empty 3D array of the same dimensions as the variable in the netcdf file
  var_adjusted <- array(NA, dim = c(nc$dim$lon$len, nc$dim$lat$len, ntimesteps))
  # get data from the netcdf file
  var <- ncvar_get(nc, varid = paste(variable, crop, irrigation, sep = "-"))
  
  if (use_hyear) {
    # apply adjust_temporal_vec function for all lon/lat combinations to extract vector indices
    for (i in 1:nc$dim$lon$len) {
      for (j in 1:nc$dim$lat$len) {
        adjust <- adjust_temporal_vec(harvyear[i, j, ], ntimesteps, sy)
        # assign the values in the correct position in the 3D array, rest is NA anyway
        var_adjusted[i, j, adjust$to] <- var[i, j, adjust$from]
      }
    }
  } else {
    # apply adjust_temporal_vec function for all lon/lat combinations to extract vector indices
    for (i in 1:nc$dim$lon$len) {
      for (j in 1:nc$dim$lat$len) {
        if(!is.finite(planting_day[i, j] * maturity_day[i, j])){
          # keep values
          var_adjusted[i, j, ] <- var[i, j, ]
          next
        }
        if(planting_day[i, j] < maturity_day[i, j]){
          # keep values
          var_adjusted[i, j, ] <- var[i, j, ]
        } else {
          # shift by one year
          var_adjusted[i, j, c(2:(dim(var_adjusted)[3]))] <- var[i, j, c(1:(dim(var_adjusted)[3]-1))]
        }
      }
    }
  }

  # create new file name for the adjusted netcdf file
  file_name_adj <- gsub(".nc", "_calendar-year-adjusted.nc", file_name)
  # create new directory for the adjusted netcdf file if it does not exist
  dir_adj <- here::here("data", "processed", "GGCMI_calendar_adjusted", crop)
   # get units, long_name and missval from the original netcdf file
  var_atts <- ncatt_get(nc, varid = paste(variable, crop, irrigation, sep = "-"))
  # get lon dimension attributes
  lon_atts <- ncatt_get(nc, varid = "lon")
  # and values
  lon_vals <- ncvar_get(nc, varid = "lon")
  # create ncvar object for lon
  dim_lon <- ncdim_def("lon", lon_atts$units, lon_vals)

  # get lat dimension attributes
  lat_atts <- ncatt_get(nc, varid = "lat")
  # and values
  lat_vals <- ncvar_get(nc, varid = "lat")
  # create ncvar object for lat
  dim_lat <- ncdim_def("lat", lat_atts$units, lat_vals)
  # get time dimension attributes
  time_atts <- ncatt_get(nc, varid = "time")
  # and values
  time_vals <- ncvar_get(nc, varid = "time")
  # create ncvar object for time
  dim_time <- ncdim_def("time", gsub("growing seasons", "years", time_atts$units), time_vals)

  # catch missing `long_name` attribute for PROMET
  if (is.null(var_atts$long_name)) {
    var_atts$long_name <- var_atts$`long name`
  }
  ncv <- ncvar_def(name = paste(variable, crop, irrigation, sep = "-"), dim = list(dim_lon, dim_lat, dim_time),
    missval = var_atts$'_FillValue', prec = "double",
    units = var_atts$units, longname = var_atts$long_name,
    compression = 9, verbose = FALSE)
  
  # create a new netcdf file with the adjusted variable
  print(paste("Writing", file_name_adj))
  nc_adj <- nc_create(paste(dir_adj, file_name_adj, sep = "/"), list(ncv), verbose = FALSE)

  ncvar_put(nc_adj, varid =ncv, var_adjusted, start = c(1, 1, 1), count = c(-1, -1, -1))


  # get global attributes from the original netcdf file
  global_attributes <- ncatt_get(nc,0)
  # loop through the list of global attributes and add them to the new netcdf file
  for (i in seq_along(global_attributes)) {
    ncatt_put(nc_adj, 0, names(global_attributes[i]), global_attributes[[i]], prec = "char")
  }
  # close the netcdf file
  nc_close(nc_adj)
  nc_close(nc)

  # conduct the test
  if (use_hyear) {
    test_that("Check NA values for i/j combinations with harvyear = 1851", {
      # replace fill values with NA
      var_adjusted[var_adjusted == var_atts$'_FillValue'] <- NA
      # check if all values in var_adjusted are NA for i/j combinations with harvyear = 1851
      expect_true(all(is.na(var_adjusted[, , 1][harvyear[, , 1] == 1851])))
    })
    test_that("Check if all values in var_adjusted are NOT NA for i/j combinations with harvyear = 1850", {
      # need to omit NAs for correct indexing
      harvyear2 <- harvyear
      harvyear2[is.na(harvyear2)] <- -9
      # Check if var_subset has NA values in the first time step
      if (all(!is.na(var_adjusted[, , 1][harvyear2[, , 1] == 1850]))) {
        expect_true(TRUE)
      } else {
        if (!all(is.na(as.vector(harvyear[is.na(as.vector(var))])))) {
          warning(paste("inconsitent NA values across `harvyear` and", variable, "for some pixels"))
        } else {
          warning("Some data points with harvyear = 1850 have NA values while NA values are consistent across `harvyear` and", variable, "for all pixels")
          expect_true(FALSE)
        }
      }
    })
  }
  
}

# loop over modelnames
# loop over gcms
# loop over rcsps
# loop over vars
# loop over crops
# loop over irrigs
for(climate_scenario in rcsps){
  # print(climate_scenario)
  print(paste("Processing", climate_scenario, "scenario"))
  # get the start year from settings
  if (climate_scenario %in% c("picontrol", "historical")) {
    start_year <- 1850
  } else if (climate_scenario %in% c("ssp126", "ssp585", "ssp370")) {
    start_year <- 2015
  } else if (climate_scenario %in% c("obsclim", "counterclim")) {
    start_year <- 1901
  } else {
    start_year <- NA
  }

  # Get target end year
  if (climate_scenario %in% c("historical")) {
    end_year1 <- end_year <- 2014
  } else if (climate_scenario %in% c("picontrol", "ssp126", "ssp585", "ssp370")) {
    end_year1 <- end_year <- 2100
  } else if (climate_scenario %in% c("obsclim", "counterclim")) {
    end_year1 <- end_year <- 2016
  } else {
    end_year <- NA
  }

  for(climate_forcing in gcms){
    for(modelname in modelnames){
      if(phase == phases[1] && modelname == "PROMET" && climate_scenario == "obsclim"){
        end_year <- 2015
      } else {
        end_year <- end_year1
      }
      if(phase == phases[1] && modelname == "LPJ-GUESS" && climate_scenario == "obsclim"){
        soc_scenario <- "histsoc"
      } else {
        soc_scenario <- soc_scenario1
      }
      for(crop in crops){
        for(variable in vars){
          for(irrigation in irrigs){
            # use harvyear if file exists, otherwise fall back to crop calendar
            use_hyear <- TRUE
            filename <- paste(
              base_path, crop,
              tolower(
                paste0(
                  modelname, "_", climate_forcing, "_",
                  if(phase == phases[2]) {paste0(bias_adjustment, "_")} else {""},
                  climate_scenario, "_", soc_scenario, "_", sens_scenario, 
                  "_harvyear-", crop, "-", irrigation, "_", region, "_", timestep, 
                  "_", start_year, "_", end_year, ".nc"
                )
              ),
              sep = "/"
            )
            if (!file.exists(filename)) {
              print(paste("File does not exist:", filename, "skipping", modelname))  
              # if file doesn't exist, use crop calendar
              use_hyear <- FALSE
            }
            print(paste("... Processing", modelname, "for", climate_forcing, "for crop", crop, "and", irrigation, "\n"))
            main(modelname, climate_forcing, bias_adjustment, climate_scenario, soc_scenario, sens_scenario, variable,
              crop, irrigation, region, timestep, start_year, end_year, use_hyear)
          }
        }
      }
    }
  }
}



# Example usage
#main(modelname, climate_forcing, bias_adjustment, climate_scenario, soc_scenario, sens_scenario, variable, crop, irrigation, region, timestep, start_year, end_year)

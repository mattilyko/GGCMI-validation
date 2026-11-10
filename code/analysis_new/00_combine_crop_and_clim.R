library(tidyverse)
library(tidync)
library(furrr)

get_climdata <- function(crop){
  clim_dat_path <- here::here(
    "data", "processed", "extremes_indicators", 
    ifelse(crop == "aggr", "crop_aggregated", "crop_specific")
  )
  
  cold <- tidync(file.path(
    clim_dat_path, paste0("cold_indicators_", crop, ".nc")
  )) |> 
    hyper_tibble(force = TRUE) |> 
    select(lat, lon, year, FCD_p05) |> 
    mutate(across(everything(), as.numeric)) |> 
    mutate(
      cold = (FCD_p05 >= quantile(FCD_p05, 0.95, na.rm = TRUE)),
      .by = c(lat, lon)
    )
  
  hot <- tidync(file.path(
    clim_dat_path, paste0("hot_indicators_", crop, ".nc")
  )) |> 
    hyper_tibble(force = TRUE) |> 
    select(lat, lon, year, FHD_p95) |> 
    mutate(across(everything(), as.numeric)) |> 
    mutate(
      hot = (FHD_p95 >= quantile(FHD_p95, 0.95, na.rm = TRUE)),
      .by = c(lat, lon)
    )
  
  temp <- cold |> 
    full_join(hot, by = c("lat", "lon", "year"))
  
  rm(cold, hot)
  
  precip <- tidync(file.path(
    clim_dat_path, paste0("SPI_", crop, ".nc")
  )) |> 
    activate("D2,D1,D0") |>
    hyper_tibble(force = TRUE) |> 
    select(lat, lon, year, SPI_min, SPI_max) |> 
    mutate(across(everything(), as.numeric)) |> 
    mutate(
      dry = (SPI_min <= -2),
      wet = (SPI_max >= 2),
      .by = c(lat, lon)
    )
  
  clim <- temp |> 
    full_join(precip) |> 
    mutate('hot&dry' = (hot & dry))
  
  rm(temp, precip)
  
  return(clim)
}

process_crop <- function(crop){
  message("Processing crop ", crop, "...")
  
  if (crop == "aggr"){
    load(here::here("data", "processed", "integrated_cropdata", "aggr_bench.RData"))
    load(here::here("data", "processed", "integrated_cropdata", "aggr_sim.RData"))
    
    dat_sim <- detr_sim |> 
      ungroup() |> 
      select(lat, lon, year, model, difftrend_sim, divtrend_sim)
    
    dat_sim <- dat_sim |> 
      summarize(
        divtrend_sim  = median(divtrend_sim, na.rm = TRUE),
        difftrend_sim = median(difftrend_sim, na.rm = TRUE),
        .by = c(lat, lon, year)
      ) |> 
      mutate(model = "ensemble") |> 
      bind_rows(dat_sim)
    
    dat <- detr_obs |> 
      ungroup() |> 
      distinct(lon, lat, year, .keep_all = TRUE) |> 
      rename("yield_obs" = aggregated_yield_obs) |> 
      select(-fitted) |> 
      inner_join(dat_sim, by = c("lat", "lon", "year"))
    
    rm(dat_sim)
  } else {
    load(here::here("data", "processed", "integrated_cropdata", "crop_specific_data.RData"))
    
    dat <- get(paste0("dat_", crop)) |> 
      ungroup() |> 
      mutate(
        total_area = rain_area + irr_area,
        .keep = "unused"
      ) |> 
      rename("yield_obs" = yield) |> 
      select(
        lat, lon, year, ctr, total_area, yield_obs, difftrend_obs, divtrend_obs,
        model, difftrend_sim, divtrend_sim
      )
    
    dat <- dat |> 
      summarize(
        divtrend_sim  = median(divtrend_sim, na.rm = TRUE),
        difftrend_sim = median(difftrend_sim, na.rm = TRUE),
        ctr = first(ctr),
        total_area = first(total_area),
        yield_obs = first(yield_obs),
        difftrend_obs = first(difftrend_obs),
        divtrend_obs = first(divtrend_obs),
        .by = c(lat, lon, year)
      ) |> 
      mutate(model = "ensemble") |> 
      bind_rows(dat)
  }
  
  dat <- dat |> 
    left_join(get_climdata(crop), by = c("lat", "lon", "year"))
  
  qs2::qs_save(dat, file = here::here(
    "data", "processed", "figure_ready_data", paste0("full_data_", crop, ".RDS")
  ))
  
  message("Saved data for crop ", crop, "!")
  
  rm(dat)
}

crops <- c("aggr", "mai", "ri1", "ri2", "soy", "swh", "wwh")

plan(multisession, workers = 7)
future_map(crops, process_crop)
plan(sequential)
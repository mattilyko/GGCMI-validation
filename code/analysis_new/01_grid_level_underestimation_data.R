library(tidyverse)
library(furrr)

compute_loss_freqs <- function(dat, q10, group_var = NULL){
  dat |> 
    summarize(
      loss_freq_obs = mean(divtrend_obs <= q10, na.rm = TRUE),
      loss_freq_sim = mean(divtrend_sim <= q10, na.rm = TRUE),
      ens = mean(model == "ensemble", na.rm = TRUE),
      .by = all_of(c("model", group_var))
    ) |> 
    summarize(
      n_underestimation = sum(!ens & (loss_freq_sim < loss_freq_obs)),
      perc_underestimation = 100 * n_underestimation / sum(!ens),
      loss_freq_obs  = 100 * sum(ens * loss_freq_obs)/sum(ens),
      loss_freq_ens  = 100 * sum(ens * loss_freq_sim)/sum(ens),
      .by = any_of(group_var)
    )
}

get_per_cell_data <- function(cell_dat, country = NULL){
  
  multiple_countries <- (length(unique(cell_dat$ctr)) > 1)
  
  if (!is.null(country)){
    cell_dat <- cell_dat |> 
      filter(ctr == country)
  }
  
  q10 <- quantile(cell_dat$divtrend_obs, 0.1, na.rm = TRUE)
  
  cell_dat_long <- cell_dat |> 
    pivot_longer(
      c(hot, cold, dry, wet, 'hot&dry'),
      names_to = "extreme",
      values_to = "occurence"
    )
  
  res_extremes <- cell_dat_long |> 
    filter(occurence == TRUE) |> 
    compute_loss_freqs(q10 = q10, group_var = "extreme")
  
  res_all <- cell_dat |> 
    compute_loss_freqs(q10 = q10) |> 
    mutate(extreme = "All climates")
  
  res <- bind_rows(res_extremes, res_all)
    
  return(list("res" = res, "multiple_countries" = multiple_countries))
}

process_cell <- function(cell_lat, cell_lon, dat, window_size = 4){
  cell_dat <- dat |> 
    filter(
      abs(lat - cell_lat) <= 0.5 * window_size,
      abs(lon - cell_lon) <= 0.5 * window_size
    )
  
  cell_ctr <- dat |> 
    filter(lat == cell_lat, lon == cell_lon) |> 
    pull(ctr) |> 
    first()
  
  cell_area <- dat |> 
    filter(lat == cell_lat, lon == cell_lon) |> 
    pull(total_area) |> 
    first()
    
  res_full <- get_per_cell_data(cell_dat, country = NULL)
  
  if (res_full$multiple_countries){
    res_ctr <- get_per_cell_data(cell_dat, country = cell_ctr)[["res"]]
  } else {
    res_ctr <- res_full[["res"]]
  }
  
  cell_res <- inner_join(
    res_full[["res"]], res_ctr, by = "extreme", suffix = c("", "_ctr")
  ) |> 
    mutate(
      lat = cell_lat, lon = cell_lon, ctr = cell_ctr, total_area = cell_area
    )
  
  return(cell_res)
}

get_underestimation_data <- function(crop, window_size = 4){
  message("Starting crop ", crop, "...")
  
  dat <- qs2::qs_read(here::here(
    "data", "processed", "figure_ready_data", paste0("full_data_", crop, ".RDS")
  ))
  
  grid <- dat |> 
    distinct(lat, lon) |>
    rename("cell_lat" = lat, "cell_lon" = lon)
  
  res <- future_pmap_dfr(
    grid, 
    process_cell, 
    dat = dat,
    window_size = window_size
  )
  
  qs2::qs_save(res, file = here::here(
    "data", "processed", "figure_ready_data", paste0("underestimation_data_", crop, ".RDS")
  ))
  
  message("Processed crop ", crop, "!")
}

plan(multisession, workers = 10)

crops <- c("aggr", "mai", "ri1", "ri2", "soy", "swh", "wwh")
map(crops, get_underestimation_data, window_size = 4)

plan(sequential)
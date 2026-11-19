library(tidyverse)
library(furrr)
theme_set(
  theme_minimal()+ 
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.box = "vertical",
      legend.box.just = "left"
    )
)

######### bin col fn ################################################

bin_col <- function(
    dat, col, n_bins, binning_method = "equal_frequency"
){
  if (binning_method == "equal_width"){
    min <- min(dat[[col]], na.rm = TRUE)
    max <- max(dat[[col]], na.rm = TRUE)
    bin_size <- (max - min)/n_bins
    breaks <- seq(min, max, length.out = n_bins + 1)
  } else if(binning_method == "equal_frequency"){
    probs <- seq(0, 1, length.out = n_bins + 1)
    breaks <- unique(quantile(dat[[col]], probs = probs, na.rm = TRUE))
  } else{
    stop("Invalid binning method!")
  }
  
  dat <- dat |> 
    mutate(
      bin = cut(.data[[col]], breaks = breaks, include.lowest = TRUE)
    )
  
  if (binning_method == "equal_width") {
    midpoints <- head(breaks, -1) + (diff(breaks) / 2)
    names(midpoints) <- levels(dat$bin)
    
    dat <- dat |> 
      mutate(bin_midpoint = unname(midpoints[as.character(bin)]))
  } else if (binning_method == "equal_frequency") {
    dat <- dat |> 
      group_by(bin, .add = TRUE) |> 
      mutate(bin_midpoint = median(.data[[col]], na.rm = TRUE)) |> 
      ungroup(bin)
  }
  
  return(dat)
}

#################### save plots fn ######################################

save_plot <- function(p, name = NULL, ...){
  if(!is.character(name)){
    name <- deparse(substitute(p))
  }
  
  ggsave(
    plot = p,
    filename = here::here("code", "analysis_bonus", "results", paste0(name, ".png")),
    ...
  )
}

#################### get data ###########################################

crop_countries <- tibble(
  crop = c("mai", "ri1", "ri2", "soy", "swh", "wwh"),
  ctr = list(
    "mai" = c("USA", "ARG", "ETH", "NGA", "UKR"),
    "ri1" = c("CHN", "BGD", "IDN", "THA", "NGA"),
    "ri2" = c("CHN", "VNM", "THA", "IDN", "NGA"),
    "soy" = c("ARG", "BRA", "USA", "IND", "NGA"),
    "swh" = c("CAN", "UKR", "CHN", "KAZ", "PER"),
    "wwh" = c("CAN", "PAK", "CHN", "FRA", "UKR")
  )
) |> 
  unnest(ctr)

weighted_avg <- function(x, weigh_by){
  sum(weigh_by * x, na.rm = TRUE) / sum(weigh_by, na.rm = TRUE)
}

get_line_plot_dat <- function(cr, crop_countries){
  dat <- qs2::qs_read(here::here(
    "data", "processed", "figure_ready_data", paste0("full_data_", cr, ".RDS")
  )) |> 
    inner_join(
      crop_countries |> filter(str_detect(cr, fixed(crop))), 
      by = "ctr"
    )
  
  dat_long <- dat |> 
    distinct(lat, lon, year, ctr, total_area, divtrend_obs, SPI_min, SPI_max, FHD_p95, FCD_p05) |> 
    mutate(model = "benchmark") |> 
    rename("divtrend" = divtrend_obs) |> 
    bind_rows(
      dat |> 
        distinct(lat, lon, year, ctr, total_area, model, divtrend_sim, SPI_min, SPI_max, FHD_p95, FCD_p05) |> 
        rename("divtrend" = divtrend_sim)
    ) |> 
    pivot_longer(c(SPI_min, SPI_max, FHD_p95, FCD_p05), names_to = "indicator") |> 
    mutate(indicator = factor(
      indicator,
      levels = c("SPI_min", "SPI_max", "FHD_p95", "FCD_p05"),
      labels = c("SPI_min", "SPI_max", "FHD", "FCD")
    ))
  
  bins <- dat_long |> 
    filter(model == "benchmark") |> 
    group_by(ctr, indicator) |> 
    group_modify(
      ~{.x |> bin_col("value", 10, "equal_frequency")}
    ) |> 
    ungroup() |> 
    select(lat, lon, year, indicator, bin_midpoint)
  
  quantiles <- dat |> 
    distinct(lat, lon, ctr, total_area) |> 
    rename(
      "cell_lat" = lat, "cell_lon" = lon, "cell_ctr" = ctr, "cell_area" = total_area
    ) |> 
    future_pmap_dfr(function(cell_lat, cell_lon, cell_ctr, cell_area){
      q10 <- dat |> 
        filter(
          abs(lat - cell_lat) <= 2,
          abs(lon - cell_lon) <= 2,
          ctr == cell_ctr
        ) |> 
        summarize(q10 = quantile(divtrend_obs, 0.1, na.rm = TRUE)) |> 
        as.numeric()
      
      return(tibble(
        lat = cell_lat, lon = cell_lon, q10 = q10
      )) 
    })
  
  out <- dat_long |> 
    left_join(quantiles, by = c("lat", "lon")) |> 
    left_join(bins, by = c("lat", "lon", "year", "indicator")) |> 
    summarize(
      loss_freq = weighted_avg(divtrend <= q10, weigh_by = total_area),
      .by = c(ctr, model, indicator, bin_midpoint)
    )
  
  return(out)
}

################### line plots #####################################

process_crop <- function(crop, crop_countries){
  plot_dat <- get_line_plot_dat(crop, crop_countries)
  
  models <- setdiff(unique(plot_dat$model), c("benchmark", "ensemble"))
  
  loss_freq_by_country <- plot_dat |> 
    ggplot(aes(
      x = bin_midpoint,
      y = loss_freq,
      color = model
    )) +
    geom_line(
      data = ~ filter(.x, !(model %in% c("benchmark", "ensemble"))),
      linewidth = 0.75, alpha = 0.5
    ) +
    geom_line(
      data = ~ filter(.x, model == "ensemble"),
      linewidth = 1.25
    ) +
    geom_line(
      data = ~ filter(.x, model == "benchmark"),
      linewidth = 1.25
    ) +
    scale_color_manual(
      values = c(
        benchmark = "black",
        ensemble  = "blue",
        setNames(
          scales::hue_pal()(length(models)),
          models
        )
      )
    ) +
    labs(
      x = "Indicator (binned)",
      y = "Loss frequency",
      color = "Model",
      title = paste0(crop, ": Loss frequency line plots by country")
    ) +
    facet_grid(rows = vars(ctr), cols = vars(indicator), scales = "free_x")+
    scale_y_continuous(labels = scales::label_percent(scale = 100))
  
  save_plot(
    loss_freq_by_country,
    name = paste0("line_plots_", crop),
    height = 10,
    width = 8,
    dpi = 300,
    unit = "in"
  )
}

plan(multisession, workers = 10)

crops <- c("mai", "ri1", "ri2", "soy", "swh", "wwh")
future_map(crops, process_crop, crop_countries = crop_countries)

plan(sequential)

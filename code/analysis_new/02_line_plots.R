library(tidyverse)
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
    filename = here::here("results", "main", paste0(name, ".png")),
    ...
  )
}

#################### line plots ###########################################

weighted_avg <- function(x, weigh_by){
  sum(weigh_by * x) / sum(weigh_by)
}

countries <- c("IND", "USA", "BRA", "FRA", "NGA", "CHN")

dat_aggr <- qs2::qs_read(here::here(
  "data", "processed", "figure_ready_data", "full_data_aggr.RDS"
)) |> 
  filter(ctr %in% countries) |> 
  mutate(
    ctr = factor(
      ctr,
      labels = c("Brazil", "China", "France", "India", "Nigeria", "USA"),
      levels = c("BRA", "CHN", "FRA", "IND", "NGA", "USA")
    )
  )

dat_long <- dat_aggr |> 
  distinct(lat, lon, year, ctr, divtrend_obs, SPI_min, SPI_max, FHD_p95, FCD_p05) |> 
  mutate(model = "benchmark") |> 
  rename("divtrend" = divtrend_obs) |> 
  bind_rows(
    dat_aggr |> 
      distinct(lat, lon, year, ctr, model, divtrend_sim, SPI_min, SPI_max, FHD, FCD) |> 
      rename("divtrend" = divtrend_sim)
  ) |> 
  pivot_longer(c(SPI_min, SPI_max, FHD_p95, FCD_p05), names_to = "indicator") |> 
  mutate(indicator = factor(
    indicator,
    levels = c("SPI_min", "SPI_max", "FHD_p95", "FCD_p05"),
    labels = c("SPI_min", "SPI_max", "FHD", "FCD")
  ))

models <- unique(dat_aggr$model) - 1

bins <- dat_long |> 
  filter(model == "benchmark") |> 
  group_by(ctr, indicator) |> 
  group_modify(
    ~{.x |> bin_col("value", 10, "equal_frequency")}
  ) |> 
  ungroup() |> 
  select(lat, lon, year, indicator, bin_midpoint)

quantiles <- dat_aggr |> 
  distinct(lat, lon, ctr, total_area) |> 
  rename(
    "cell_lat" = lat, "cell_lon" = lon, "cell_ctr" = ctr, "cell_area" = total_area
  ) |> 
  pmap_dfr(function(cell_lat, cell_lon, cell_ctr, cell_area){
    q10 <- dat_aggr |> 
      filter(
        abs(lat - cell_lat) <= 0.5 * window_size,
        abs(lon - cell_lon) <= 0.5 * window_size,
        ctr == cell_ctr
      ) |> 
      summarize(q10 = quantile(divtrend_obs, 0.1, na.rm = TRUE)) |> 
      as_numeric()
     
    return(tibble(
      lat = cell_lat, lon = cell_lon, q10 = q10
    )) 
  })

loss_freq_by_country <- dat_long |> 
  left_join(quantiles, by = c("lat", "lon")) |> 
  left_join(bins, by = c("lat", "lon", "year")) |> 
  summarize(
    loss_freq = sum(total_area * (divtrend <= q10))/sum(total_area),
    .by = c(ctr, model, indicator, bin_midpoint)
  ) |> 
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
    title = "Severe loss frequencies by country"
  ) +
  facet_grid(rows = vars(ctr), cols = vars(indicator), scales = "free_x")+
  scale_y_continuous(labels = scales::label_percent(scale = 100))

save_plot(
  loss_freq_by_country,
  height = 10,
  width = 8,
  dpi = 300,
  unit = "in"
)
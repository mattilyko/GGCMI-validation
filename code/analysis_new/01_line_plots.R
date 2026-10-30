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

#################### get data ###########################################

crop_and_climate_aggr <- qs2::qs_read(here::here(
  "data", "processed", "figure_ready_data", "full_data_aggr.RDS"
))

crop_and_climate_long <- crop_and_climate_aggr |> 
  distinct(lat, lon, year, total_area, ctr, divtrend_obs, SPI_min, SPI_max, FHD, FCD) |> 
  mutate(model = "benchmark") |> 
  rename("divtrend" = divtrend_obs) |> 
  bind_rows(
    crop_and_climate_aggr |> 
      distinct(lat, lon, year, total_area, ctr, model, divtrend_sim, SPI_min, SPI_max, FHD, FCD) |> 
      rename("divtrend" = divtrend_sim)
  ) |> 
  pivot_longer(c(SPI_min, SPI_max, FHD, FCD), names_to = "indicator") |> 
  mutate(indicator = factor(
    indicator,
    levels = c("SPI_min", "SPI_max", "FHD", "FCD")
  ))

models <- unique(crop_and_climate_aggr$model)

############ producer categories ##################################

prod_dat <- crop_and_climate_aggr |> 
  distinct(lat, lon, year, yield_obs, total_area) |> 
  summarize(
    tot_prod = first(total_area) * mean(yield_obs, na.rm = TRUE),
    cell_area = pi * 6371^2 * 0.5 * (1/180) * (
      sin(2*pi*(first(lat) + 0.25)/360) - sin(2*pi*(first(lat) - 0.25)/360)
    ),
    norm_prod = tot_prod / cell_area,
    .by = c(lat, lon)
  ) |> 
  mutate(
    prod_cat = case_when(
      norm_prod >= quantile(norm_prod, 0.9, na.rm = TRUE) ~ "large",
      norm_prod <= quantile(norm_prod, 0.1, na.rm = TRUE) ~ "small",
      .default  =  "medium"
    )
  ) |> 
  select(lat, lon, prod_cat)

prod_quantiles <- crop_and_climate_aggr |> 
  distinct(lat, lon, year, divtrend_obs) |> 
  left_join(prod_dat, by = c("lat", "lon")) |> 
  filter(prod_cat != "medium") |> 
  summarize(
    q10 = quantile(divtrend_obs, 0.1, na.rm = TRUE),
    .by = prod_cat
  )

loss_freq_by_prod_cat <- crop_and_climate_long |>
  left_join(prod_dat, by = c("lat", "lon")) |> 
  left_join(prod_quantiles, by = "prod_cat") |> 
  group_by(prod_cat, indicator, model) |> 
  group_modify(
    ~{.x |> bin_col("value", 10, "equal_frequency")}
  ) |> 
  group_by(bin_midpoint, .add = TRUE) |> 
  summarize(
    loss_freq = mean(divtrend <= q10),
    .groups = "drop"
  ) |> 
  mutate(
    prod_cat = factor(
      prod_cat,
      levels = c("large", "small"),
      labels = c("Large producers", "Small producers")
    )
  ) |> 
  drop_na() |> 
  ggplot(aes(
    x = bin_midpoint,
    y = loss_freq,
    color = model,
    linewidth = model == "benchmark"
  )) +
  geom_line(
    data = ~ filter(.x, model != "benchmark"),
    linewidth = 0.75
  ) +
  geom_line(
    data = ~ filter(.x, model == "benchmark"),
    color = "black",
    linewidth = 1.25
  ) +
  scale_color_manual(
    values = c(
      benchmark = "black",
      setNames(
        scales::hue_pal()(length(models)),
        models
      )
    )
  )+
  labs(
    x = "Indicator (binned)",
    y = "Loss frequency",
    color = "Model",
    title = "Severe loss frequencies by producer category"
  ) +
  facet_grid(rows = vars(prod_cat), cols = vars(indicator), scales = "free_x")+
  scale_y_continuous(labels = scales::label_percent(scale = 100))

save_plot(
  loss_freq_by_prod_cat,
  height = 5,
  width = 8,
  dpi = 300,
  unit = "in"
)

############### countries ##############################

countries <- c("IND", "USA", "BRA", "FRA", "NGA")

country_quantiles <- crop_and_climate_aggr |> 
  distinct(lat, lon, year, ctr, total_area, divtrend_obs) |> 
  filter(ctr %in% countries) |> 
  summarize(
    q10 = modi::weighted.quantile(divtrend_obs, total_area/sum(total_area), 0.1),
    .by = ctr
  )

loss_freq_by_country <- crop_and_climate_long |>
  filter(ctr %in% countries) |> 
  left_join(country_quantiles, by = "ctr") |> 
  group_by(ctr, indicator, model) |> 
  group_modify(
    ~{.x |> bin_col("value", 10, "equal_frequency")}
  ) |> 
  group_by(bin_midpoint, .add = TRUE) |> 
  summarize(
    loss_freq = sum(total_area * (divtrend <= q10))/sum(total_area),
    .groups = "drop"
  ) |> 
  mutate(
    ctr = factor(
      ctr,
      labels = c("USA", "Brazil", "India", "France", "Nigeria"),
      levels = c("USA", "BRA", "IND", "FRA", "NGA")
    )
  ) |> 
  ggplot(aes(
    x = bin_midpoint,
    y = loss_freq,
    color = model,
    linewidth = model == "benchmark"
  )) +
  geom_line(
    data = ~ filter(.x, model != "benchmark"),
    linewidth = 0.75
  ) +
  geom_line(
    data = ~ filter(.x, model == "benchmark"),
    color = "black",
    linewidth = 1.25
  ) +
  scale_color_manual(
    values = c(
      benchmark = "black",
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
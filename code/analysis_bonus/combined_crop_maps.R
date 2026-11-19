library(tidyverse)
library(sf)
library(rnaturalearth)
theme_set(
  theme_minimal()+ 
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.box = "vertical",
      legend.box.just = "left"
    )
)

############ plot fns #####################################################

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

coast <- ne_coastline(scale = "medium")

make_grid_cell <- function(lon, lat) {
  st_polygon(x = list(matrix(c(lon-0.25,lat-0.25,
                               lon-0.25,lat+0.25,
                               lon+0.25,lat+0.25,
                               lon+0.25,lat-0.25,
                               lon-0.25,lat-0.25),
                             nrow = 5, byrow = TRUE)), dim = "XY")
}

global_map <- function(
    dat, col, threshold = Inf, label = waiver(), unit = "", 
    midpoint = 0, title = waiver()
){
  p <- dat |>
    mutate(square = st_sfc(
      mapply(make_grid_cell, lon, lat, SIMPLIFY = FALSE),
      crs = st_crs("WGS84")
    )) |>
    mutate(!!sym(col) := pmin(.data[[col]], threshold)) |> 
    st_as_sf(sf_column_name = "square") |> 
    ggplot() +
    geom_sf(aes(fill = .data[[col]]), lwd = 0) +
    geom_sf(data = coast, color = "black", size = 0.2) + 
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    theme(axis.text = element_blank()) +
    labs(
      title = title,
      fill  = label
    ) + 
    facet_grid(rows = vars(extreme), cols = vars(crop))
  
  if (max(abs(dat[[col]] - midpoint), na.rm = TRUE) < threshold){
    p <- p + 
      scico::scale_fill_scico(
        palette = "vanimo",
        midpoint = midpoint,
        breaks = c(
          min(dat[[col]], na.rm = TRUE), 
          midpoint, 
          max(dat[[col]], na.rm = TRUE)
        ),
        labels = function(x) {
          if (any(abs(x) >= 1000, na.rm = TRUE)) {
            label_vals <- format(x, scientific = TRUE, digits = 2)
          } else {
            label_vals <- round(x, digits = 1)
          }
          paste0(label_vals, unit)
        }
      )
    } else{
      if (min(dat[[col]], na.rm = TRUE) < midpoint - threshold){
        lo <- midpoint - threshold
      } else{
        lo <- min(dat[[col]], na.rm = TRUE)
      }
      
      if (max(dat[[col]], na.rm = TRUE) > midpoint + threshold){
        hi <- midpoint + threshold
      } else{
        hi <- max(dat[[col]], na.rm = TRUE)
      }
      
      p <- p + 
        scico::scale_fill_scico(
          palette = "vanimo",
          midpoint = midpoint,
          limits = c(lo, hi),
          oob = scales::squish,
          breaks = c(lo, midpoint, hi),
          labels = function(x) {
            ifelse(
              x - midpoint == -threshold, 
              paste0("<", -round(threshold, digits = 1), unit), 
              ifelse(
                x - midpoint == threshold,
                paste0(">", round(threshold, digits = 1), unit),
                paste0(round(x, digits = 1), unit)
              )
            )
          }
        )
    }
  
  return(p)
}

######## all crops plot ###########################

crops <- c("mai", "ri1", "ri2", "soy", "swh", "wwh")

combined_dat <- map_dfr(crops, function(crop){
  qs2::qs_read(here::here(
    "data", "processed", "figure_ready_data", paste0("underestimation_data_", crop, ".RDS")
  )) |> 
    mutate(crop = crop)
}) |> 
  filter(extreme %in% c("hot", "dry", "wet", "hot&dry")) |> 
  mutate(
    extreme = factor(
      extreme,
      levels = c("dry", "wet", "hot", "hot&dry"),
      labels = c("Dry", "Wet", "Hot", "Hot & dry")
    )
  )

combined_dat |> 
  global_map(
    "perc_underestimation", unit = "%", midpoint = 50, label = "Percentage of models",
    title = "Percentage of models underestimating frequency of losses in climate extremes"
  ) |> 
  save_plot(
    name = "model_perc_map_combined_crops",
    dpi = 500, height = 9, width = 16, unit = "in"
  )

combined_dat |> 
  mutate(
    ens_underestimation = loss_freq_obs - loss_freq_ens
  ) |> 
  global_map(
    "ens_underestimation", unit = "%", threshold = 25, midpoint = 0, label = "Difference of loss frequencies",
    title = "Loss frequency difference in climate extremes between benchmark and ensemble median"
  ) |> 
  save_plot(
    name = "ens_median_map_combined_crops",
    dpi = 500, height = 9, width = 16, unit = "in"
  )

############# reduced crops #######################

crops_red <- c("mai", "ri1", "soy", "wwh")

red_dat <- map_dfr(crops_red, function(crop){
  qs2::qs_read(here::here(
    "data", "processed", "figure_ready_data", paste0("underestimation_data_", crop, ".RDS")
  )) |> 
    mutate(crop = crop)
}) |> 
  filter(extreme %in% c("hot", "dry", "wet", "hot&dry")) |> 
  mutate(
    extreme = factor(
      extreme,
      levels = c("dry", "wet", "hot", "hot&dry"),
      labels = c("Dry", "Wet", "Hot", "Hot & dry")
    )
  )

red_dat |> 
  global_map(
    "perc_underestimation", unit = "%", midpoint = 50, label = "Percentage of models",
    title = "Percentage of models underestimating frequency of losses in climate extremes"
  ) |> 
  save_plot(
    name = "model_perc_map_reduced_crops",
    dpi = 500, height = 9, width = 16, unit = "in"
  )

red_dat |> 
  mutate(
    ens_underestimation = loss_freq_obs - loss_freq_ens
  ) |> 
  global_map(
    "ens_underestimation", unit = "%", threshold = 25, midpoint = 0, label = "Difference of loss frequencies",
    title = "Loss frequency difference in climate extremes between benchmark and ensemble median"
  ) |> 
  save_plot(
    name = "ens_median_map_reduced_crops",
    dpi = 500, height = 9, width = 16, unit = "in"
  )
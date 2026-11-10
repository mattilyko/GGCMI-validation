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
    filename = here::here("results", "main", paste0(name, ".png")),
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
    dat, col, threshold = Inf, label = waiver(), unit = "", facet_by = NULL, ncol = 1,
    bi_color_scale = FALSE, midpoint = 0, title = waiver()
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
    coord_sf(expand = FALSE) +
    scale_x_continuous(breaks = seq(-180, 180, by = 30)) +
    scale_y_continuous(breaks = seq(-90, 90, by = 30))+
    labs(
      title = title,
      fill  = label
    )
  
  if (!is.null(facet_by)) {
    p <- p + facet_wrap(facets = facet_by, ncol = ncol)
  }
  
  if (bi_color_scale){
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
  } else{
    if (threshold < max(dat[[col]], na.rm = TRUE)){
      p <- p+
        scale_fill_viridis_c(
          option = "magma",
          breaks = c(
            min(dat[[col]], na.rm = TRUE), 
            (min(dat[[col]], na.rm = TRUE) + threshold) / 2, 
            threshold
          ),
          labels = function(x) {
            ifelse(
              x == threshold, 
              paste0(">", round(threshold, digits = 1), unit), 
              paste0(round(x, digits = 1), unit)
            )
          }
        )
    } else{
      threshold <- if(is.finite(threshold)) threshold else max(dat[[col]], na.rm = TRUE)
      p <- p+
        scale_fill_viridis_c(
          option = "magma",
          breaks = c(
            min(dat[[col]], na.rm = TRUE), 
            (min(dat[[col]], na.rm = TRUE) + threshold) / 2, 
            max(threshold)
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
    }
  }
  
  return(p)
}

############## create plots ################################################

get_maps_per_crop <- function(c){
  dat <- qs2::qs_read(here::here(
    "data", "processed", "figure_ready_data", paste0("underestimation_data_", c, ".RDS")
  )) |> 
    mutate(
      extreme = factor(
        extreme,
        levels = c("dry", "wet", "hot", "cold", "hot&dry", "All climates"),
        labels = c("Dry", "Wet", "Hot", "Cold", "Hot & dry", "All climates")
      )
    )
  
  dat |> 
    global_map(
      "perc_underestimation", facet_by = "extreme", ncol = 2, unit = "%",
      bi_color_scale = TRUE, midpoint = 50, label = "Percentage of models",
      title = paste0(c, ": Percentage of models underestimating frequency of severe losses in climate extremes")
    ) |> 
    save_plot(
      name = paste0("model_perc_map_", c),
      dpi = 500, height = 10, width = 10, unit = "in"
    )
  
  dat |> 
    mutate(
      ens_underestimation = loss_freq_obs - loss_freq_ens
    ) |> 
    global_map(
      "ens_underestimation", facet_by = "extreme", ncol = 2, unit = "%", threshold = 25,
      bi_color_scale = TRUE, midpoint = 0, label = "Difference of loss frequencies",
      title = paste0(c, ": Severe loss frequency difference in climate extremes between benchmark and ensemble median")
    ) |> 
    save_plot(
      name = paste0("ens_median_map_", c),
      dpi = 500, height = 10, width = 10, unit = "in"
    )
  
  return()
}

crops <- c("aggr", "mai", "ri1", "ri2", "soy", "swh", "wwh")
map(crops, get_maps_per_crop)

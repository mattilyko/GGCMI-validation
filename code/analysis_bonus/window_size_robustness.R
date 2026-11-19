library(tidyverse)
library(furrr)
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

#### get no. of events per cell dependent on window size ##############

process_cell <- function(cell_lat, cell_lon, dat, window_size = 4){
  dat |> 
    filter(
      abs(lat - cell_lat) <= 0.5 * window_size,
      abs(lon - cell_lon) <= 0.5 * window_size
    ) |> 
    summarize(n = sum(occurence), .by = extreme) |> 
    mutate(lat = cell_lat, lon = cell_lon)
}

check_windowsize <- function(window_size, dat){
  grid <- dat |> 
    distinct(lat, lon) |>
    rename("cell_lat" = lat, "cell_lon" = lon)
  
  res <- future_pmap_dfr(
    grid, 
    process_cell, 
    dat = dat,
    window_size = window_size
  ) |> 
    mutate(window_size = window_size)
  
  return(res)
}

dat <- qs2::qs_read(here::here(
  "data", "processed", "figure_ready_data", "full_data_aggr.RDS"
)) |> 
  distinct(lat, lon, year, .keep_all = TRUE) |> 
  pivot_longer(
    c(hot, cold, dry, wet, 'hot&dry'),
    names_to = "extreme",
    values_to = "occurence"
  ) |>  
  select(lat, lon, year, extreme, occurence)

plan(multisession, workers = availableCores() - 1)

occurence_map_data <- future_pmap_dfr(
  tibble(window_size = 1:4),
  check_windowsize,
  dat = dat
) |> 
  mutate(window_size = case_when(
    window_size == 1 ~ "3x3",
    window_size == 2 ~ "5x5",
    window_size == 3 ~ "7x7",
    .default = "9x9"
  ))

plan(sequential)

qs2::qs_save(occurence_map_data, file = here::here(
  "code", "analysis_bonus", "data", "occurence_map_data.RDS"
))

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

################# histogram #####################################

histogram <- occurence_map_data |> 
  ggplot(aes(x=n))+
  geom_histogram(bins=50)+
  facet_grid(rows = vars(window_size), cols = vars(extreme), scales = "free_x")

save_plot(
  histogram, "extreme_occurences_histogram", width = 8, height = 7, unit = "in", dpi = 300
)

################ maps ############################################

maps <- occurence_map_data |>
  mutate(
    n = case_when(
      n < 5    ~ "<5",
      n <= 10  ~ "6 to 10",
      n <= 25  ~ "11 to 25",
      n <= 50  ~ "26 to 50",
      n <= 100 ~ "51 to 100",
      .default = ">100"
    )
  ) |> 
  mutate(n = factor(
    n, levels = c("<5", "6 to 10", "11 to 25", "26 to 50", "51 to 100", ">100")
  )) |> 
  mutate(square = st_sfc(
    mapply(make_grid_cell, lon, lat, SIMPLIFY = FALSE),
    crs = st_crs("WGS84")
  )) |>
  st_as_sf(sf_column_name = "square") |> 
  ggplot() +
  geom_sf(aes(fill = n), lwd = 0) +
  geom_sf(data = coast, color = "black", size = 0.2) + 
  coord_sf(expand = FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, by = 30)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 30))+
  labs(
    title = "Number of extreme events dependent on neighborhood size",
    fill  = "Number of events"
  ) +
  facet_grid(rows = vars(window_size), cols = vars(extreme))

save_plot(
  maps, "extreme_occurences_maps", width = 10, height = 9, unit = "in", dpi = 500
)

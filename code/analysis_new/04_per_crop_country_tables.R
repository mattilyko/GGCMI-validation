library(tidyverse)
library(gt)

crop_countries <- tibble(
  crop = c("mai", "ri", "soy", "wh"),
  ctr = list(
    "mai" = c("USA", "ARG", "ETH", "KEN", "ZWE", "ZMB", "MWI", "ZAF", "NGA", "UKR"),
    "ri"  = c("CHN", "BGD", "IDN", "VNM", "THA", "NGA"),
    "soy" = c("ARG", "BRA", "USA", "UKR", "IND", "NGA"),
    "wh"  = c("CAN", "PAK", "ETH", "KEN", "FRA", "GER", "UKR", "AUS")
  )
) |> 
  unnest(ctr)

weighted_avg <- function(x, weigh_by){
  sum(weigh_by * x) / sum(weigh_by)
}

weighted_type1_median <- function(x, weigh_by){
  ggdist::weighted_quantile(
    x, probs = 0.5, weights = weigh_by, na.rm = TRUE, type = 1
  )
}

get_crop_country_data <- function(c, crop_countries){
  dat <- qs2::qs_read(here::here(
    "data", "processed", "figure_ready_data", paste0("underestimation_data_", c, ".RDS")
  )) |> 
    inner_join(
      crop_countries |> filter(str_detect(c, fixed(crop))), 
      by = "ctr"
    )
  
  res <- dat |> 
    summarize(
      n_underestimation = weighted_type1_median(n_underestimation_ctr, total_area),
      perc_underestimation = weighted_avg(perc_underestimation_ctr, total_area),
      ens_underestimation  = weighted_avg(
        loss_freq_obs_ctr - loss_freq_ens_ctr, total_area
      ),
      .by = c(ctr, extreme)
    ) |> 
    mutate(crop = c)
  
  return(res)
}

crops <- c("mai", "ri1", "ri2", "soy", "swh", "wwh")
underestimation_dat <- map_dfr(
  crops, get_crop_country_data, crop_countries = crop_countries
)

table_dat <- underestimation_dat |>  
  pivot_wider(
    id_cols = c(ctr, crop),
    names_from = extreme,
    values_from = c(perc_underestimation, ens_underestimation, n_underestimation),
    names_glue = "{.value}_{extreme}"
  ) |>
  arrange(crop, ctr)

###### percentages table ######################

perc_cols <- table_dat |> 
  select(-crop, -ctr, -starts_with("n_")) |> 
  colnames()

country_table_perc <- table_dat |> 
  select(!starts_with("n_")) |> 
  gt(rowname_col = "ctr") |> 
  tab_spanner(
    label = "Percentage of model underestimation",
    columns = starts_with("perc_")
  ) |> 
  tab_spanner(
    label = "Ensemble median underestimation",
    columns = starts_with("ens_")
  ) |> 
  cols_label(
    .list = setNames(
      sub("^(?:[^_]*_){2}", "", perc_cols),
      perc_cols
    )
  ) |> 
  fmt_number(
    columns = -ctr,
    decimals = 1,
    pattern = "{x}%"
  ) |> 
  data_color(
    columns = starts_with("perc_"),
    fn = function(x) scales::col_numeric(
      palette = scico::scico(256, palette = "vanimo"),
      domain = c(0, 100)
    )(pmin(pmax(x, 0), 100))
  ) |>
  data_color(
    columns = starts_with("ens_"),
    fn = function(x) scales::col_numeric(
      palette = scico::scico(256, palette = "vanimo"),
      domain = c(-30, 30)
    )(pmin(pmax(x, -30), 30))
  ) |>
  cols_hide(columns = crop) |>
  tab_style(
    style = cell_borders(sides = "left", weight = px(2), color = "gray50"),
    locations = cells_body(columns = starts_with("ens_")[1])
  ) |>
  tab_style(
    style = cell_borders(sides = "left", weight = px(2), color = "gray50"),
    locations = cells_column_labels(columns = starts_with("ens_")[1])
  )

for (c in rev(sort(unique(table_dat$crop)))) {
  country_table_perc <- country_table_perc |>
    tab_row_group(
      label = c,
      rows = table_dat$crop == c
    )
}

gtsave(
  country_table_perc, 
  filename = "per_crop_table_percentages.png", 
  path = here::here("results", "main")
)

###### count table ##########################

count_cols <- table_dat |> 
  select(-crop, -ctr, -starts_with("perc_")) |> 
  colnames()

country_table_count <- table_dat |> 
  select(!starts_with("perc_")) |> 
  gt(rowname_col = "ctr") |> 
  cols_move_to_start(columns = starts_with("n_")) |>
  tab_spanner(
    label = "Number of underestimating models",
    columns = starts_with("n_")
  ) |> 
  tab_spanner(
    label = "Ensemble median underestimation",
    columns = starts_with("ens_")
  ) |> 
  cols_label(
    .list = setNames(
      sub("^(?:[^_]*_){2}", "", count_cols),
      count_cols
    )
  ) |> 
  fmt_number(
    columns = starts_with("ens_"),
    decimals = 1,
    pattern = "{x}%"
  ) |>
  data_color(
    columns = starts_with("ens_"),
    fn = function(x) scales::col_numeric(
      palette = scico::scico(256, palette = "vanimo"),
      domain = c(-30, 30)
    )(pmin(pmax(x, -30), 30))
  ) |>
  cols_hide(columns = crop) |>
  tab_style(
    style = cell_borders(sides = "left", weight = px(2), color = "gray50"),
    locations = cells_body(columns = starts_with("ens_")[1])
  ) |>
  tab_style(
    style = cell_borders(sides = "left", weight = px(2), color = "gray50"),
    locations = cells_column_labels(columns = starts_with("ens_")[1])
  )

for (crop_name in rev(sort(unique(table_dat$crop)))) {
  crop_max <- table_dat |>
    filter(crop == crop_name) |>
    select(starts_with("n_")) |>
    unlist() |>
    max(na.rm = TRUE)
  
  country_table_count <- country_table_count |>
    data_color(
      columns = starts_with("n_"),
      rows = crop == crop_name,
      fn = local({
        upper <- crop_max
        function(x) scales::col_numeric(
          palette = scico::scico(256, palette = "vanimo"),
          domain = c(0, upper)
        )(pmin(pmax(x, 0), upper))
      })
    ) |>
    tab_row_group(
      label = crop_name,
      rows = table_dat$crop == crop_name
    )
}

gtsave(
  country_table_count, 
  filename = "per_crop_table_counts.png", 
  path = here::here("results", "main")
)

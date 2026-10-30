library(tidyverse)
library(gt)

##################### get data #############################

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

underestimation_dat <- pmap_dfr(
  crop_countries$crop,
  function(c){
    label <- ifelse(c == "ri", "ri1", ifelse(c == "wh", "wwh", c))
    
    crop_dat <- qs2::qs_read(here::here(
      "dat", "processed", "figure_ready_data", paste0("full_data_", label, ".RDS")
    )) |> 
      filter(ctr %in% filter(crop_countries, crop == c)$ctr)
    
    quantiles <- crop_dat |> 
      distinct(lat, lon, year, divtrend_obs) |>
      summarize(
        q10 = modi::weighted.quantile(divtrend_obs, total_area/sum(total_area), 0.1),
        .by = ctr
      )
    
    extremes <- crop_dat |> 
      left_join(quantiles, by = "ctr") |> 
      pivot_longer(
        contains("hot", "cold", "dry", "wet"), 
        names_to = "extreme",
        values_to = "occurence"
      )
      filter(occurence == TRUE) |> 
      summarize(
        loss_freq_obs = mean(divtrend_obs <= q10, na.rm = TRUE),
        loss_freq_sim = mean(divtrend_sim <= q10, na.rm = TRUE),
        ens = mean(model == "ensemble", na.rm = TRUE),
        .by = c(model, ctr, extreme)
      ) |> 
      summarize(
        perc_underestimation = 100 * (
          sum(!ens & loss_freq_sim < loss_freq_obs, na.rm = TRUE)/sum(!ens, na.rm = TRUE)
        ),
        loss_freq_obs  = 100 * (
          sum(ens * loss_freq_obs, na.rm = TRUE)/sum(ens, na.rm = TRUE)
        ),
        loss_freq_ens  = 100 * (
          sum(ens * loss_freq_sim, na.rm = TRUE)/sum(ens, na.rm = TRUE)
        ),
        .by = c(extreme, ctr)
      ) |> 
      mutate(crop = c)
      
      all_climates <- crop_dat |> 
        left_join(quantiles, by = "ctr") |> 
        summarize(
          loss_freq_obs = mean(divtrend_obs <= q10, na.rm = TRUE),
          loss_freq_sim = mean(divtrend_sim <= q10, na.rm = TRUE),
          ens = mean(model == "ensemble", na.rm = TRUE),
          .by = c(model, ctr)
        ) |> 
        summarize(
          perc_underestimation = 100 * (
            sum(!ens & loss_freq_sim < loss_freq_obs, na.rm = TRUE)/sum(!ens, na.rm = TRUE)
          ),
          loss_freq_obs  = 100 * (
            sum(ens * loss_freq_obs, na.rm = TRUE)/sum(ens, na.rm = TRUE)
          ),
          loss_freq_ens  = 100 * (
            sum(ens * loss_freq_sim, na.rm = TRUE)/sum(ens, na.rm = TRUE)
          ),
          .by = c(ctr)
        ) |> 
        mutate(crop = c, extreme = "All climates")
      
      output <- bind_rows(extremes, all_climates)
      
      return(output)
  }
)

###################### make table ######################################

table_dat <- underestimation_dat |> 
  mutate(
    ctr = factor(
      ctr,
      levels = c("USA", "BRA", "ARG", "RUS", "CHN", "IND", "BGD", "EGY", "ZAF", "AUS"),
      labels = c("USA", "Brazil", "Argentina", "Russia", "China", "India", "Bangladesh",
                 "Egypt", "South Africa", "Australia")
    ),
    ens_underestimation = loss_freq_obs - loss_freq_ens,
    .keep = "unused"
  ) |> 
  pivot_wider(
    id_cols = ctr,
    names_from = extreme,
    values_from = c(perc_underestimation, ens_underestimation),
    names_glue = "{.value}_{extreme}"
  )

extreme_levels <- setdiff(names(table_dat), "ctr")

country_table <- table_dat |> 
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
      sub("^(?:[^_]*_){2}", "", extreme_levels),
      extreme_levels
    )
  ) |> 
  fmt_number(
    columns = -ctr,
    decimals = 1,
    pattern = "{x}%"
  ) |> 
  data_color(
    columns = starts_with("perc_"),
    fn = scales::col_numeric(
      palette = scico::scico(256, palette = "vanimo"),
      domain = c(0, 100)
    )
  ) |>
  data_color(
    columns = starts_with("ens_"),
    fn = scales::col_numeric(
      palette = scico::scico(256, palette = "vanimo"),
      domain = c(-30, 30)
    )
  )

for (c in unique(table_dat$crop)) {
  country_table <- country_table |>
    tab_row_group(
      label = c,
      rows = table_dat$crop == c
    )
}

gtsave(
  country_table, 
  filename = "per_crop_table.png", 
  path = here::here("results", "main")
)
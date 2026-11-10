## COMPUTATION OF SPI (STANDARDIZED PRECIPITATION INDEX) FOR CROP GROWING SEASONS
#
# Follows the Copernicus EDO/GDO factsheet (2025):
#   https://drought.emergency.copernicus.eu/data/factsheets/factsheet_spi.pdf
#   and Edwards & McKee (1997) for the per-calendar-month SPI-3 methodology.
#
# Key methodological choices:
#   - SPI-3: backwards-looking 3-month accumulation assigned to each month
#     (e.g. May value = March + April + May totals)
#   - One Gamma distribution fitted per calendar month using all ~39 years
#     of that month's 3-month accumulation → 12 Gamma fits per grid cell
#   - Gamma distribution fitted via L-moments (lmoments3)
#   - Zero precipitation defined as 3-month accumulation < ZERO_THRESHOLD (0.01 mm)
#   - If fewer than n_min (default 10) non-zero baseline points are available,
#     or the Gamma fit fails (error/warning), fall back to an empirical SPI
#     (Weibull-type plotting position, Cunnane constants 0.44 / 0.12) instead of NaN
#   - Non-zero values: p0-adjusted CDF (Eq. 2 in Copernicus factsheet):
#     F_p0(x) = p0 + (1 - p0) * F(x)
#   - Zero values: probability (n_zero + 1) / (2 * (n + 1)) (zero correction)
#   - Probit transform (inverse normal CDF) applied to the probabilities
#
# The per-month SPI logic (_compute_spi_for_month) mirrors the R functions
# new_spi() / gamma_spi() / emp_spi() (gamma_spi with zero_correction=TRUE and
# compute_empirical_on_failure=TRUE).
#
# Adaptation for growing seasons:
#   - Growing season months are identified from the union of crop calendars per cell
#   - For each growing season month, SPI-3 is computed (backwards-looking window)
#   - Per year outputs:
#       SPI_min       : most extreme dry monthly SPI in the growing season
#       SPI_max       : most extreme wet monthly SPI in the growing season
#       SPI_min_month : calendar month (1-12) of SPI_min
#       SPI_max_month : calendar month (1-12) of SPI_max
#       SPI_min_phase : season phase of SPI_min (1=early, 2=mid, 3=late)
#       SPI_max_phase : season phase of SPI_max (1=early, 2=mid, 3=late)
#       SPI_all       : full array of monthly SPI values across growing season months
#                       stored as a (year, month_index, lat, lon) DataArray
#                       where month_index runs 0..max_season_months-1
#
# Output: xr.Dataset saved to
#   GGCMI-validation/data/processed/extremes_indicators/crop_aggregated/

import warnings
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
import pyreadr
from scipy.special import ndtri
from scipy.stats import gamma as sp_gamma
from lmoments3 import distr as lmom_distr


# Identify growing season months for a grid cell

def _get_season_months(start_doy: int, end_doy: int) -> list:
    """
    Return the list of calendar months (1-12) that overlap with the growing
    season defined by start_doy and end_doy (day-of-year).

    Handles cross-year seasons (start_doy > end_doy, Southern Hemisphere).

    Parameters
    ----------
    start_doy : int
        First day of growing season (day-of-year, 1-365).
    end_doy : int
        Last day of growing season (day-of-year, 1-365).

    Returns
    -------
    months : list of int
        Ordered list of calendar months overlapping the growing season.
    """
    # Approximate month boundaries (non-leap year)
    month_start_doy = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
    month_end_doy   = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

    months = []
    if start_doy <= end_doy:
        # Within-year season
        for m in range(1, 13):
            ms = month_start_doy[m - 1]
            me = month_end_doy[m - 1]
            if ms <= end_doy and me >= start_doy:
                months.append(m)
    else:
        # Cross-year season (e.g. Oct -> Mar)
        for m in range(1, 13):
            ms = month_start_doy[m - 1]
            me = month_end_doy[m - 1]
            if me >= start_doy or ms <= end_doy:
                months.append(m)
        # Sort so season order is preserved: start month first
        start_month = next(m for m in range(1, 13)
                           if month_end_doy[m-1] >= start_doy)
        months = ([m for m in months if m >= start_month] +
                  [m for m in months if m < start_month])

    return months


def get_crop_spi(crop, ds_spi):

    crop_names = ['mai', 'ri1', 'ri2', 'soy', 'swh', 'wwh']

    rdata = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")

    if crop == "aggr":
        cropdat_list = []

        for c in crop_names:
            dat = rdata[f"dat_{c}"]
            dat = dat[(dat["rain_area"] > 0) | (dat["irr_area"] > 0)].copy()
            dat["crop"] = c
            cropdat_list.append(dat)

        cropdat = pd.concat(cropdat_list, ignore_index=True)
        cropdat = cropdat.drop_duplicates(
            subset=["lon", "lat", "crop"]
        )

    else:
        cropdat = rdata[f"dat_{crop}"]
        cropdat = cropdat[(cropdat["rain_area"] > 0) | (cropdat["irr_area"] > 0)].copy()

        cropdat = cropdat.drop_duplicates(subset=["lon", "lat"])
        cropdat["crop"] = crop

    season_firr_dict = {}
    season_noirr_dict = {}

    crops_to_process = crop_names if crop == "aggr" else [crop]

    for c in crops_to_process:
        season_firr_dict[c] = xr.open_dataset(
            repo_path /
            f"data/raw/other/ggcmi-crop-calendar-phase3_2015soc_{c}_firr_landonly.nc"
        )

        season_noirr_dict[c] = xr.open_dataset(
            repo_path /
            f"data/raw/other/ggcmi-crop-calendar-phase3_2015soc_{c}_noirr_landonly.nc"
        )

    years = ds_spi["year"].values
    latitudes = np.unique(cropdat["lat"])
    longitudes = np.unique(cropdat["lon"])

    n_years = len(years)
    n_lat = len(latitudes)
    n_lon = len(longitudes)

    SPI_min = np.full((n_years, n_lat, n_lon), np.nan)
    SPI_max = np.full((n_years, n_lat, n_lon), np.nan)

    lat_idx = {lat: i for i, lat in enumerate(latitudes)}
    lon_idx = {lon: i for i, lon in enumerate(longitudes)}

    for (lat, lon), cell_dat in cropdat.groupby(["lat", "lon"]):

        i_lat = lat_idx[lat]
        i_lon = lon_idx[lon]

        start_doy = np.inf
        end_doy = -np.inf

        for _, row in cell_dat.iterrows():

            c = row["crop"]

            rain_area = row["rain_area"]
            irr_area = row["irr_area"]

            if rain_area > 0 and irr_area == 0:
                s = season_noirr_dict[c]["planting_day"].sel(
                    lat=lat, lon=lon
                ).item()

                e = season_noirr_dict[c]["maturity_day"].sel(
                    lat=lat, lon=lon
                ).item()

            elif irr_area > 0 and rain_area == 0:
                s = season_firr_dict[c]["planting_day"].sel(
                    lat=lat, lon=lon
                ).item()

                e = season_firr_dict[c]["maturity_day"].sel(
                    lat=lat, lon=lon
                ).item()

            else:
                s = min(
                    season_firr_dict[c]["planting_day"].sel(
                        lat=lat, lon=lon
                    ).item(),
                    season_noirr_dict[c]["planting_day"].sel(
                        lat=lat, lon=lon
                    ).item()
                )

                e = max(
                    season_firr_dict[c]["maturity_day"].sel(
                        lat=lat, lon=lon
                    ).item(),
                    season_noirr_dict[c]["maturity_day"].sel(
                        lat=lat, lon=lon
                    ).item()
                )

            if not np.isfinite(s) or not np.isfinite(e):
                continue

            start_doy = min(start_doy, s)
            end_doy = max(end_doy, e)

        if start_doy == np.inf or end_doy == -np.inf:
            continue

        # Determine growing season months
        season_months = _get_season_months(
            int(start_doy), int(end_doy)
        )

        # Extract SPI for each year
        for yi, year in enumerate(years):

            row_spi = np.array([
                ds_spi["SPI_all"].sel(
                    year=year,
                    month_index=month - 1,
                    lat=lat,
                    lon=lon
                ).item()
                for month in season_months
            ])

            if np.all(np.isnan(row_spi)):
                continue

            SPI_min[yi, i_lat, i_lon] = np.nanmin(row_spi)
            SPI_max[yi, i_lat, i_lon] = np.nanmax(row_spi)

    # Assemble output dataset
    ds = xr.Dataset(
        {
            "SPI_min": (
                ["year", "lat", "lon"],
                SPI_min
            ),
            "SPI_max": (
                ["year", "lat", "lon"],
                SPI_max
            )
        },
        coords={
            "year": years,
            "lat": latitudes,
            "lon": longitudes
        },
        attrs={
            "description": (
                "Minimum and maximum monthly SPI-6 values "
                "during each calendar year and crop growing season."
            ),
            "crop": crop
        }
    )

    # Close calendar datasets
    for c in crops_to_process:
        season_firr_dict[c].close()
        season_noirr_dict[c].close()

    return ds


# Script entry point

repo_path = Path(__file__).resolve().parent.parent.parent

## 3. Growing season calendars
crop_names = ['aggr', 'mai', 'ri1', 'ri2', 'soy', 'swh', 'wwh']

## 4. Compute SPI
...


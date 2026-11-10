## COMPUTATION OF SPI (STANDARDIZED PRECIPITATION INDEX) FOR CROP GROWING SEASONS
#
# Follows the Copernicus EDO/GDO factsheet (2025):
#   https://drought.emergency.copernicus.eu/data/factsheets/factsheet_spi.pdf
#   and Edwards & McKee (1997) for the per-calendar-month SPI-6 methodology.
#
# Key methodological choices:
#   - SPI-6: backwards-looking 6-month accumulation assigned to each month
#     (e.g. June value = January + Febraury + March + April + May + June totals)
#   - One Gamma distribution fitted per calendar month using all ~39 years
#     of that month's 3-month accumulation → 12 Gamma fits per grid cell
#   - Gamma distribution fitted via L-moments (lmoments3)
#   - Zero precipitation defined as 6-month accumulation < ZERO_THRESHOLD (0.01 mm)
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
#   - For each growing season month, SPI-6 is computed (backwards-looking window)
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

# Accumulations below this value [mm] are treated as "zero" precipitation
ZERO_THRESHOLD = 0.01

# SPI values are capped at +/- SPI_CAP 
SPI_CAP = 3.0

# Timescale of SPI (backwards-looking window in months)
X_MONTHS = 6


# L-moments Gamma fitting

def _lmoments_gamma(data: np.ndarray):
    """
    Fit a two-parameter Gamma distribution using L-moments via lmoments3.

    lmoments3.distr.gam.lmom_fit() computes sample L-moments from the data
    and uses them to estimate the Gamma shape (α) and scale (θ) parameters,
    following Hosking (1990). This is the same methodology used by the
    Copernicus EDO/GDO for SPI computation.

    Parameters
    ----------
    data : np.ndarray
        Non-zero, positive precipitation accumulations.

    Returns
    -------
    alpha : float
        Gamma shape parameter (NaN if fitting fails).
    theta : float
        Gamma scale parameter (NaN if fitting fails).
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            paras = lmom_distr.gam.lmom_fit(data)
            alpha, theta = paras['a'], paras['scale']
        if not (np.isfinite(alpha) and np.isfinite(theta)
                and alpha > 0 and theta > 0):
            return np.nan, np.nan
        return alpha, theta
    except Exception:
        return np.nan, np.nan


# Empirical SPI (fallback)

def _empirical_spi(values: np.ndarray) -> np.ndarray:
    """
    Empirical SPI for a 1-D array of (non-NaN) accumulations.

    Non-zero values: plotting position (n * cume_dist(x) - 0.44) / (n + 0.12),
    where cume_dist(x) = fraction of values <= x (ties included).
    Zero values: (n_zero + 1) / (2 * (n + 1)).
    The probabilities are passed through the probit transform.
    """
    n = len(values)
    zero_mask = values < ZERO_THRESHOLD
    n_zero = np.sum(zero_mask)

    cume_dist = np.searchsorted(np.sort(values), values, side="right") / n
    prob = np.where(
        zero_mask,
        (n_zero + 1) / (2 * (n + 1)),
        (n * cume_dist - 0.44) / (n + 0.12),
    )
    return np.clip(ndtri(prob), -SPI_CAP, SPI_CAP)


# Core SPI calculation for a 1-D array of accumulations (one per year)
# for a single calendar month

def _compute_spi_for_month(accumulations: np.ndarray, n_min: int = 10) -> np.ndarray:
    """
    Given a 1-D array of x-month precipitation accumulations (one value per
    year, all for the same calendar month), return the SPI for each year.

    Gamma-based SPI with zero correction, falling
    back to an empirical SPI when there are too few non-zero values or the
    Gamma fit fails.

    Steps:
      1. Identify zero values (< ZERO_THRESHOLD mm).
      2. If < n_min non-zero values: empirical SPI (no NaN output).
      3. Fit Gamma to non-zero values via L-moments; if the fit fails
         (error, warning or invalid parameters): empirical SPI.
      4. Non-zero values: p0 + (1 - p0) * Gamma CDF.
      5. Zero values: (n_zero + 1) / (2 * (n + 1)).
      6. Apply probit transform -> SPI.

    Parameters
    ----------
    accumulations : np.ndarray, shape (n_years,)
        x-month precipitation totals for one calendar month, one per year [mm].
    n_min : int
        Minimum number of non-zero values required for the Gamma fit.

    Returns
    -------
    spi : np.ndarray, shape (n_years,)
        SPI values; NaN only where the input accumulation is NaN.
    """
    n_years = len(accumulations)
    spi = np.full(n_years, np.nan)

    valid_mask = np.isfinite(accumulations)
    valid_vals = accumulations[valid_mask]

    if len(valid_vals) == 0:
        return spi

    zero_mask = valid_vals < ZERO_THRESHOLD
    nonzero_vals = valid_vals[~zero_mask]

    # Too few non-zero values -> empirical SPI
    if len(nonzero_vals) < n_min:
        spi[valid_mask] = _empirical_spi(valid_vals)
        return spi

    n = len(valid_vals)
    n_zero = np.sum(zero_mask)
    p0 = n_zero / n

    alpha, theta = _lmoments_gamma(nonzero_vals)

    # Gamma fit failed -> empirical SPI
    if np.isnan(alpha) or np.isnan(theta):
        spi[valid_mask] = _empirical_spi(valid_vals)
        return spi

    F_p0 = p0 + (1.0 - p0) * sp_gamma.cdf(valid_vals, a=alpha, scale=theta) # center-of-mass correction
    prob = np.where(zero_mask, (n_zero + 1) / (2 * (n + 1)), F_p0)

    spi[valid_mask] = ndtri(prob)
    spi[valid_mask] = np.clip(spi[valid_mask], -SPI_CAP, SPI_CAP)
    return spi


# Compute x-month backwards-looking precipitation totals for a grid cell

def _get_xmonth_totals(climate_loc: xr.Dataset,
                       calendar_month: int,
                       unique_years: np.ndarray) -> np.ndarray:
    """
    For a single grid cell and a given calendar month, compute the
    backwards-looking x-month precipitation total for each year:
        total(year, month) = sum of daily pr in (month-x_months+1, ..., month)

    Month wrapping is handled (e.g. January looks back to Nov-Dec of prior year).

    Parameters
    ----------
    climate_loc : xr.Dataset
        Daily precipitation at one grid cell, variable "pr" [mm/day].
    calendar_month : int
        Target calendar month (1=Jan ... 12=Dec).
    unique_years : np.ndarray
        Array of years to compute totals for.

    Returns
    -------
    totals : np.ndarray, shape (n_years,)
        x-month precipitation totals [mm]. NaN if data is incomplete.
    """
    totals = np.full(len(unique_years), np.nan)

    # The x months in the backwards window
    months = [((calendar_month - 1 - k) % 12) + 1 for k in range(X_MONTHS - 1, -1, -1)]

    for i, year in enumerate(unique_years):
        total = 0.0
        valid = True
        for j, m in enumerate(months):
            # Months numerically greater than calendar_month must belong
            # to the previous year (they wrapped around the year boundary)
            y = year - 1 if m > calendar_month else year
            data_m = climate_loc.sel(
                time=(climate_loc["time.year"] == y) &
                     (climate_loc["time.month"] == m)
            )["pr"]
            if len(data_m) == 0:
                valid = False
                break
            total += float(data_m.sum().values)
        if valid:
            totals[i] = total

    return totals


# Main function: stat_spi

def stat_spi(climate, grid, n_min=10):
    """
    Compute monthly SPI-6 values over crop growing seasons and derive
    summary indicators per growing season year.

    For each grid cell:
      1. Identify growing season months from the union of crop calendars.
      2. For each calendar month in the growing season, compute backwards-
         looking 6-month precipitation totals across all years.
      3. Fit one Gamma distribution per calendar month (via L-moments)
         using all years as baseline.
      4. Compute SPI-6 for each growing season month x year.
      5. Summarize per year:
           - SPI_min / SPI_max : most extreme dry / wet monthly SPI
           - SPI_min_month / SPI_max_month : calendar month of extremes
           - SPI_min_phase / SPI_max_phase : season phase (1=early,2=mid,3=late)
           - SPI_all : full array of monthly SPI values (year x month_index)

    Parameters
    ----------
    climate : xr.Dataset
        Daily precipitation dataset with variable "pr" [mm/day] and
        dimension "time". Convert from kg m-2 s-1 before calling:
            climate["pr"] = climate["pr"] * 86400
    grid : pd.DataFrame
        Grid of unique (lat,lon) combinations

    Returns
    -------
    ds : xr.Dataset
        Variable (year, month_index, lat, lon):
            SPI_all  — NaN-padded; month_index runs 0..max_season_months-1
    """

    climate['time'] = pd.to_datetime(climate['time'].values)
    unique_years = np.unique(climate['time.year'].values)
    n_years = len(unique_years)
    latitudes  = np.unique(grid["lat"])
    longitudes = np.unique(grid["lon"])
    n_lat = len(latitudes)
    n_lon = len(longitudes)

    # Maximum possible season months (all 12); SPI_all padded with NaN
    n_months = 12

    # Pre-allocate output array
    SPI_all = np.full((n_years, n_months, n_lat, n_lon), np.nan)

    lat_idx = {lat: i for i, lat in enumerate(latitudes)}
    lon_idx = {lon: i for i, lon in enumerate(longitudes)}

    for _, row in grid.iterrows():
        lat = row['lat']
        lon = row['lon']
        i_lat = lat_idx[lat]
        i_lon = lon_idx[lon]

        # Extract precipitation at this grid cell
        climate_loc = climate.sel(lat=lat, lon=lon)

        # Compute SPI-6 per calendar month across all years 
        # spi_matrix: (n_years, n_months)
        spi_matrix = np.full((n_years, n_months), np.nan)

        for j, cal_month in enumerate(range(1, 13)):
            totals = _get_xmonth_totals(climate_loc, cal_month, unique_years)
            spi_matrix[:, j] = _compute_spi_for_month(totals, n_min=n_min)

        # Store full SPI array
        SPI_all[:, :n_months, i_lat, i_lon] = spi_matrix

    #  Package into xr.Dataset 
    coords = {"year": unique_years,
              "month_index": np.arange(12),
              "lat": latitudes, "lon": longitudes}

    def _da(arr, dims, coords, long_name, units=""):
        return xr.DataArray(arr, dims=dims, coords=coords,
                            attrs={"long_name": long_name, "units": units})

    ds = xr.Dataset({
        "SPI_all": _da(
            SPI_all, ["year", "month_index", "lat", "lon"], coords,
            "Monthly SPI-6 for each growing season month (NaN-padded to 12)"),
    }, attrs={
        "description": (
            "SPI-6 computed per growing season month, backwards-looking window. "
            "Gamma fitted per calendar month (empirical fallback). "
            "Follows Copernicus EDO SPI factsheet and Edwards & McKee (1997)."
        ),
        "min_nonzero_baseline_points": n_min,
        "zero_threshold_mm": ZERO_THRESHOLD,
    })

    return ds


# Script entry point

repo_path = Path(__file__).resolve().parent.parent.parent

## 1. Precipitation data (convert kg m-2 s-1 → mm day-1) 
precip1 = xr.open_dataset(repo_path / "data/raw/climdata/gswp3-w5e5_obsclim_pr_landonly_daily_1981_1990.nc", engine='netcdf4')
precip2 = xr.open_dataset(repo_path / "data/raw/climdata/gswp3-w5e5_obsclim_pr_landonly_daily_1991_2000.nc", engine='netcdf4')
precip3 = xr.open_dataset(repo_path / "data/raw/climdata/gswp3-w5e5_obsclim_pr_landonly_daily_2001_2010.nc", engine='netcdf4')
precip4 = xr.open_dataset(repo_path / "data/raw/climdata/gswp3-w5e5_obsclim_pr_landonly_daily_2011_2019.nc", engine='netcdf4')
precip_days = xr.concat([precip1, precip2, precip3, precip4], dim='time')
precip_days["pr"] = precip_days["pr"] * 86400  # kg m-2 s-1 → mm day-1

## 2. Crop data
cropdat_aggr = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/aggr_bench.RData")
cropdat_aggr = next(iter(cropdat_aggr.values()))
grid = cropdat_aggr.drop_duplicates(subset=['lon', 'lat'])[['lon', 'lat']]

## 3. Compute SPI
ds_spi = stat_spi(precip_days, grid)

## 4. Save
out_dir = repo_path / "data/processed/extremes_indicators/crop_aggregated"
ds_spi.to_netcdf(out_dir / "SPI_all.nc")
print("Saved SPI_all.nc")

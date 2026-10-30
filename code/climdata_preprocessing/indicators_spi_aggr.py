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

    valid_mask = ~np.isnan(accumulations)
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


# Phase assignment: early / mid / late within the growing season

def _assign_phase(month_index: int, n_months: int) -> int:
    """
    Assign a season phase (1=early, 2=mid, 3=late) based on the position
    of a month within the growing season.

    Parameters
    ----------
    month_index : int
        0-based index of the month within the ordered season months list.
    n_months : int
        Total number of months in the growing season.

    Returns
    -------
    phase : int
        1 (early), 2 (mid), or 3 (late).
    """
    third = n_months / 3.0
    if month_index < third:
        return 1
    elif month_index < 2 * third:
        return 2
    else:
        return 3


# Main function: season_stat_spi

def season_stat_spi(climate, season_firr_dict, season_noirr_dict, cropdat, n_min=10):
    """
    Compute monthly SPI-3 values over crop growing seasons and derive
    summary indicators per growing season year.

    For each grid cell:
      1. Identify growing season months from the union of crop calendars.
      2. For each calendar month in the growing season, compute backwards-
         looking 3-month precipitation totals across all years.
      3. Fit one Gamma distribution per calendar month (via L-moments)
         using all years as baseline.
      4. Compute SPI-3 for each growing season month x year.
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
    season_firr_dict : dict
        {crop_name: xr.Dataset} with full-irrigation crop calendars.
        Each dataset must contain "planting_day" and "maturity_day".
    season_noirr_dict : dict
        {crop_name: xr.Dataset} with rainfed crop calendars.
    cropdat : pd.DataFrame
        One row per unique (lat, lon) cell with columns:
            "lat", "lon", "rain_area_{crop}", "irr_area_{crop}"
        for each crop in crop_names.

    Returns
    -------
    ds : xr.Dataset
        Variables (year, lat, lon):
            SPI_min, SPI_max, SPI_min_month, SPI_max_month,
            SPI_min_phase, SPI_max_phase
        Variable (year, month_index, lat, lon):
            SPI_all  — NaN-padded; month_index runs 0..max_season_months-1
    """
    crop_names = ['mai', 'ri1', 'ri2', 'soy', 'swh', 'wwh']

    climate['time'] = pd.to_datetime(climate['time'].values)
    unique_years = np.unique(climate['time.year'].values)
    n_years = len(unique_years)
    latitudes  = np.unique(cropdat["lat"])
    longitudes = np.unique(cropdat["lon"])
    n_lat = len(latitudes)
    n_lon = len(longitudes)

    # Maximum possible season months (all 12); SPI_all padded with NaN
    max_months = 12

    # Pre-allocate output arrays
    SPI_min       = np.full((n_years, n_lat, n_lon), np.nan)
    SPI_max       = np.full((n_years, n_lat, n_lon), np.nan)
    SPI_min_month = np.full((n_years, n_lat, n_lon), np.nan)
    SPI_max_month = np.full((n_years, n_lat, n_lon), np.nan)
    SPI_min_phase = np.full((n_years, n_lat, n_lon), np.nan)
    SPI_max_phase = np.full((n_years, n_lat, n_lon), np.nan)
    SPI_all       = np.full((n_years, max_months, n_lat, n_lon), np.nan)

    lat_idx = {lat: i for i, lat in enumerate(latitudes)}
    lon_idx = {lon: i for i, lon in enumerate(longitudes)}

    for _, row in cropdat.iterrows():
        lat = row['lat']
        lon = row['lon']
        i_lat = lat_idx[lat]
        i_lon = lon_idx[lon]

        # Union of growing season days across all crops
        start_doy = np.inf
        end_doy   = -np.inf

        for crop in crop_names:
            rain_area = row[f"rain_area_{crop}"]
            irr_area  = row[f"irr_area_{crop}"]
            if rain_area == 0 and irr_area == 0:
                continue

            s_firr = season_firr_dict[crop]['planting_day'].sel(lat=lat, lon=lon).values
            e_firr = season_firr_dict[crop]['maturity_day'].sel(lat=lat, lon=lon).values
            s_noirr = season_noirr_dict[crop]['planting_day'].sel(lat=lat, lon=lon).values
            e_noirr = season_noirr_dict[crop]['maturity_day'].sel(lat=lat, lon=lon).values

            if rain_area > 0 and irr_area == 0:
                s, e = s_noirr, e_noirr
            elif irr_area > 0 and rain_area == 0:
                s, e = s_firr, e_firr
            else:
                s = min(s_firr, s_noirr)
                e = max(e_firr, e_noirr)

            start_doy = min(start_doy, s)
            end_doy   = max(end_doy,   e)

        if start_doy == np.inf or end_doy == -np.inf:
            continue

        # Growing season months
        season_months = _get_season_months(int(start_doy), int(end_doy))
        n_months = len(season_months)

        # Extract precipitation at this grid cell
        climate_loc = climate.sel(lat=lat, lon=lon)

        # Compute SPI-3 per calendar month across all years 
        # spi_matrix: (n_years, n_months)
        spi_matrix = np.full((n_years, n_months), np.nan)

        for j, cal_month in enumerate(season_months):
            totals = _get_xmonth_totals(climate_loc, cal_month, unique_years)
            spi_matrix[:, j] = _compute_spi_for_month(totals, n_min=n_min)

        # Store full SPI array
        SPI_all[:, :n_months, i_lat, i_lon] = spi_matrix

        # Summary indicators per year
        for yi in range(n_years):
            row_spi = spi_matrix[yi, :]
            if np.all(np.isnan(row_spi)):
                continue

            # Dry extreme (minimum SPI)
            j_min = np.nanargmin(row_spi)
            SPI_min[yi, i_lat, i_lon]       = row_spi[j_min]
            SPI_min_month[yi, i_lat, i_lon] = season_months[j_min]
            SPI_min_phase[yi, i_lat, i_lon] = _assign_phase(j_min, n_months)

            # Wet extreme (maximum SPI)
            j_max = np.nanargmax(row_spi)
            SPI_max[yi, i_lat, i_lon]       = row_spi[j_max]
            SPI_max_month[yi, i_lat, i_lon] = season_months[j_max]
            SPI_max_phase[yi, i_lat, i_lon] = _assign_phase(j_max, n_months)

    #  Package into xr.Dataset 
    coords_3d = {"year": unique_years, "lat": latitudes, "lon": longitudes}
    coords_4d = {"year": unique_years,
                 "month_index": np.arange(max_months),
                 "lat": latitudes, "lon": longitudes}

    def _da(arr, dims, coords, long_name, units=""):
        return xr.DataArray(arr, dims=dims, coords=coords,
                            attrs={"long_name": long_name, "units": units})

    ds = xr.Dataset({
        "SPI_min": _da(
            SPI_min, ["year", "lat", "lon"], coords_3d,
            "Most extreme dry monthly SPI-3 during growing season"),
        "SPI_max": _da(
            SPI_max, ["year", "lat", "lon"], coords_3d,
            "Most extreme wet monthly SPI-3 during growing season"),
        "SPI_min_month": _da(
            SPI_min_month, ["year", "lat", "lon"], coords_3d,
            "Calendar month (1-12) of most extreme dry SPI-3"),
        "SPI_max_month": _da(
            SPI_max_month, ["year", "lat", "lon"], coords_3d,
            "Calendar month (1-12) of most extreme wet SPI-3"),
        "SPI_min_phase": _da(
            SPI_min_phase, ["year", "lat", "lon"], coords_3d,
            "Season phase of most extreme dry SPI-3 (1=early, 2=mid, 3=late)"),
        "SPI_max_phase": _da(
            SPI_max_phase, ["year", "lat", "lon"], coords_3d,
            "Season phase of most extreme wet SPI-3 (1=early, 2=mid, 3=late)"),
        "SPI_all": _da(
            SPI_all, ["year", "month_index", "lat", "lon"], coords_4d,
            "Monthly SPI-3 for each growing season month (NaN-padded to 12)"),
    }, attrs={
        "description": (
            "SPI-3 computed per growing season month, backwards-looking window. "
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
cropdat_swh = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")["dat_swh"]
cropdat_wwh = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")["dat_wwh"]
cropdat_soy = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")["dat_soy"]
cropdat_ri2 = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")["dat_ri2"]
cropdat_ri1 = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")["dat_ri1"]
cropdat_mai = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")["dat_mai"]

cropdat_swh['crop'] = 'swh'
cropdat_wwh['crop'] = 'wwh'
cropdat_soy['crop'] = 'soy'
cropdat_ri2['crop'] = 'ri2'
cropdat_ri1['crop'] = 'ri1'
cropdat_mai['crop'] = 'mai'

cropdat_all = pd.concat([cropdat_swh, cropdat_wwh, cropdat_soy, cropdat_ri2, cropdat_ri1, cropdat_mai])
reshaped_data = cropdat_all.pivot_table(
    index=['lat', 'lon'],
    columns='crop',
    values=['rain_area', 'irr_area'],
    fill_value=0
)
reshaped_data.columns = ['_'.join(col).strip() for col in reshaped_data.columns.values]
reshaped_data.reset_index(inplace=True)
cropdat_unique = reshaped_data.drop_duplicates(subset=['lon', 'lat'])

## 3. Growing season calendars
crop_names = ['mai', 'ri1', 'ri2', 'soy', 'swh', 'wwh']
season_firr_list = [
    xr.open_dataset(repo_path / f"data/raw/other/ggcmi-crop-calendar-phase3_2015soc_{c}_firr_landonly.nc")
    for c in crop_names
]
season_noirr_list = [
    xr.open_dataset(repo_path / f"data/raw/other/ggcmi-crop-calendar-phase3_2015soc_{c}_noirr_landonly.nc")
    for c in crop_names
]
season_firr_dict  = {crop: season for crop, season in zip(crop_names, season_firr_list)}
season_noirr_dict = {crop: season for crop, season in zip(crop_names, season_noirr_list)}

## 4. Compute SPI
ds_spi = season_stat_spi(precip_days, season_firr_dict, season_noirr_dict, cropdat_unique)

## 5. Save
out_dir = repo_path / "data/processed/extremes_indicators/crop_aggregated"
ds_spi.to_netcdf(out_dir / "SPI_aggr.nc")
print("Saved SPI_aggr.nc")

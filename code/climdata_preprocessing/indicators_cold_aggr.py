## COMPUTATION OF COLD EXTREME CLIMATE INDICATORS

# This script processes daily temperature (tasmin) data from the GSWP3-W5E5 dataset to compute 
# annual indicators of cold climate extremes.

# Indicators are calculated:
# - Per calendar year
# - At each grid cell where crops are grown
# - Only using the subset of days within the year when at least one of the five main crops 
#   (maize, soy, spring wheat, winter wheat, or rice) is typically grown at that location

# The growing season periods are derived from crop- and land-use-specific calendars.

# Note: Indicators are not computed per crop. Instead, for each grid cell and year, 
# we consider the union of growing season days across all crops present. These data 
# are further used for the crop aggregated analysis (see code/notebooks/main.qmd) presented in the main paper.

# Output: a single xarray Dataset with all FCD and LCS variants saved to 
# GGCMI-validation/data/processed/extremes_indicators/crop_aggregated/cold_indicators_aggr.nc
#
# Variables in the dataset follow naming conventions:
#   FCD_p{q}   — Frequency of Cold Days using the q-th percentile as threshold (e.g. FCD_p05)
#   LCS_p{q}   — Longest Cold Spell using the q-th percentile as threshold (e.g. LCS_p05)
#   FCD_{T}C   — Frequency of Cold Days using an absolute temperature threshold in Celsius (e.g. FCD_0C)
#   LCS_{T}C   — Longest Cold Spell using an absolute temperature threshold in Celsius (e.g. LCS_0C)
#
# Absolute thresholds cover 0°C, -5°C and -10°C.
# Percentile thresholds: 1st, 5th and 10th percentile.

import numpy as np
import xarray as xr
import pandas as pd 
from datetime import datetime, timedelta 
import pyreadr 
from pathlib import Path

# User-defined base path to the repository
repo_path = Path("/gpfs/home3/kdm27/GGCMI") #change accordingly!

## -- Threshold definitions --------------------------------------------------------------

# Percentile thresholds (lower tail: cold days)
PERCENTILE_THRESHOLDS = [1, 5, 10]  # will produce variables FCD_p01, FCD_p05, FCD_p10, LCS_p01, LCS_p05, LCS_p10

# Absolute temperature thresholds in Celsius; labels used in variable names
# 0°C (frost), -5°C, -10°C
ABSOLUTE_THRESHOLDS = {
    "0C": 0,     # 0°C  — frost
    "-5C": -5,   # -5°C — frost with more severe impacts
    "-10C": -10  # -10°C — severe frost
}

## 1. Get the temperature data and join them

temp1 = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/climdata/gswp3-w5e5_obsclim_tasmin_global_daily_1981_1990.nc", engine='netcdf4')
temp2 = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/climdata/gswp3-w5e5_obsclim_tasmin_global_daily_1991_2000.nc", engine='netcdf4')
temp3 = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/climdata/gswp3-w5e5_obsclim_tasmin_global_daily_2001_2010.nc", engine='netcdf4')
temp4 = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/climdata/gswp3-w5e5_obsclim_tasmin_global_daily_2011_2019.nc", engine='netcdf4')
temp_days = xr.concat([temp1, temp2, temp3, temp4], dim='time')

## 2. Bring in and process crop data and growing season data 

## 2.1 crop data

# load crop data (earlier processed in R as RData)
cropdat_swh = pyreadr.read_r(repo_path / "GGCMI-validation/data/processed/integrated_cropdata/crop_specific_data.RData")["dat_swh"] 
cropdat_wwh = pyreadr.read_r(repo_path / "GGCMI-validation/data/processed/integrated_cropdata/crop_specific_data.RData")["dat_wwh"]
cropdat_ri1 = pyreadr.read_r(repo_path / "GGCMI-validation/data/processed/integrated_cropdata/crop_specific_data.RData")["dat_soy"]
cropdat_ri2 = pyreadr.read_r(repo_path / "GGCMI-validation/data/processed/integrated_cropdata/crop_specific_data.RData")["dat_ri1"]
cropdat_soy = pyreadr.read_r(repo_path / "GGCMI-validation/data/processed/integrated_cropdata/crop_specific_data.RData")["dat_ri2"]
cropdat_mai = pyreadr.read_r(repo_path / "GGCMI-validation/data/processed/integrated_cropdata/crop_specific_data.RData")["dat_mai"] 

# Add a column to each DataFrame indicating the crop name
cropdat_swh['crop'] = 'swh'
cropdat_wwh['crop'] = 'wwh'
cropdat_soy['crop'] = 'soy'
cropdat_ri2['crop'] = 'ri2'
cropdat_ri1['crop'] = 'ri1'
cropdat_mai['crop'] = 'mai'

# Concatenate all dataframes into one, keeping the crop column
cropdat_all = pd.concat([cropdat_swh, cropdat_wwh, cropdat_soy, cropdat_ri2, cropdat_ri1, cropdat_mai])

# Pivot the data to get one row per lat-lon, and separate columns for each crop's rain and irr areas
reshaped_data = cropdat_all.pivot_table(
    index=['lat', 'lon'],
    columns='crop',
    values=['rain_area', 'irr_area'],
    fill_value=0
)

# Flatten multi-level columns
reshaped_data.columns = ['_'.join(col).strip() for col in reshaped_data.columns.values]

# Reset index to make 'lat' and 'lon' regular columns again
reshaped_data.reset_index(inplace=True)

# Drop duplicates to get unique grid cells
cropdat_unique = reshaped_data.drop_duplicates(subset=['lon', 'lat'])

# Crop names for reference
crop_names = ['mai', 'ri1', 'ri2', 'soy', 'swh', 'wwh']

## 2.2 Growing season data

# Load growing seasons
season_firr_swh  = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/other/ggcmi-crop-calendar-phase3_2015soc_swh_firr_landonly.nc") 
season_noirr_swh = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/other/ggcmi-crop-calendar-phase3_2015soc_swh_noirr_landonly.nc")
season_firr_wwh  = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/other/ggcmi-crop-calendar-phase3_2015soc_wwh_firr_landonly.nc") 
season_noirr_wwh = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/other/ggcmi-crop-calendar-phase3_2015soc_wwh_noirr_landonly.nc") 
season_firr_soy  = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/other/ggcmi-crop-calendar-phase3_2015soc_soy_firr_landonly.nc") 
season_noirr_soy = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/other/ggcmi-crop-calendar-phase3_2015soc_soy_noirr_landonly.nc") 
season_firr_mai  = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/other/ggcmi-crop-calendar-phase3_2015soc_mai_firr_landonly.nc") 
season_noirr_mai = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/other/ggcmi-crop-calendar-phase3_2015soc_mai_noirr_landonly.nc")
season_firr_ri1  = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/other/ggcmi-crop-calendar-phase3_2015soc_ri1_firr_landonly.nc") 
season_noirr_ri1 = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/other/ggcmi-crop-calendar-phase3_2015soc_ri1_noirr_landonly.nc") 
season_firr_ri2  = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/other/ggcmi-crop-calendar-phase3_2015soc_ri2_firr_landonly.nc") 
season_noirr_ri2 = xr.open_dataset(repo_path / "GGCMI-validation/data/raw/other/ggcmi-crop-calendar-phase3_2015soc_ri2_noirr_landonly.nc")  

season_firr_list  = [season_firr_mai,  season_firr_ri1,  season_firr_ri2,  season_firr_soy,  season_firr_swh,  season_firr_wwh]
season_noirr_list = [season_noirr_mai, season_noirr_ri1, season_noirr_ri2, season_noirr_soy, season_noirr_swh, season_noirr_wwh]

# Create dictionaries mapping crop names to season data
season_firr_dict  = {crop: season for crop, season in zip(crop_names, season_firr_list)}
season_noirr_dict = {crop: season for crop, season in zip(crop_names, season_noirr_list)}

## 3. Helper function — longest consecutive sequence of 1s

def extreme_length(array):
    # INPUT:
    # - Array where we want to find maximum length of consecutive cold days
    # OUTPUT:
    # - Maximum length of cold days
    max_length = 0
    current_length = 0
    for value in array:
        if value == 1:
            current_length += 1
            max_length = max(max_length, current_length)
        else:
            current_length = 0
    return max_length


## 4. Main computation function
#
# For every grid cell × growing-season year, and for every threshold definition
# (percentile-based and absolute-temperature-based), we compute:
#
#   FCD — Frequency of Cold Days (fraction of growing-season days below threshold)
#   LCS — Longest Cold Spell (max consecutive cold days)
#
# Percentile thresholds are computed from the full multi-year growing-season
# climatology at each grid cell.

def season_stat(climate, season_firr_dict, season_noirr_dict, cropdat,
                percentile_thresholds=PERCENTILE_THRESHOLDS,
                absolute_thresholds=ABSOLUTE_THRESHOLDS):
    # INPUT: 
    # - climate: xarray dataset with climate variable (tasmin)
    # - season_firr_dict: dictionary mapping crop name: full-irrigation crop calendar dataset
    # - season_noirr_dict: dictionary mapping crop name: rainfed crop calendar dataset
    # - cropdat: dataframe with unique locations and land-use info (rain_area / irr_area per crop)
    # - percentile_thresholds: list of lower-tail percentiles defining "cold" 
    # - absolute_thresholds: dict of {label: temp_C} for absolute temperature thresholds
    # OUTPUT:
    # - xr.Dataset with variables FCD_p{q}, LCS_p{q} (percentile-based) and
    #   FCD_{label}, LCS_{label} (absolute-threshold-based)

    # Coordinates 
    climate['time'] = pd.to_datetime(climate['time'].values)
    unique_years = np.unique(climate['time.year'].values)
    latitudes    = np.unique(cropdat["lat"])
    longitudes   = np.unique(cropdat["lon"])

    dims   = ["year", "lat", "lon"]
    coords = {"year": unique_years, "lat": latitudes, "lon": longitudes}
    shape  = (len(unique_years), len(latitudes), len(longitudes))

    # Initialise one DataArray per (indicator × threshold) combination
    # Keys follow the final variable names used in the output Dataset.
    arrays = {}

    thresh_coords = {"lat": latitudes, "lon": longitudes}
    thresh_shape  = (len(latitudes), len(longitudes))

    for q in percentile_thresholds:
        label = f"p{q:02d}"
        arrays[f"FCD_{label}"]    = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        arrays[f"LCS_{label}"]    = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        arrays[f"thresh_{label}"] = xr.DataArray(np.full(thresh_shape, np.nan), coords=thresh_coords, dims=["lat", "lon"])

    for label in absolute_thresholds:
        arrays[f"FCD_{label}"] = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        arrays[f"LCS_{label}"] = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)

    # Helper: fill one year's values for a given binary cold-day series 
    def _fill_year(binary_series, year, lat, lon, total_days):
        """Compute FCD and LCS from a binary (0/1) cold-day DataArray and store
        results for every threshold key whose binary series is provided."""
        fcd_val = binary_series.astype(int).sum(dim="time").values / total_days
        lcs_val = xr.apply_ufunc(
            extreme_length,
            binary_series,
            input_core_dims=[['time']],
            vectorize=True,
            dask='allowed',
            output_dtypes=[int]
        ).values
        return float(fcd_val), float(lcs_val)

    # Main loop over grid cells
    for index, row in cropdat.iterrows():
        lat = row['lat']
        lon = row['lon']

        # Determine the union growing season across all crops at this location
        start_day_loc = np.inf
        end_day_loc   = -np.inf

        for crop in crop_names:
            rain_area = row[f"rain_area_{crop}"]
            irr_area  = row[f"irr_area_{crop}"]

            if rain_area == 0 and irr_area == 0:
                continue

            season_firr  = season_firr_dict[crop]
            season_noirr = season_noirr_dict[crop]

            start_day_firr  = season_firr['planting_day']
            end_day_firr    = season_firr['maturity_day']
            start_day_noirr = season_noirr['planting_day']
            end_day_noirr   = season_noirr['maturity_day']

            if rain_area > 0 and irr_area == 0:        # rainfed only
                start_day_crop = start_day_noirr.sel(lat=lat, lon=lon).values
                end_day_crop   = end_day_noirr.sel(lat=lat, lon=lon).values

            elif irr_area > 0 and rain_area == 0:      # irrigated only
                start_day_crop = start_day_firr.sel(lat=lat, lon=lon).values
                end_day_crop   = end_day_firr.sel(lat=lat, lon=lon).values

            elif rain_area > 0 and irr_area > 0:       # both → extend season
                rainfed_start  = start_day_firr.sel(lat=lat, lon=lon).values
                rainfed_end    = end_day_firr.sel(lat=lat, lon=lon).values
                irrigated_start = start_day_noirr.sel(lat=lat, lon=lon).values
                irrigated_end   = end_day_noirr.sel(lat=lat, lon=lon).values
                start_day_crop = min(rainfed_start, irrigated_start)
                end_day_crop   = max(rainfed_end,   irrigated_end)

            start_day_loc = min(start_day_loc, start_day_crop)
            end_day_loc   = max(end_day_loc,   end_day_crop)

        if start_day_loc == np.inf or end_day_loc == -np.inf:
            continue

        #  Retrieve climate series for this grid cell 
        climate_loc = climate.sel(lat=lat, lon=lon)
        climate_loc['time'] = pd.to_datetime(climate_loc['time'].values)
        dayofyear = climate_loc["time"].dt.dayofyear

        # Branch: growing season within the same calendar year
        if start_day_loc <= end_day_loc:

            climate_loc_growing = climate_loc.where(
                (dayofyear >= start_day_loc) & (dayofyear <= end_day_loc), drop=True
            )
            total_days = end_day_loc - start_day_loc + 1

            # Percentile thresholds — computed from full climatology
            percentile_values = {
                q: climate_loc_growing["tasmin"].quantile(q / 100, dim="time", skipna=True)
                for q in percentile_thresholds
            }

            # Store the threshold value (in K) for each percentile at this grid cell
            for q in percentile_thresholds:
                arrays[f"thresh_p{q:02d}"].loc[lat, lon] = float(percentile_values[q].values)

            for year in unique_years:
                climate_yr = climate_loc.sel(time=climate_loc["time.year"] == year)
                doy_yr     = climate_yr["time"].dt.dayofyear
                climate_gs = climate_yr.sel(
                    time=(doy_yr >= start_day_loc) & (doy_yr <= end_day_loc)
                )
                tasmin_gs  = climate_gs["tasmin"]

                # Percentile-based FCD / LCS
                for q in percentile_thresholds:
                    label  = f"p{q:02d}"
                    binary = (tasmin_gs <= percentile_values[q]).astype(int)
                    fcd, lcs = _fill_year(binary, year, lat, lon, total_days)
                    arrays[f"FCD_{label}"].loc[year, lat, lon] = fcd
                    arrays[f"LCS_{label}"].loc[year, lat, lon] = lcs

                # Absolute-threshold-based FCD / LCS
                for label, thresh_C in absolute_thresholds.items():
                    binary = (tasmin_gs < thresh_C + 273.15).astype(int)
                    fcd, lcs = _fill_year(binary, year, lat, lon, total_days)
                    arrays[f"FCD_{label}"].loc[year, lat, lon] = fcd
                    arrays[f"LCS_{label}"].loc[year, lat, lon] = lcs

        # Branch: growing season spanning two calendar years
        if end_day_loc < start_day_loc:

            climate_loc_growing = climate_loc.where(
                (dayofyear >= start_day_loc) | (dayofyear <= end_day_loc), drop=True
            )
            total_days = (365 - start_day_loc + 1) + end_day_loc  # +1 to include start day

            # Percentile thresholds — computed from full climatology
            percentile_values = {
                q: climate_loc_growing["tasmin"].quantile(q / 100, dim="time", skipna=True)
                for q in percentile_thresholds
            }

            # Store the threshold value (in K) for each percentile at this grid cell
            for q in percentile_thresholds:
                arrays[f"thresh_p{q:02d}"].loc[lat, lon] = float(percentile_values[q].values)

            for year in unique_years:
                current_year = year + 1

                if year == 2019:  # Last year has an incomplete growing season
                    for q in percentile_thresholds:
                        label = f"p{q:02d}"
                        arrays[f"FCD_{label}"].loc[year, lat, lon] = 0
                        arrays[f"LCS_{label}"].loc[year, lat, lon] = 0
                    for label in absolute_thresholds:
                        arrays[f"FCD_{label}"].loc[year, lat, lon] = 0
                        arrays[f"LCS_{label}"].loc[year, lat, lon] = 0
                    continue

                climate_curr  = climate_loc.sel(time=climate_loc["time.year"] == current_year)
                climate_prev  = climate_loc.sel(time=climate_loc["time.year"] == year)
                doy_curr      = climate_curr["time"].dt.dayofyear
                doy_prev      = climate_prev["time"].dt.dayofyear

                # Split growing season across the two calendar years
                climate_gs_curr = climate_curr.sel(time=doy_curr <= end_day_loc)
                climate_gs_prev = climate_prev.sel(time=doy_prev >= start_day_loc)
                tasmin_curr     = climate_gs_curr["tasmin"]
                tasmin_prev     = climate_gs_prev["tasmin"]

                # Percentile-based FCD / LCS
                for q in percentile_thresholds:
                    label    = f"p{q:02d}"
                    thresh   = percentile_values[q]
                    bin_curr = (tasmin_curr <= thresh).astype(int)
                    bin_prev = (tasmin_prev <= thresh).astype(int)
                    fcd_val  = (bin_curr.sum(dim="time") + bin_prev.sum(dim="time")).values / total_days
                    lcs_val  = xr.apply_ufunc(
                        extreme_length,
                        xr.concat([bin_prev, bin_curr], dim='time'),
                        input_core_dims=[['time']],
                        vectorize=True,
                        dask='allowed',
                        output_dtypes=[int]
                    ).values
                    arrays[f"FCD_{label}"].loc[current_year, lat, lon] = float(fcd_val)
                    arrays[f"LCS_{label}"].loc[current_year, lat, lon] = float(lcs_val)

                # Absolute-threshold-based FCD / LCS
                for label, thresh_C in absolute_thresholds.items():
                    bin_curr = (tasmin_curr < thresh_C + 273.15).astype(int)
                    bin_prev = (tasmin_prev < thresh_C + 273.15).astype(int)
                    fcd_val  = (bin_curr.sum(dim="time") + bin_prev.sum(dim="time")).values / total_days
                    lcs_val  = xr.apply_ufunc(
                        extreme_length,
                        xr.concat([bin_prev, bin_curr], dim='time'),
                        input_core_dims=[['time']],
                        vectorize=True,
                        dask='allowed',
                        output_dtypes=[int]
                    ).values
                    arrays[f"FCD_{label}"].loc[current_year, lat, lon] = float(fcd_val)
                    arrays[f"LCS_{label}"].loc[current_year, lat, lon] = float(lcs_val)

    # Assemble Dataset with metadata
    data_vars = {}

    for q in percentile_thresholds:
        label = f"p{q:02d}"
        data_vars[f"FCD_{label}"] = arrays[f"FCD_{label}"]
        data_vars[f"FCD_{label}"].attrs = {
            "long_name": f"Frequency of Cold Days (p{q:02d} threshold)",
            "description": (
                f"Fraction of growing-season days with tasmin ≤ {q}th percentile "
                f"of the local growing-season climatology"
            ),
            "units": "fraction",
        }
        data_vars[f"LCS_{label}"] = arrays[f"LCS_{label}"]
        data_vars[f"LCS_{label}"].attrs = {
            "long_name": f"Longest Cold Spell (p{q:02d} threshold)",
            "description": (
                f"Maximum consecutive growing-season days with tasmin ≤ {q}th percentile "
                f"of the local growing-season climatology"
            ),
            "units": "days",
        }
        data_vars[f"thresh_{label}"] = arrays[f"thresh_{label}"]
        data_vars[f"thresh_{label}"].attrs = {
            "long_name": f"Cold threshold temperature (p{q:02d})",
            "description": (
                f"{q}th percentile of tasmin across all growing-season days in the full climatology, "
                f"used as the cold extreme threshold at each grid cell"
            ),
            "units": "K",
        }

    for label, thresh_C in absolute_thresholds.items():
        data_vars[f"FCD_{label}"] = arrays[f"FCD_{label}"]
        data_vars[f"FCD_{label}"].attrs = {
            "long_name": f"Frequency of Cold Days ({label} / {thresh_C:.0f}°C threshold)",
            "description": (
                f"Fraction of growing-season days with tasmin < {thresh_C}  ({thresh_C:.0f}°C)"
            ),
            "units": "fraction",
            "threshold_C": thresh_C,
        }
        data_vars[f"LCS_{label}"] = arrays[f"LCS_{label}"]
        data_vars[f"LCS_{label}"].attrs = {
            "long_name": f"Longest Cold Spell ({label} / {thresh_C:.0f}°C threshold)",
            "description": (
                f"Maximum consecutive growing-season days with tasmin < {thresh_C} C ({thresh_C:.0f}°C)"
            ),
            "units": "days",
            "threshold_C": thresh_C,
        }

    ds = xr.Dataset(data_vars)
    ds.attrs = {
        "title": "Cold extreme climate indicators — crop-aggregated growing seasons",
        "source": "GSWP3-W5E5 obsclim tasmin, GGCMI Phase 3 crop calendars",
        "percentile_thresholds": str(percentile_thresholds),
        "absolute_thresholds_C": str(absolute_thresholds),
        "created": pd.Timestamp.now().isoformat(),
    }

    return ds


## 5. Perform main computation
ds_cold = season_stat(temp_days, season_firr_dict, season_noirr_dict, cropdat_unique)

## 6. Save to a single NetCDF file
out_path = repo_path / "GGCMI-validation/data/processed/extremes_indicators/crop_aggregated/cold_indicators_aggr.nc"
ds_cold.to_netcdf(out_path)
print(f"Saved dataset to {out_path}")
print(ds_cold)

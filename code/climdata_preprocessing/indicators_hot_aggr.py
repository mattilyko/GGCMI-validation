## COMPUTATION OF HOT EXTREME CLIMATE INDICATORS

# This script processes daily maximum temperature (tasmax) data from the GSWP3-W5E5 dataset to compute 
# annual indicators of heat climate extremes.

# Indicators are calculated:
# - Per calendar year
# - At each grid cell where crops are grown
# - Only using the subset of days within the year when at least one of the five main crops 
#   (maize, soy, spring wheat, winter wheat, or rice) is typically grown at that location

# The growing season periods are derived from crop- and land-use-specific calendars.

# Note: Indicators are not computed per crop. Instead, for each grid cell and year, 
# we consider the union of growing season days across all crops present. These data 
# are further used for the crop aggregated analysis (see code/notebooks/main.qmd) presented in the main paper.

# Output: a single xarray Dataset with all FHD and LHS variants saved to 
# GGCMI-validation/data/processed/extremes_indicators/crop_aggregated/hot_indicators_aggr.nc
#
# Variables in the dataset follow naming conventions:
#   FHD_p{q}   — Frequency of Hot Days using the q-th percentile as threshold (e.g. FHD_p90)
#   LHS_p{q}   — Longest Hot Spell using the q-th percentile as threshold (e.g. LHS_p90)
#   FHD_{T}C  — Frequency of Hot Days using an absolute temperature threshold in Celsius (e.g. FHD_30C)
#   LHS_{T}C   — Longest Hot Spell using an absolute temperature threshold in Celsius (e.g. LHS_30C)
#
# Absolute thresholds cover 30°C, 35°C, and 40°C:
# Percentile thresholds: 90th, 95th, 99th percentile.

import numpy as np
import xarray as xr
import pandas as pd 
from datetime import datetime, timedelta 
import pyreadr 
from pathlib import Path

repo_path = Path(__file__).resolve().parent.parent.parent

## Threshold definitions 

# Percentile thresholds (upper tail → hot days)
PERCENTILE_THRESHOLDS = [90, 95, 99]  # will produce variables FHD_p90, FHD_p95, FHD_p99, LHS_p90, LHS_p95, LHS_p99

# Absolute temperature thresholds in Celsius; labels used in variable names
ABSOLUTE_THRESHOLDS = {
    "30C": 30,   
    "35C": 35,   
    "40C": 40,   
}

## 1. Get the temperature data and join them

temp1 = xr.open_dataset(repo_path / "data/raw/climdata/gswp3-w5e5_obsclim_tasmax_landonly_daily_1981_1990.nc", engine='netcdf4')
temp2 = xr.open_dataset(repo_path / "data/raw/climdata/gswp3-w5e5_obsclim_tasmax_landonly_daily_1991_2000.nc", engine='netcdf4')
temp3 = xr.open_dataset(repo_path / "data/raw/climdata/gswp3-w5e5_obsclim_tasmax_landonly_daily_2001_2010.nc", engine='netcdf4')
temp4 = xr.open_dataset(repo_path / "data/raw/climdata/gswp3-w5e5_obsclim_tasmax_landonly_daily_2011_2019.nc", engine='netcdf4')
temp_days = xr.concat([temp1, temp2, temp3, temp4], dim='time')

## 2. Bring in and process crop data and growing season data 

## 2.1 crop data

# load crop data (earlier processed in R as RData)
cropdat_swh = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")["dat_swh"] 
cropdat_wwh = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")["dat_wwh"]
cropdat_soy = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")["dat_soy"]
cropdat_ri1 = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")["dat_ri1"]
cropdat_ri2 = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")["dat_ri2"]
cropdat_mai = pyreadr.read_r(repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData")["dat_mai"] 

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
season_firr_swh  = xr.open_dataset(repo_path / "data/raw/other/ggcmi-crop-calendar-phase3_2015soc_swh_firr_landonly.nc") 
season_noirr_swh = xr.open_dataset(repo_path / "data/raw/other/ggcmi-crop-calendar-phase3_2015soc_swh_noirr_landonly.nc")
season_firr_wwh  = xr.open_dataset(repo_path / "data/raw/other/ggcmi-crop-calendar-phase3_2015soc_wwh_firr_landonly.nc") 
season_noirr_wwh = xr.open_dataset(repo_path / "data/raw/other/ggcmi-crop-calendar-phase3_2015soc_wwh_noirr_landonly.nc") 
season_firr_soy  = xr.open_dataset(repo_path / "data/raw/other/ggcmi-crop-calendar-phase3_2015soc_soy_firr_landonly.nc") 
season_noirr_soy = xr.open_dataset(repo_path / "data/raw/other/ggcmi-crop-calendar-phase3_2015soc_soy_noirr_landonly.nc") 
season_firr_mai  = xr.open_dataset(repo_path / "data/raw/other/ggcmi-crop-calendar-phase3_2015soc_mai_firr_landonly.nc") 
season_noirr_mai = xr.open_dataset(repo_path / "data/raw/other/ggcmi-crop-calendar-phase3_2015soc_mai_noirr_landonly.nc")
season_firr_ri1  = xr.open_dataset(repo_path / "data/raw/other/ggcmi-crop-calendar-phase3_2015soc_ri1_firr_landonly.nc") 
season_noirr_ri1 = xr.open_dataset(repo_path / "data/raw/other/ggcmi-crop-calendar-phase3_2015soc_ri1_noirr_landonly.nc") 
season_firr_ri2  = xr.open_dataset(repo_path / "data/raw/other/ggcmi-crop-calendar-phase3_2015soc_ri2_firr_landonly.nc") 
season_noirr_ri2 = xr.open_dataset(repo_path / "data/raw/other/ggcmi-crop-calendar-phase3_2015soc_ri2_noirr_landonly.nc")  

season_firr_list  = [season_firr_mai,  season_firr_ri1,  season_firr_ri2,  season_firr_soy,  season_firr_swh,  season_firr_wwh]
season_noirr_list = [season_noirr_mai, season_noirr_ri1, season_noirr_ri2, season_noirr_soy, season_noirr_swh, season_noirr_wwh]

# Create dictionaries mapping crop names to season data
season_firr_dict  = {crop: season for crop, season in zip(crop_names, season_firr_list)}
season_noirr_dict = {crop: season for crop, season in zip(crop_names, season_noirr_list)}

## 3. Helper function — longest consecutive sequence of 1s

def extreme_length(array):
    # INPUT:
    # - Array where we want to find maximum length of consecutive hot days
    # OUTPUT:
    # - Maximum length of hot days
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
#   FHD — Frequency of Hot Days (fraction of growing-season days above threshold)
#   LHS — Longest Hot Spell (max consecutive hot days)
#
# Percentile thresholds are computed from the full multi-year growing-season
# climatology at each grid cell (upper tail: days exceeding the threshold).

def season_stat(climate, season_firr_dict, season_noirr_dict, cropdat,
                percentile_thresholds=PERCENTILE_THRESHOLDS,
                absolute_thresholds=ABSOLUTE_THRESHOLDS):
    # INPUT: 
    # - climate: xarray dataset with climate variable (tasmax)
    # - season_firr_dict: dictionary mapping crop name → full-irrigation crop calendar dataset
    # - season_noirr_dict: dictionary mapping crop name → rainfed crop calendar dataset
    # - cropdat: dataframe with unique locations and land-use info (rain_area / irr_area per crop)
    # - percentile_thresholds: list of upper-tail percentiles defining "hot" 
    # - absolute_thresholds: dict of {label: temp_C} for absolute temperature thresholds
    # OUTPUT:
    # - xr.Dataset with variables FHD_p{q}, LHS_p{q} (percentile-based) and
    #   FHD_{label}, LHS_{label} (absolute-threshold-based)

    # Coordinates 
    climate['time'] = pd.to_datetime(climate['time'].values)
    unique_years = np.unique(climate['time.year'].values)
    latitudes    = np.unique(cropdat["lat"])
    longitudes   = np.unique(cropdat["lon"])

    dims   = ["year", "lat", "lon"]
    coords = {"year": unique_years, "lat": latitudes, "lon": longitudes}
    shape  = (len(unique_years), len(latitudes), len(longitudes))

    # Initialise one DataArray per (indicator × threshold) combination
    arrays = {}

    thresh_coords = {"lat": latitudes, "lon": longitudes}
    thresh_shape  = (len(latitudes), len(longitudes))

    for q in percentile_thresholds:
        label = f"p{q:02d}"
        arrays[f"FHD_{label}"]    = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        arrays[f"LHS_{label}"]    = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        arrays[f"thresh_{label}"] = xr.DataArray(np.full(thresh_shape, np.nan), coords=thresh_coords, dims=["lat", "lon"])

    for label in absolute_thresholds:
        arrays[f"FHD_{label}"] = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        arrays[f"LHS_{label}"] = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)

    # Helper: compute FHD and LHS from a binary hot-day series
    def _fill_year(binary_series, total_days):
        fhd_val = float(binary_series.astype(int).sum(dim="time").values) / total_days
        lhs_val = float(xr.apply_ufunc(
            extreme_length,
            binary_series,
            input_core_dims=[['time']],
            vectorize=True,
            dask='allowed',
            output_dtypes=[int]
        ).values)
        return fhd_val, lhs_val

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
                rainfed_start   = start_day_firr.sel(lat=lat, lon=lon).values
                rainfed_end     = end_day_firr.sel(lat=lat, lon=lon).values
                irrigated_start = start_day_noirr.sel(lat=lat, lon=lon).values
                irrigated_end   = end_day_noirr.sel(lat=lat, lon=lon).values
                start_day_crop  = min(rainfed_start, irrigated_start)
                end_day_crop    = max(rainfed_end,   irrigated_end)

            start_day_loc = min(start_day_loc, start_day_crop)
            end_day_loc   = max(end_day_loc,   end_day_crop)

        if start_day_loc == np.inf or end_day_loc == -np.inf:
            continue

        # Retrieve climate series for this grid cell
        climate_loc = climate.sel(lat=lat, lon=lon)
        climate_loc['time'] = pd.to_datetime(climate_loc['time'].values)
        dayofyear = climate_loc["time"].dt.dayofyear

        # Branch: growing season within the same calendar year
        if start_day_loc <= end_day_loc:

            climate_loc_growing = climate_loc.where(
                (dayofyear >= start_day_loc) & (dayofyear <= end_day_loc), drop=True
            )
            total_days = end_day_loc - start_day_loc + 1

            # Percentile thresholds — computed from full climatology (upper tail)
            percentile_values = {
                q: climate_loc_growing["tasmax"].quantile(q / 100, dim="time", skipna=True)
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
                tasmax_gs = climate_gs["tasmax"]

                # Percentile-based FHD / LHS
                for q in percentile_thresholds:
                    label  = f"p{q:02d}"
                    binary = (tasmax_gs >= percentile_values[q]).astype(int)
                    fhd, lhs = _fill_year(binary, total_days)
                    arrays[f"FHD_{label}"].loc[year, lat, lon] = fhd
                    arrays[f"LHS_{label}"].loc[year, lat, lon] = lhs

                # Absolute-threshold-based FHD / LHS
                for label, thresh_C in absolute_thresholds.items():
                    binary = (tasmax_gs > thresh_C + 273.15).astype(int)
                    fhd, lhs = _fill_year(binary, total_days)
                    arrays[f"FHD_{label}"].loc[year, lat, lon] = fhd
                    arrays[f"LHS_{label}"].loc[year, lat, lon] = lhs

        # Branch: growing season spanning two calendar years
        if end_day_loc < start_day_loc:

            climate_loc_growing = climate_loc.where(
                (dayofyear >= start_day_loc) | (dayofyear <= end_day_loc), drop=True
            )
            total_days = (365 - start_day_loc + 1) + end_day_loc  # +1 to include start day

            # Percentile thresholds — computed from full climatology (upper tail)
            percentile_values = {
                q: climate_loc_growing["tasmax"].quantile(q / 100, dim="time", skipna=True)
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
                        arrays[f"FHD_{label}"].loc[year, lat, lon] = 0
                        arrays[f"LHS_{label}"].loc[year, lat, lon] = 0
                    for label in absolute_thresholds:
                        arrays[f"FHD_{label}"].loc[year, lat, lon] = 0
                        arrays[f"LHS_{label}"].loc[year, lat, lon] = 0
                    continue

                climate_curr = climate_loc.sel(time=climate_loc["time.year"] == current_year)
                climate_prev = climate_loc.sel(time=climate_loc["time.year"] == year)
                doy_curr     = climate_curr["time"].dt.dayofyear
                doy_prev     = climate_prev["time"].dt.dayofyear

                # Split growing season across the two calendar years
                climate_gs_curr = climate_curr.sel(time=doy_curr <= end_day_loc)
                climate_gs_prev = climate_prev.sel(time=doy_prev >= start_day_loc)
                tasmax_curr     = climate_gs_curr["tasmax"]
                tasmax_prev     = climate_gs_prev["tasmax"]

                # Percentile-based FHD / LHS
                for q in percentile_thresholds:
                    label    = f"p{q:02d}"
                    thresh   = percentile_values[q]
                    bin_curr = (tasmax_curr >= thresh).astype(int)
                    bin_prev = (tasmax_prev >= thresh).astype(int)
                    fhd_val  = (bin_curr.sum(dim="time") + bin_prev.sum(dim="time")).values / total_days
                    lhs_val  = xr.apply_ufunc(
                        extreme_length,
                        xr.concat([bin_prev, bin_curr], dim='time'),
                        input_core_dims=[['time']],
                        vectorize=True,
                        dask='allowed',
                        output_dtypes=[int]
                    ).values
                    arrays[f"FHD_{label}"].loc[current_year, lat, lon] = float(fhd_val)
                    arrays[f"LHS_{label}"].loc[current_year, lat, lon] = float(lhs_val)

                # Absolute-threshold-based FHD / LHS
                for label, thresh_C in absolute_thresholds.items():
                    bin_curr = (tasmax_curr > thresh_C + 273.15).astype(int)
                    bin_prev = (tasmax_prev > thresh_C + 273.15).astype(int)
                    fhd_val  = (bin_curr.sum(dim="time") + bin_prev.sum(dim="time")).values / total_days
                    lhs_val  = xr.apply_ufunc(
                        extreme_length,
                        xr.concat([bin_prev, bin_curr], dim='time'),
                        input_core_dims=[['time']],
                        vectorize=True,
                        dask='allowed',
                        output_dtypes=[int]
                    ).values
                    arrays[f"FHD_{label}"].loc[current_year, lat, lon] = float(fhd_val)
                    arrays[f"LHS_{label}"].loc[current_year, lat, lon] = float(lhs_val)

    # Assemble Dataset with metadata
    data_vars = {}

    for q in percentile_thresholds:
        label = f"p{q:02d}"
        data_vars[f"FHD_{label}"] = arrays[f"FHD_{label}"]
        data_vars[f"FHD_{label}"].attrs = {
            "long_name": f"Frequency of Hot Days (p{q:02d} threshold)",
            "description": (
                f"Fraction of growing-season days with tasmax ≥ {q}th percentile "
                f"of the local growing-season climatology"
            ),
            "units": "fraction",
        }
        data_vars[f"LHS_{label}"] = arrays[f"LHS_{label}"]
        data_vars[f"LHS_{label}"].attrs = {
            "long_name": f"Longest Hot Spell (p{q:02d} threshold)",
            "description": (
                f"Maximum consecutive growing-season days with tasmax ≥ {q}th percentile "
                f"of the local growing-season climatology"
            ),
            "units": "days",
        }
        data_vars[f"thresh_{label}"] = arrays[f"thresh_{label}"]
        data_vars[f"thresh_{label}"].attrs = {
            "long_name": f"Hot threshold temperature (p{q:02d})",
            "description": (
                f"{q}th percentile of tasmax across all growing-season days in the full climatology, "
                f"used as the heat extreme threshold at each grid cell"
            ),
            "units": "K",
        }

    for label, thresh_C in absolute_thresholds.items():
        data_vars[f"FHD_{label}"] = arrays[f"FHD_{label}"]
        data_vars[f"FHD_{label}"].attrs = {
            "long_name": f"Frequency of Hot Days ({label} / {thresh_C:.0f}°C threshold)",
            "description": (
                f"Fraction of growing-season days with tasmax > {thresh_C} K ({thresh_C:.0f}°C)"
            ),
            "units": "fraction",
            "threshold_C": thresh_C,
        }
        data_vars[f"LHS_{label}"] = arrays[f"LHS_{label}"]
        data_vars[f"LHS_{label}"].attrs = {
            "long_name": f"Longest Hot Spell ({label} / {thresh_C:.0f}°C threshold)",
            "description": (
                f"Maximum consecutive growing-season days with tasmax > {thresh_C} K ({thresh_C:.0f}°C)"
            ),
            "units": "days",
            "threshold_C": thresh_C,
        }

    ds = xr.Dataset(data_vars)
    ds.attrs = {
        "title": "Hot extreme climate indicators — crop-aggregated growing seasons",
        "source": "GSWP3-W5E5 obsclim tasmax, GGCMI Phase 3 crop calendars",
        "percentile_thresholds": str(percentile_thresholds),
        "absolute_thresholds_C": str(absolute_thresholds),
        "created": pd.Timestamp.now().isoformat(),
    }

    return ds


## 5. Perform main computation
ds_hot = season_stat(temp_days, season_firr_dict, season_noirr_dict, cropdat_unique)

## 6. Save to a single NetCDF file
out_path = repo_path / "data/processed/extremes_indicators/crop_aggregated/hot_indicators_aggr.nc"
ds_hot.to_netcdf(out_path)
print(f"Saved dataset to {out_path}")
print(ds_hot)

## COMPUTATION OF COLD EXTREME CLIMATE INDICATORS — CROP-SPECIFIC

# This script processes daily minimum temperature (tasmin) data from the GSWP3-W5E5 dataset
# to compute annual indicators of cold climate extremes for a single specified crop.
#
# Unlike the crop-aggregated version, growing season days here correspond strictly to
# the target crop's own calendar (firr or noirr depending on land use), not the union
# across all crops at a location.
#
# All crops (mai, ri1, ri2, soy, swh, wwh) are processed in one run.
#
# Output: one NetCDF per crop saved to
#   GGCMI-validation/data/processed/extremes_indicators/crop_specific/cold_indicators_{crop}.nc
#
# Variables follow the same naming conventions as the aggregated version:
#   FCD_p{q}      — Frequency of Cold Days (percentile threshold)
#   LCS_p{q}      — Longest Cold Spell (percentile threshold)
#   thresh_p{q}   — Climatological percentile threshold at each grid cell (K)
#   FCD_{T}C      — Frequency of Cold Days (absolute threshold)
#   LCS_{T}C      — Longest Cold Spell (absolute threshold)

import numpy as np
import xarray as xr
import pandas as pd
import pyreadr
from pathlib import Path

##  Paths and thresholds 

repo_path = Path(__file__).resolve().parent.parent.parent

crop_names = ['mai', 'ri1', 'ri2', 'soy', 'swh', 'wwh']

PERCENTILE_THRESHOLDS = [1, 5, 10]

ABSOLUTE_THRESHOLDS = {
    "0C":   0,    # 0°C  — frost
    "-5C":  -5,   # -5°C — severe frost
    "-10C": -10,  # -10°C — extreme frost
}

## 1. Load temperature data

temp1 = xr.open_dataset(repo_path / "data/raw/climdata/gswp3-w5e5_obsclim_tasmin_landonly_daily_1981_1990.nc", engine='netcdf4')
temp2 = xr.open_dataset(repo_path / "data/raw/climdata/gswp3-w5e5_obsclim_tasmin_landonly_daily_1991_2000.nc", engine='netcdf4')
temp3 = xr.open_dataset(repo_path / "data/raw/climdata/gswp3-w5e5_obsclim_tasmin_landonly_daily_2001_2010.nc", engine='netcdf4')
temp4 = xr.open_dataset(repo_path / "data/raw/climdata/gswp3-w5e5_obsclim_tasmin_landonly_daily_2011_2019.nc", engine='netcdf4')
temp_days = xr.concat([temp1, temp2, temp3, temp4], dim='time')

## 2. Crop data (read the RData file once; select per crop in the loop)

# Map crop name to the correct RData key
rdata_key_map = {
    "swh": "dat_swh",
    "wwh": "dat_wwh",
    "soy": "dat_soy",   
    "ri1": "dat_ri1",   
    "ri2": "dat_ri2",   
    "mai": "dat_mai",
}

rdata_path = repo_path / "data/processed/integrated_cropdata/crop_specific_data.RData"
rdata = pyreadr.read_r(rdata_path)

cal_base = repo_path / "data/raw/other"

## 4. Helper function — longest consecutive sequence of 1s

def extreme_length(array):
    max_length = 0
    current_length = 0
    for value in array:
        if value == 1:
            current_length += 1
            max_length = max(max_length, current_length)
        else:
            current_length = 0
    return max_length

## 5. Main computation function

def season_stat(climate, season_firr, season_noirr, cropdat,
                percentile_thresholds=PERCENTILE_THRESHOLDS,
                absolute_thresholds=ABSOLUTE_THRESHOLDS):
    # INPUT:
    # - climate: xarray dataset with tasmin
    # - season_firr: full-irrigation crop calendar xarray dataset for the target crop
    # - season_noirr: rainfed crop calendar xarray dataset for the target crop
    # - cropdat: dataframe with grid cells where the target crop is grown,
    #            with columns rain_area, irr_area, lat, lon
    # - percentile_thresholds: list of lower-tail percentiles
    # - absolute_thresholds: dict of {label: temp_C}
    # OUTPUT:
    # - xr.Dataset with FCD, LCS, and thresh variables for all threshold definitions

    # Coordinates
    climate['time'] = pd.to_datetime(climate['time'].values)
    unique_years = np.unique(climate['time.year'].values)
    latitudes    = np.unique(cropdat["lat"])
    longitudes   = np.unique(cropdat["lon"])

    dims   = ["year", "lat", "lon"]
    coords = {"year": unique_years, "lat": latitudes, "lon": longitudes}
    shape  = (len(unique_years), len(latitudes), len(longitudes))

    thresh_coords = {"lat": latitudes, "lon": longitudes}
    thresh_shape  = (len(latitudes), len(longitudes))

    # Initialise arrays 
    arrays = {}

    for q in percentile_thresholds:
        label = f"p{q:02d}"
        arrays[f"FCD_{label}"]    = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        arrays[f"LCS_{label}"]    = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        arrays[f"thresh_{label}"] = xr.DataArray(np.full(thresh_shape, np.nan), coords=thresh_coords, dims=["lat", "lon"])

    for label in absolute_thresholds:
        arrays[f"FCD_{label}"] = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)
        arrays[f"LCS_{label}"] = xr.DataArray(np.full(shape, np.nan), coords=coords, dims=dims)

    # Helper 
    def _fill_year(binary_series, total_days):
        fcd_val = float(binary_series.astype(int).sum(dim="time").values) / total_days
        lcs_val = float(xr.apply_ufunc(
            extreme_length,
            binary_series,
            input_core_dims=[['time']],
            vectorize=True,
            dask='allowed',
            output_dtypes=[int]
        ).values)
        return fcd_val, lcs_val

    # Main loop over grid cells 
    for _, row in cropdat.iterrows():
        lat       = row['lat']
        lon       = row['lon']
        rain_area = row['rain_area']
        irr_area  = row['irr_area']

        # Select the growing season for this crop and land-use type
        if rain_area > 0 and irr_area == 0:        # rainfed only
            start_day = season_noirr['planting_day'].sel(lat=lat, lon=lon).values
            end_day   = season_noirr['maturity_day'].sel(lat=lat, lon=lon).values

        elif irr_area > 0 and rain_area == 0:      # irrigated only
            start_day = season_firr['planting_day'].sel(lat=lat, lon=lon).values
            end_day   = season_firr['maturity_day'].sel(lat=lat, lon=lon).values

        elif rain_area > 0 and irr_area > 0:       # both → extend to cover both calendars
            start_day = min(
                season_firr['planting_day'].sel(lat=lat, lon=lon).values,
                season_noirr['planting_day'].sel(lat=lat, lon=lon).values,
            )
            end_day = max(
                season_firr['maturity_day'].sel(lat=lat, lon=lon).values,
                season_noirr['maturity_day'].sel(lat=lat, lon=lon).values,
            )

        else:
            continue  # no area, skip (shouldn't happen after the filter above)

        if np.isnan(start_day) or np.isnan(end_day):
            continue

        # Climate data for this grid cell
        climate_loc = climate.sel(lat=lat, lon=lon)
        climate_loc['time'] = pd.to_datetime(climate_loc['time'].values)
        dayofyear = climate_loc["time"].dt.dayofyear

        # Branch: growing season within the same calendar year
        if start_day <= end_day:

            climate_loc_growing = climate_loc.where(
                (dayofyear >= start_day) & (dayofyear <= end_day), drop=True
            )
            total_days = end_day - start_day + 1

            percentile_values = {
                q: climate_loc_growing["tasmin"].quantile(q / 100, dim="time", skipna=True)
                for q in percentile_thresholds
            }

            for q in percentile_thresholds:
                arrays[f"thresh_p{q:02d}"].loc[lat, lon] = float(percentile_values[q].values)

            for year in unique_years:
                climate_yr = climate_loc.sel(time=climate_loc["time.year"] == year)
                doy_yr     = climate_yr["time"].dt.dayofyear
                tasmin_gs  = climate_yr.sel(
                    time=(doy_yr >= start_day) & (doy_yr <= end_day)
                )["tasmin"]

                for q in percentile_thresholds:
                    label  = f"p{q:02d}"
                    binary = (tasmin_gs <= percentile_values[q]).astype(int)
                    fcd, lcs = _fill_year(binary, total_days)
                    arrays[f"FCD_{label}"].loc[year, lat, lon] = fcd
                    arrays[f"LCS_{label}"].loc[year, lat, lon] = lcs

                for label, thresh_C in absolute_thresholds.items():
                    binary = (tasmin_gs < thresh_C + 273.15).astype(int)
                    fcd, lcs = _fill_year(binary, total_days)
                    arrays[f"FCD_{label}"].loc[year, lat, lon] = fcd
                    arrays[f"LCS_{label}"].loc[year, lat, lon] = lcs

        # Branch: growing season spanning two calendar years
        if end_day < start_day:

            climate_loc_growing = climate_loc.where(
                (dayofyear >= start_day) | (dayofyear <= end_day), drop=True
            )
            total_days = (365 - start_day + 1) + end_day

            percentile_values = {
                q: climate_loc_growing["tasmin"].quantile(q / 100, dim="time", skipna=True)
                for q in percentile_thresholds
            }

            for q in percentile_thresholds:
                arrays[f"thresh_p{q:02d}"].loc[lat, lon] = float(percentile_values[q].values)

            for year in unique_years:
                current_year = year + 1

                if year == 2019:  # incomplete final growing season
                    for q in percentile_thresholds:
                        label = f"p{q:02d}"
                        arrays[f"FCD_{label}"].loc[year, lat, lon] = 0
                        arrays[f"LCS_{label}"].loc[year, lat, lon] = 0
                    for label in absolute_thresholds:
                        arrays[f"FCD_{label}"].loc[year, lat, lon] = 0
                        arrays[f"LCS_{label}"].loc[year, lat, lon] = 0
                    continue

                climate_curr = climate_loc.sel(time=climate_loc["time.year"] == current_year)
                climate_prev = climate_loc.sel(time=climate_loc["time.year"] == year)
                doy_curr     = climate_curr["time"].dt.dayofyear
                doy_prev     = climate_prev["time"].dt.dayofyear

                tasmin_curr = climate_curr.sel(time=doy_curr <= end_day)["tasmin"]
                tasmin_prev = climate_prev.sel(time=doy_prev >= start_day)["tasmin"]

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
                f"Fraction of {crop} growing-season days with tasmin ≤ {q}th percentile "
                f"of the local growing-season climatology"
            ),
            "units": "fraction",
        }
        data_vars[f"LCS_{label}"] = arrays[f"LCS_{label}"]
        data_vars[f"LCS_{label}"].attrs = {
            "long_name": f"Longest Cold Spell (p{q:02d} threshold)",
            "description": (
                f"Maximum consecutive {crop} growing-season days with tasmin ≤ {q}th percentile "
                f"of the local growing-season climatology"
            ),
            "units": "days",
        }
        data_vars[f"thresh_{label}"] = arrays[f"thresh_{label}"]
        data_vars[f"thresh_{label}"].attrs = {
            "long_name": f"Cold threshold temperature (p{q:02d})",
            "description": (
                f"{q}th percentile of tasmin across all {crop} growing-season days "
                f"in the full climatology, used as the cold extreme threshold at each grid cell"
            ),
            "units": "K",
        }

    for label, thresh_C in absolute_thresholds.items():
        data_vars[f"FCD_{label}"] = arrays[f"FCD_{label}"]
        data_vars[f"FCD_{label}"].attrs = {
            "long_name": f"Frequency of Cold Days ({label} / {thresh_C:.0f}°C threshold)",
            "description": (
                f"Fraction of {crop} growing-season days with tasmin < {thresh_C}°C"
            ),
            "units": "fraction",
            "threshold_C": thresh_C,
        }
        data_vars[f"LCS_{label}"] = arrays[f"LCS_{label}"]
        data_vars[f"LCS_{label}"].attrs = {
            "long_name": f"Longest Cold Spell ({label} / {thresh_C:.0f}°C threshold)",
            "description": (
                f"Maximum consecutive {crop} growing-season days with tasmin < {thresh_C}°C"
            ),
            "units": "days",
            "threshold_C": thresh_C,
        }

    ds = xr.Dataset(data_vars)
    ds.attrs = {
        "title": f"Cold extreme climate indicators — {crop} growing season",
        "crop": crop,
        "source": "GSWP3-W5E5 obsclim tasmin, GGCMI Phase 3 crop calendars",
        "percentile_thresholds": str(percentile_thresholds),
        "absolute_thresholds_C": str(absolute_thresholds),
        "created": pd.Timestamp.now().isoformat(),
    }

    return ds


## 6. Loop over all crops, run and save

out_dir = repo_path / "data/processed/extremes_indicators/crop_specific"
out_dir.mkdir(parents=True, exist_ok=True)

# which is the loop variable below.
for crop in crop_names:
    print(f"Processing cold indicators for crop: {crop}")

    # Crop data for the target crop only
    cropdat = rdata[rdata_key_map[crop]]

    # Keep only grid cells where this crop is actually grown (rain or irr area > 0)
    cropdat = cropdat[(cropdat["rain_area"] > 0) | (cropdat["irr_area"] > 0)].copy()
    cropdat = cropdat.drop_duplicates(subset=["lon", "lat"])

    print(f"  Grid cells with {crop} area: {len(cropdat)}")

    # Growing season calendars for the target crop only
    season_firr  = xr.open_dataset(cal_base / f"ggcmi-crop-calendar-phase3_2015soc_{crop}_firr_landonly.nc")
    season_noirr = xr.open_dataset(cal_base / f"ggcmi-crop-calendar-phase3_2015soc_{crop}_noirr_landonly.nc")

    ds_cold = season_stat(temp_days, season_firr, season_noirr, cropdat)

    out_path = out_dir / f"cold_indicators_{crop}.nc"
    ds_cold.to_netcdf(out_path)
    print(f"Saved: {out_path}")

    season_firr.close()
    season_noirr.close()

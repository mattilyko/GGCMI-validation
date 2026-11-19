# Raw input data: ISIMIP3a crop yield simulations

This folder is intended to hold the raw gridded crop model output data from the ISIMIP3a simulation project, which are required to reproduce the analysis in this repository.

These data files are **not included** in the GitHub repository. Instead, users must download the data manually from the official ISIMIP data portal.

## Data access instructions 

Please access the required data files here:
[https://data.isimip.org/search/tree/ISIMIP3a/OutputData/agriculture/](https://data.isimip.org/search/tree/ISIMIP3a/OutputData/agriculture/)

All files are available in **NetCDF format** for a resolution of 0.5°. 

1. Select all available crop models
2. Select the correct climate input data: **GSWP3-W5E5**
3. Select the correct experiment specifier: **obsclim**
4. Select the correct output variable: **yield**
5. Select the correct crops: **maize, rice (first growing season), rice (second growing season), soy, spring wheat, winter wheat**
6. Per specific model and crop there should be two distinct data files (one with **firr** in the name, the other with **noirr**). Select both. 
7. Click on **configure download**
8. Make sure all required files are selected. Note that not all crops are available for each model. 
9. Select **Mask only land data** in the MASK AREA box
10. Click on **download file** which creates a zip file with all data in your downloads
11. After downloading, organize the files into the designated folder structure (one subfolder per crop) inside the repository at `data/raw/GGCMI_yields`:

## File naming convention

Each file name follows the standardized ISIMIP naming scheme:

`<model>_<climate>_<experiment>_<socioeconomic>_<management>_yield-<crop>-<irrigation>_global_annual-gs_<start>_<end>.nc`

Underscores separate different specifiers of the corresponding simulation run

e.g. `lpjml_gswp3-w5e5_obsclim_2015soc_default_yield-mai-firr_global_annual-gs_1901_2016.nc`

This indicates:
- `lpjml`: GGCMI crop model.
- `gswp3-w5e5`: climate forcing input.
- `obsclim`: observational-based climate related forcing.
- `2015soc`: fixed 2015-level direct human influences (e.g. land use and fertilizers).
- `default`: For all experiments other than the sensitivity experiments.
- `yield-mai-firr`: maize yield with full irrigation.
- `global_annual-gs`: global annual output based on growing-season.
- `1901_2016`: simulation time span.

For full ISIMIP3a simulation experiment details and explanation of the specifiers, see the [ISIMIP3a protocol](https://protocol.isimip.org/#/ISIMIP3a/agriculture) 
and the corresponding paper: [Frieler et al. (2024)](https://doi.org/10.5194/gmd-17-1-2024)

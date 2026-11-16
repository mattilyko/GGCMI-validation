# GGCMI dataframes

This folder contains the restructed R dataframes extracted from the calendar adjusted NetCDF files of the GGCMI yield model runs. `02_ISIMIP3a_dataprep.qmd` is the notebook where these files can be created. 

Each RData file in the folder contains several dataframes for the specific model and crop combination: 
- **yld_long_irr**: irrigated yields per gridcell
- **yld_long_rain**: rainfed yields per gridcell
- **cntr_prod_irr**: national production values for irrigated
- **cntr_prod_rain**: national production values for rainfed
- **cntr_prod**: total national production
- **area_irr**: irrigated crop area per country
- **area_rain**: rainfed crop area per country
- **area**: total crop area per country
- **yield**: total national yields

The rows in the dataframes represent different locations (gridcells or countries), the columns represent different years. 

This folder is structured per crop. 

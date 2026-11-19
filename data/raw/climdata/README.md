# Climate Data

This folder is intended to store the ISIMIP3a climate forcing dataset **GSWP3-W5E5** (precipitation and tasmax), used to compute climate extremes indicators in our study. 

The data is **not included** in the GitHub repository. To reproduce the results of the analysis the required data files need to be downloaded manually. 

## Download Instructions

The required files can be downloaded from the ISIMIP repository: [https://data.isimip.org/search/tree/ISIMIP3a/](https://data.isimip.org/search/tree/ISIMIP3a/InputData/)

### Steps:

1. Select first the [tasmax data](https://data.isimip.org/search/tree/ISIMIP3a/InputData/climate/atmosphere/gswp3-w5e5/obsclim/tasmax/).
2. Then also select the [precipitation data](https://data.isimip.org/search/tree/ISIMIP3a/InputData/climate/atmosphere/gswp3-w5e5/obsclim/pr/).
3. Click on **configure download**
4. Select **Mask only land data** in the MASK AREA box
5. Click on **download file** which creates a zip file with all data in your downloads
6. After downloading, organize the files inside the repository at `data/raw/climdata`
   
## Dataset Description

The dataset describes daily (maximum)temperature/precipitation reanalysis data from 1901 to 2019. Outputs are always provided in **NetCDF format** with a resolution of 0.5°. The same data were used as climate forcing for the ISIMIP3a simulation runs analysed in this analysis (obsclim). More details on the dataset can be found on the [ISIMIP website](https://www.isimip.org/gettingstarted/input-data-bias-adjustment/details/110/). 


# Other required data files

This folder is intended to store additional ISIMIP3a input data used to process the crop data in our study, including: **a country mask**, **a land-use dataset** and **crop-specific growing season calendars**. 

The land-use and crop calendar data are **not included** in the GitHub repository. To reproduce the results of the analysis the required data files need to be downloaded manually. 

## Download Instructions

The required files can be downloaded from the ISIMIP repository: [https://data.isimip.org/search/tree/ISIMIP3a/](https://data.isimip.org/search/tree/ISIMIP3a/InputData/)

### Steps:

1. Select first the [land-use data](https://data.isimip.org/search/tree/ISIMIP3a/InputData/socioeconomic/landuse/2015soc/landuse-15crops/).
2. Then also select the [crop calendars](https://data.isimip.org/search/tree/ISIMIP3a/InputData/socioeconomic/crop_calendar/2015soc/ggcmi-crop-calendar-phase3/). You only need to select the calendars for **mai, ri1, ri2, wwh, swh and soy**. Note that for each crop there are two distinct calendars needed, one for full irrigation simulations (**firr**) and one for no irrigation runs (**noirr**).  
3. Click on **configure download**
4. Select **Mask only land data** in the MASK AREA box
5. Click on **download file** which creates a zip file with all data in your downloads
6. After downloading, organize the files inside the repository at `data/raw/other`
   
## Dataset Description

Outputs are always provided in **NetCDF format** with a resolution of 0.5° on a global scale.

The land-use data describes gridded harvest areas for the corresponding crops and irrigation management fixed for the year of 2015. We use these data to integrate the two different simulations (noirr and firr) resulting in one yield measure taking into account the relative share of irrigated and no-irrigated harvest areas for each crop. For more information on the land-use data we refer to the [ISIMIP website](https://data.isimip.org/10.48364/ISIMIP.571261.3).

The crop calendar data describes for the corresponding crops and irrigation management the typical plantation and harvesting dates for each gridcell fixed for the year of 2015. We use these data for the computation of the climate extremes indicators during the corresponding growing season. For more information on the crop calendar data we refer to [zenodo](https://zenodo.org/records/5062513).



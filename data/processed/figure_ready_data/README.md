# Figure ready data

 In this folder all the data files ready for the compilation of the figures are stored. The data files are computed for each crop, including also the crop-aggregated data (aggr). The code generating these files can be found in `code/analysis/filtering_extremes.qmd`. 

 
The folder contains:
- `extremes_<crop>`: data for analysing performance under climatic extremes (filtered for extremes)
- `general_<crop>`: data for analysing general performance for all gridcell-years (not filtered for extremes, but prefiltered for relevant gridcells). For all crops these contain both simulated (including ensemble median) and benchmark data, for aggr this only contains benchmark data. 
- `general_sim_aggr`: simulation data (including ensemble median) for analyzing general performance. 

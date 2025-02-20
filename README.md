## DWS_data_aggregations
The repository contains scripts for postprocessing the data which is provided as part of the NWO/ENW project: LOCO-EX (OCENW.KLEIN.138). The data which includes the Dutch Wadden Sea (DWS) related variables and is stored in NetCDF data format.

*Postprocessing includes:*

- Creation of mask to estabilish the area of Dutch Wadden Sea in [create_boundary.py](https://github.com/LOCO-EX/Spatial_and_15days_aggregations/blob/main/data_processing_scripts/create_boundary.py).
- Filtering so as to exclude the dry areas based on the bathymetry and sea surface elevation difference for given treshold in the file [depth_process_data.py](https://github.com/LOCO-EX/Spatial_and_15days_aggregations/blob/main/data_processing_scripts/depth_process_data.py). The script creates only temporary files, which are later used in 15-day and spatial aggregations.
- Spatial aggregations (average and standard deviation) of salinity and temperature variables for the whole area of Dutch Wadden Sea each hour over a given date range in the file [spatial_aggregates.py](https://github.com/LOCO-EX/Spatial_and_15days_aggregations/blob/main/data_processing_scripts/spatial_aggregates.py).
- 15-day aggregations (average and standard deviation) of salinity and temperature variables and creation of expousure percentage for 15-day periods in the file [15day_aggregates.py](https://github.com/LOCO-EX/DWS_data_aggregations/blob/main/data_processing_scripts/15day_aggregates.py).
- Calculation of the volume of the Dutch Wadden Sea each hour over a given date range in the file [timeseries_dws_volume.py](https://github.com/LOCO-EX/Spatial_and_15days_aggregations/blob/main/data_processing_scripts/timeseries_dws_volume.py).
- The volume flux from rivers into the Dutch Wadden Sea five times an hour in the file [timeseries_rivers_flux.py](https://github.com/LOCO-EX/Spatial_and_15days_aggregations/blob/main/data_processing_scripts/timeseries_rivers_flux.py).
- The volume flux and the salinity flux trough predetermined transects in the file [timeseries_transects_flux.py](https://github.com/LOCO-EX/DWS_data_aggregations/blob/main/data_processing_scripts/timeseries_transects_flux.py).


### Software
The environment employed for the analysis is based on Python v3.12 and can be found in the file [environment.yml](https://github.com/LOCO-EX/Spatial_and_15days_aggregations/blob/main/environment.yml).


### Information about raw numerical data

The netCDF which were postprocessed are provided as part of the NWO/ENW project: LOCO-EX (OCENW.KLEIN.138). The numerical simulations were done thanks to the North-German Supercomputing Alliance (HLRN).

### Information about postprocessed data

The doi to location of postprocessed data at https://data.4tu.nl/ will be provided.

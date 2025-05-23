#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import date, datetime, timezone
from pathlib import Path
from time import time

import numpy as np
import xarray as xr
from numpy.typing import NDArray

PATH_ROOT = Path(__file__).parent
START_DATE = date()
END_DATE = date()


def timeseries_dws_volume(
    date_start: date,
    date_end: date,
    save=False,
    path_boundary: str | Path = "",
    path_combined: str | Path = "",
    path_output: str | Path = "",
):
    """Creates a netCDF with the volume of the Dutch Wadden Sea each hour over a given date range.

    Parameters
    ----------
    date_start : date
        The date from which to consider the volume; is not sensitive to different days
    date_end : date
        The date until which to consider the volume; is not sensitive to different days
    save : bool, optional
        Whether to save the file, by default False
    path_boundary : str | Path, optional
        The relative path to the boundary mask file, by default ""
    path_combined : str | Path, optional
        The relative path to the combined data directory, by default ""
    path_output : str | Path, optional
        The relative path to the output directory, by default ""
    """

    # Get DWS area mask
    if not Path(PATH_ROOT / path_boundary).is_file():
        print(f"{PATH_ROOT / path_boundary} does not exist")
        quit()
    ds_bounds = xr.open_dataset(PATH_ROOT / path_boundary)
    mask_dws: NDArray[np.bool] = ~ds_bounds["mask_dws"].values.copy()
    ds_bounds.close()

    # Initialise DataArray list
    das_volume: list[xr.DataArray] = list()

    for year in np.arange(date_start.year, date_end.year + 1, 1):
        for month in np.arange(date_start.month, date_end.month + 1, 1):
            print(f"Started {year}{str(month).zfill(2)}01")
            time_start = time()

            # File path
            path_data_combined = (
                PATH_ROOT
                / path_combined
                / f"RE.TS.vertmean.{year}{str(month).zfill(2)}{"01"}.treshold0.15.nc"
            )
            if not Path(path_data_combined).is_file():
                print(f"{path_data_combined} does not exist")
                quit()

            # Open trace file
            ds_trace_date = xr.open_dataset(path_data_combined)

            # Apply DWS mask
            timesteps = ds_trace_date["time"].shape[0]
            for t in np.arange(timesteps):
                ds_trace_date["water_depth"].values[t, mask_dws] = np.nan
            print(f"DataSet is masked")

            # Create volume DataArray
            da_volume = ds_trace_date["water_depth"].sum(dim=["xc", "yc"]) * (
                200**2
            )  # Multipy by area
            das_volume.append(da_volume)
            print(f"volume is added")

            print(f"Took {time()-time_start} s")

            print("")

        ds_volume = xr.concat(das_volume, dim="time").rename("volume").to_dataset()

        # Add meta data
        ds_volume["volume"].attrs = {
            "units": "m3",
            "cell_methods": "area: sum",
            "long_name": "volume",
        }
        ds_volume["time"].attrs = {
            "standard_name": "time",
            "long_name": "time of measurement",
        }
        ds_volume.attrs = {
            "title": "Volume of the Dutch Wadden Sea per hour.",
            "Conventions": "CF-1.12",
            "institution": "www.tue.nl; www.nioz.nl; www.io-warnemuende.de",
            "email": "m.duran.matute@tue.nl; theo.gerkema@nioz.nl; ulf.graewe@io-warnemuende.de",
            "source": "GETM (www.getm.eu)",
            "comment": "This data is provided as part of the NWO/ENW project: LOCO-EX (OCENW.KLEIN.138). The numerical simulations were done thanks to the North-German Supercomputing Alliance (HLRN).",
            "history": f"Created {datetime.now().replace(tzinfo=timezone.utc).isoformat(timespec="minutes")}.",
        }

        # Define the encoding for the dataset
        encoding_format = {
            "zlib": True,
            "complevel": 4,
            "shuffle": True,
        }
        encoding = {
            "volume": encoding_format,
            "time": encoding_format,
        }

        # Save the dataset to a new file
        if save:
            ds_volume.to_netcdf(
                (PATH_ROOT / path_output / f"DWS.volume.nc").resolve(),
                "w",
                format="NETCDF4",
                encoding=encoding,
            )
            print("Saved netCDF file")


if __name__ == "__main__":
    start_date = START_DATE
    end_date = END_DATE
    timeseries_dws_volume(
        start_date,
        end_date,
        save=True,
    )

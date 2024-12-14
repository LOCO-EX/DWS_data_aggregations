#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import date, datetime, timezone
from pathlib import Path
from time import time

import numpy as np
import xarray as xr
from numpy.typing import NDArray

PATH_ROOT = Path(__file__).parent


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
                PATH_ROOT / path_combined / f"RE.DWS200m.combined.{year}{str(month).zfill(2)}{"01"}.nc"
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
            da_volume = ds_trace_date["water_depth"].sum(dim=["xc", "yc"]) * (200**2)  # Multipy by area
            das_volume.append(da_volume)
            print(f"volume is added")

            print(f"Took {time()-time_start} s")

            print("")

        ds_volume = xr.concat(das_volume, dim="time").rename("volume").to_dataset()

        ds_volume["volume"].attrs = {
            "units": "m3",
            "cell_methods": "area: sum",
            "long_name": "volume",
        }
        ds_volume.attrs = {
            "title": "Volume of the Dutch Wadden Sea per hour.",
            "Conventions": "CF-1.12",
            "institution": "www.tue.nl; www.nioz.nl",
            "email": "m.duran.matute@tue.nl; theo.gerkema@nioz.nl",
            "source": "GETM",
            "comment": "Produced as part of the LOCO-ex project",
            "history": f"Created {datetime.now().replace(tzinfo=timezone.utc).isoformat(timespec="minutes")} using timeseries_dws_volume.py",
        }

        if save:
            ds_volume.to_netcdf(
                PATH_ROOT
                / path_output
                / f"DWS.volume.{date_start.strftime("%Y%m")}01-{date_end.strftime("%Y%m")}01.nc",
                "w",
                format="NETCDF4",
            )
            print("Saved netCDF file")


if __name__ == "__main__":
    timeseries_dws_volume(
        date(2000, 1, 1),
        date(2000, 1, 1),
        save=True,
    )

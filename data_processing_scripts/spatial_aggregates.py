#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import date, datetime, timezone
from pathlib import Path
from time import time

import numpy as np
import xarray as xr
from dateutil.relativedelta import relativedelta
from numpy.typing import NDArray

PATH_ROOT = Path(__file__).parent


def spatial_aggregates(
    date_start: date,
    date_end: date,
    variable_to_summarise: str,
    save=False,
    path_boundary: str | Path = "",
    path_combined: str | Path = "",
    path_aggregates: str | Path = "",
):
    """Creates a netCDF with the spatial mean and standard deviation of a given variable of the Dutch Wadden Sea each hour over a given date range.

    Parameters
    ----------
    date_start : date
        The date from which to consider the volume; is not sensitive to different days
    date_end : date
        The date until which to consider the volume; is not sensitive to different days
    variable_to_summarise : str
        The name of the variable in the netCDF file to aggregate
    save : bool, optional
        Whether to save the file, by default False
    path_boundary : str | Path, optional
        The relative path to the boundary mask file, by default ""
    path_combined : str | Path, optional
        The relative path to the combined data directory, by default ""
    path_aggregates : str | Path, optional
        The relative path to the aggregate output directory, by default ""
    """

    # Get DWS area mask
    if not Path(PATH_ROOT / path_boundary).is_file():
        print(f"{PATH_ROOT / path_boundary} does not exist")
        quit()
    ds_bounds = xr.open_dataset(PATH_ROOT / path_boundary)
    mask_dws: NDArray[np.bool] = ~ds_bounds["mask_dws"].values.copy()
    ds_bounds.close()

    # Initialise aggregates list
    means: list[np.float32] = list()
    stds: list[np.float32] = list()

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
                ds_trace_date[variable_to_summarise].values[t, mask_dws] = np.nan
            print(f"DataSet is masked")

            # Add total_volume variable per time
            ds_trace_date = ds_trace_date.assign(
                total_volume=ds_trace_date["water_depth"].sum(dim=["xc", "yc"])
            )
            print(f"total_volume is added")

            for t in np.arange(timesteps):
                # Calculate weights; water_depth (local volume) / total_volume
                weights = (
                    ds_trace_date["water_depth"].values[t, :, :] / ds_trace_date["total_volume"].values[t]
                )

                # Calculate weighted mean and std
                mean = np.nansum(ds_trace_date[variable_to_summarise].values[t, :, :] * weights)
                std = np.sqrt(
                    np.nansum((ds_trace_date[variable_to_summarise].values[t, :, :] - mean) ** 2 * weights)
                )

                means.append(mean)
                stds.append(std)

            print(f"Variable is summarised")
            print(f"Took {time()-time_start} s")

            print("")

    ds = create_ds(means, stds, date_start, date_end, variable_to_summarise)
    print("Created DataSet")

    if save:
        ds.to_netcdf(
            PATH_ROOT
            / path_aggregates
            / f"DWS200m.aggregates.{variable_to_summarise}.{date_start.strftime("%Y%m%d")}-{date_end.strftime("%Y%m%d")}.nc",
            "w",
            format="NETCDF4",
        )
        print("Saved netCDF file")


def create_ds(
    means: list[np.float32], stds: list[np.float32], date_start: date, date_end: date, variable: str
):
    # End date with one month added to get a date range
    date_end_inclusive = date_end + relativedelta(months=1)

    match variable:
        case "T":
            var_mean_name = "T_mean"
            var_std_name = "T_std"
            var_units = "degC"
            var_name = "temperature"
        case "S":
            var_mean_name = "S_mean"
            var_std_name = "S_std"
            var_units = "g kg-1"
            var_name = "salinity"
        case _:
            print("Error")
            quit()

    data_vars_dict = {
        var_mean_name: (
            ["time"],
            means,
            {
                "units": var_units,
                "units_metadata": f"{var_name}: on_scale",
                "long_name": f"spatially averaged {var_name}",
                "cell_methods": "area: mean",
            },
        ),
        var_std_name: (
            ["time"],
            stds,
            {
                "units": var_units,
                "units_metadata": f"{var_name}: difference",
                "long_name": f"spatial standard deviation {var_name}",
                "cell_methods": "area: standard_deviation",
            },
        ),
    }
    coords_dict = dict(
        time=(
            ["time"],
            np.arange(
                np.datetime64(f"{date_start.year}-{str(date_start.month).zfill(2)}-01T01"),
                np.datetime64(f"{date_end_inclusive.year}-{str(date_end_inclusive.month).zfill(2)}-01T01"),
                np.timedelta64(1, "h"),
                dtype="datetime64[ns]",
            ),
            {"standard_name": "time", "long_name": "time of measurement"},
        )
    )
    attrs_dict = dict(
        title=f"Mean and standard deviation of the {var_name} over the Dutch Wadden Sea at each hour.",
        Conventions="CF-1.12",
        institution="www.tue.nl; www.nioz.nl; www.io-warnemuende.de",
        email="m.duran.matute@tue.nl; theo.gerkema@nioz.nl; ulf.graewe@io-warnemuende.de",
        source="GETM (www.getm.eu)",
        comment="This data is provided as part of the NWO/ENW project: LOCO-EX (OCENW.KLEIN.138). The numerical simulations were done thanks to the North-German Supercomputing Alliance (HLRN).",
        history=f"Created {datetime.now().replace(tzinfo=timezone.utc).isoformat(timespec="minutes")} using spatial_aggragates.py",
    )
    ds = xr.Dataset(
        data_vars=data_vars_dict,
        coords=coords_dict,
        attrs=attrs_dict,
    )
    ds["time"].encoding["calendar"] = "standard"

    return ds


if __name__ == "__main__":
    spatial_aggregates(
        date(2000, 1, 1),
        date(2000, 1, 1),
        "S",
        save=True,
    )
    spatial_aggregates(
        date(2000, 1, 1),
        date(2000, 1, 1),
        "T",
        save=True,
    )

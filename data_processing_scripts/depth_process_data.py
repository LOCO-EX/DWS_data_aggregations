#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
from datetime import datetime, timezone

import numpy as np
import xarray as xr

# Global params
# Define the treshold to exclude dry land
TRESHOLD = 0
# Define the path to the files
PATH_DISC = ""  # args.file_tracer_uvz  # "D:\loco-ex\PP"
PATH_SAVE_DISC = ""  # args.file_tracer_depth  # "D:\loco-ex\dry_exclude"
# Indicate the first date to process
START_DATE = 0
# Indicate part of the names of the files with data
Z_FILE_NAME = ""
TRACER_FILE_NAME = ""
# Indicate how many years to process
YEARS = 0


# Function to process the files
def process_file_based_on_water_depth_and_treshold(
    uvz_path: str,
    tracer_path: str,
    diff_path: str,
    include_water_depth: bool,
    include_mask: bool,
):
    """
    Processes the given file, create a new one and save it.

    Parameters:
    uvz_path : str
        The path to the uvz file from which the sea surface elevation and bathymetry will be extracted.
    tracer_path : str
        The path to the tracer file to be processed.
    diff_path : str
        The destination path where the copy of tracer file with added differences between sea surface elevation and bathymetry,
        dropped variables ['Denoever', 'Kornwerderzand'], and modified salinity and temperature variables
        (modification based on the differences between sea surface elevation and bathymetry and TRESHOLD) will be saved.
    include_water_depth : bool
        Whether to include the water depth as variable in the output netCDF file
    include_mask : bool
        Whether to include the mask as variable in the output netCDF file
    """
    # Open the uvz file where bathymetry and sea surface elevation are stored
    ds_uvz = xr.open_dataset(uvz_path)

    # Open the tracer file
    ds_tracer = xr.open_dataset(tracer_path)

    # Read the sea surface elevation and bathymetry
    sea_level = ds_uvz["z"].values  # Shape: (time, yc, xc)
    bathymetry = ds_uvz["bathymetry"].values  # Shape: (yc, xc)

    # Read the xc, yc and time as DataArrays
    da_xc, da_yc, da_time = ds_uvz["xc"], ds_uvz["yc"], ds_uvz["time"]

    # Select only the full hours
    time_minutes = da_time.dt.minute
    full_hour_indices = np.where(time_minutes == 0)[0]
    sea_level = sea_level[full_hour_indices, :, :]

    # Open the uvz file to check if the time is the same, readable only if decode_times=False
    ds_uvz_ = xr.open_dataset(uvz_path, decode_times=False)
    ds_tracer_ = xr.open_dataset(tracer_path, decode_times=False)
    # Check if the files have the same start time
    if ds_uvz_["time"].attrs["units"] != ds_tracer_["time"].attrs["units"]:
        raise ValueError(
            f"The files {uvz_path} and {tracer_path} do not have the same start time. Please, check the files."
        )
    ds_tracer_.close()

    # Compute the difference between the sea surface elevation and the bathymetry
    water_depth = bathymetry + sea_level

    # Create a mask for points where the bathymetry is nan or -10.0 (missing_value/FillValue)
    mask_bathymetry = np.logical_or((bathymetry == -10.0), (np.isnan(bathymetry)))
    # Create a mask for points where the depth is < 0.2 (chosen treshold) or the sea level is -9999 (missing_value/FillValue)
    mask_over_time = np.logical_or((water_depth < TRESHOLD), (sea_level == -9999.0))
    # Combine the masks
    mask = np.empty_like(mask_over_time)
    for t in np.arange(sea_level.shape[0]):
        mask[t, :, :] = np.logical_or(mask_over_time[t, :, :], mask_bathymetry)

    # Read the S and T variables
    S = ds_tracer["salt"].values
    T = ds_tracer["temp"].values

    time_minutes_ST = ds_tracer["time"].dt.minute
    full_hour_indices_ST = np.where(time_minutes_ST == 0)[0]
    S = S[full_hour_indices_ST, :, :]
    T = T[full_hour_indices_ST, :, :]

    # Replace S, T and diff variables with nan at the mask
    S[mask] = np.nan
    T[mask] = np.nan
    water_depth[mask] = np.nan

    # Create new DataSet
    ds = xr.Dataset(
        data_vars={
            "S": (
                ["time", "yc", "xc"],
                S,
                {
                    "units": "g kg-1",
                    "coordinates": "xc yc",
                    "long_name": "depth-avg. salinity",
                    "standard_name": "sea_water_absolute_salinity",
                },
            ),
            "T": (
                ["time", "yc", "xc"],
                T,
                {
                    "units": "degC",
                    "units_metadata": "temperature: on_scale",
                    "coordinates": "xc yc",
                    "long_name": "depth-avg. temperature",
                    "standard_name": "sea_water_temperature",
                },
            ),
            "lonc": (
                ["yc", "xc"],
                ds_uvz["lonc"].values,
                {"units": "degrees_east", "long_name": "longitude"},
            ),
            "latc": (
                ["yc", "xc"],
                ds_uvz["latc"].values,
                {"units": "degrees_north", "long_name": "latitude"},
            ),
        },
        coords=dict(
            # time=(
            #     ["time"],
            #     da_time[full_hour_indices],
            #     {"standard_name": "time", "long_name": "time of measurement"},
            # ),
            time=ds_uvz_["time"],
            xc=da_xc,
            yc=da_yc,
        ),
        attrs=dict(
            title="Dutch Wadden Sea - 200 m resolution : vUERRAL02",
            Conventions="CF-1.12",
            institution="www.tue.nl; www.nioz.nl; www.io-warnemuende.de",
            email="m.duran.matute@tue.nl; theo.gerkema@nioz.nl; ulf.graewe@io-warnemuende.de",
            source="GETM (www.getm.eu)",
            comment="This data is provided as part of the NWO/ENW project: LOCO-EX (OCENW.KLEIN.138). The numerical simulations were done thanks to the North-German Supercomputing Alliance (HLRN).",
            history=f"Created {datetime.now().replace(tzinfo=timezone.utc).isoformat(timespec="minutes")} using depth_process_data.py",
        ),
    )
    ds_uvz_.close()

    if include_water_depth:
        ds["water_depth"] = (
            ("time", "yc", "xc"),
            water_depth,
            {
                "units": "m",
                "coordinates": "xc yc",
                "cell_methods": "time: point",
                "long_name": "difference between sea surface elevation and bathymetry",
            },
        )
    if include_mask:
        ds["mask"] = (("time", "yc", "xc"), mask)
        # ds["land"].attrs = {"long_name": "dry land"}

    # Define the encoding for the dataset
    encoding_format = {
        "zlib": True,
        "complevel": 4,
        "shuffle": True,
    }
    # Specify compression for the variable
    encoding = {
        "S": encoding_format,
        "T": encoding_format,
    }
    if include_water_depth:
        encoding["water_depth"] = encoding_format
    if include_mask:
        encoding["mask"] = encoding_format

    ds.to_netcdf(diff_path, format="NETCDF4", encoding=encoding)
    print(ds)

    print("File processed successfully.")

    ds_uvz.close()
    ds_tracer.close()


def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Process a .nc file.")

    parser.add_argument(
        "-d",
        "--water_depth",
        action="store_true",
        help="Include the water depth in the netCDF file",
    )
    parser.add_argument(
        "-m", "--mask", action="store_true", help="Include the mask in the netCDF file"
    )

    args = parser.parse_args()
    print("Treshold: ", TRESHOLD)

    for j in range(0, YEARS):
        date = START_DATE + j * 10000
        for i in range(0, 12):  # for processing the whole year
            print(date + i * 100)
            process_file_based_on_water_depth_and_treshold(
                os.path.join(
                    PATH_DISC + "\\" + str(date + i * 100),
                    Z_FILE_NAME + str(date + i * 100) + ".nc",
                ),
                os.path.join(
                    PATH_DISC + "\\" + str(date + i * 100),
                    TRACER_FILE_NAME + str(date + i * 100) + ".nc",
                ),
                os.path.join(
                    PATH_SAVE_DISC,
                    f"{TRACER_FILE_NAME + str(date + i * 100)}.treshold{TRESHOLD}.nc",
                ),
                args.water_depth,
                args.mask,
            )


if __name__ == "__main__":
    main()

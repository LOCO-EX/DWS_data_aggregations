#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from datetime import datetime, timezone
from netCDF4 import date2num

import numpy as np
import xarray as xr
import os


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
    ds_tracer = xr.open_dataset(tracer_path, decode_times=False)

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
    ds_uvz_units = xr.open_dataset(uvz_path, decode_times=False)
    # Check if the files have the same start time
    if (ds_uvz_units["time"].attrs["units"] != ds_tracer["time"].attrs["units"]):
        raise ValueError(
            f"The files {uvz_path} and {tracer_path} do not have the same start time. Please, check the files."
        )
    ds_uvz_units.close()

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
    S = ds_tracer["S"].values
    T = ds_tracer["T"].values

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
        },
        coords=dict(
            # time=(
            #     ["time"],
            #     da_time[full_hour_indices],
            #     {"standard_name": "time", "long_name": "time of measurement"},
            # ),
            time=ds_tracer["time"],
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
        "file_uvz",
        type=str,
        help="Path to the folder with the files from which the sea level and bathymetry will be extracted.",
    )
    parser.add_argument(
        "file_tracer",
        type=str,
        help="Path to the folder with the files to be processed.",
    )
    parser.add_argument(
        "file_tracer_depth",
        type=str,
        help="Path to the folder with the files where to save the modified copies of the files from file_tracer.",
    )
    parser.add_argument(
        "-t",
        "--treshold",
        type=float,
        default=0.15,
        help="The treshold for the difference between sea surface elevation and bathymetry",
    )
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

    # Assign treshold value, default 15.0
    global TRESHOLD
    TRESHOLD = args.treshold
    print("Treshold: ", TRESHOLD)

    # Get list of files from the directories and sort them alphabetically
    files_in_uvz = os.listdir(args.file_uvz)
    files_in_uvz_sorted = sorted(files_in_uvz)

    files_in_tracer = os.listdir(args.file_tracer)
    files_in_tracer_sorted = sorted(files_in_tracer)

    for i in range(len(files_in_uvz_sorted)):
        # Process the file
        print(f"Input files: {files_in_uvz_sorted[i]}, {files_in_tracer_sorted[i]}")
        process_file_based_on_water_depth_and_treshold(
            os.path.join(args.file_uvz, files_in_uvz_sorted[i]),
            os.path.join(args.file_tracer, files_in_tracer_sorted[i]),
            os.path.join(
                args.file_tracer_depth,
                f"{files_in_tracer_sorted[i][:-3]}.treshold{TRESHOLD}.nc",
            ),
            args.water_depth,
            args.mask,
        )


if __name__ == "__main__":
    main()

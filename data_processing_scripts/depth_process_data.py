#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from datetime import datetime, timezone

import numpy as np
import xarray as xr

# Define the threshold for the difference between sea surface elevation and bathymetry
TRESHOLD = 0.15


def process_file_based_on_water_depth_and_treshold(
    uvz_path: str, tracer_path: str, diff_path: str, include_water_depth: bool, include_mask: bool
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

    # Read the sea surface elevation and bathymetry
    sea_level = ds_uvz["z"].values  # Shape: (time, yc, xc)
    bathymetry = ds_uvz["bathymetry"].values  # Shape: (yc, xc)

    # Read the xc, yc and time as DataArrays
    da_xc, da_yc, da_time = ds_uvz["xc"], ds_uvz["yc"], ds_uvz["time"]

    # Select only the full hours
    time_minutes = da_time.dt.minute
    full_hour_indices = np.where(time_minutes == 0)[0]
    sea_level = sea_level[full_hour_indices, :, :]

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

    # Open the tracer file
    ds_tracer = xr.open_dataset(tracer_path)

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
            time=(
                ["time"],
                da_time[full_hour_indices],
                {"standard_name": "time", "long_name": "time of measurement"},
            ),
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
        ds["land"].attrs = {"long_name": "dry land"}

    # Specify compression for the variable
    encoding = {
        "S": {
            "zlib": True,  # Enable compression
            "complevel": 4,  # Compression level (1-9, where 9 is maximum)
            "shuffle": True,  # Apply the shuffle filter (helps with compression)
        },
        "T": {
            "zlib": True,
            "complevel": 4,
            "shuffle": True,
        },
    }
    if include_water_depth:
        encoding["water_depth"] = {
            "zlib": True,
            "complevel": 4,
            "shuffle": True,
        }
    if include_mask:
        encoding["mask"] = {
            "zlib": True,
            "complevel": 4,
            "shuffle": True,
        }

    ds.to_netcdf(diff_path, format="NETCDF4", encoding=encoding)
    print(ds)

    print("File processed successfully.")


def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Process a .nc file.")
    parser.add_argument(
        "file_uvz",
        type=str,
        help="Path to the file from which the sea level and bathymetry will be extracted.",
    )
    parser.add_argument("file_tracer", type=str, help="Path to the file to be processed.")
    parser.add_argument(
        "file_tracer_diff",
        type=str,
        help="Path to the file where to save the modified copy of the file_tracer.",
    )
    parser.add_argument(
        "-d", "--water_depth", action="store_true", help="Include the water depth in the netCDF file"
    )
    parser.add_argument("-m", "--mask", action="store_true", help="Include the mask in the netCDF file")
    args = parser.parse_args()

    # Process the file
    process_file_based_on_water_depth_and_treshold(
        args.file_uvz, args.file_tracer, args.file_tracer_diff, args.water_depth, args.mask
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
from datetime import datetime, timedelta, timezone

import numpy as np
import xarray as xr
from netCDF4 import date2num

PATH_ROOT = ""
# Dates only to name the file properly - it doesn't affect the processing
START_DATE = 0
END_DATE = 0


def process_files(
    data_dir: str | str = "",
    processed_data_file: str | str = "",
):
    """
    Processes the files from the given directory.
    It calculates the 15-day average and standard deviation of the salinity and temperature.
    It also calculates the ratio of dry measurements to all measurements over a 15-day period.

    Parameters:
    data_dir : str
        The path to the directory where the data to process is stored.
    processed_data_file : str
        The path to the file where the processed files will be stored.
    """
    data_dir = PATH_ROOT + data_dir
    # processed_data_file = (
    #     PATH_ROOT
    #     + processed_data_file
    #     + "15day.aggregates."
    #     + str(FIRST_DATE)
    #     + "-"
    #     + str(END_DATE)
    #     + ".nc",
    # )
    processed_data_file = PATH_ROOT + processed_data_file
    print(processed_data_file)

    # Get list of files from the directory and sort them alphabetically
    files_in_dir1 = os.listdir(data_dir)
    files_in_dir1_sorted = sorted(files_in_dir1)

    # Create empty lists to store the data
    S_merge_ = []
    T_merge_ = []
    last_date = None
    S_15day_avg = None

    for file1 in files_in_dir1_sorted:

        # Load the file
        ds = xr.open_dataset(data_dir + "/" + file1)

        # The check of the end and start date from the previous and following files
        if last_date is not None:
            print(f"Delta: {ds['time'].values[0] - last_date}")
            # Check if the time is continuous between the files (1 hour difference)
            if (ds["time"].values[0] - last_date) != np.timedelta64(
                3600000000000, "ns"
            ):
                print(
                    f"Error: There is not a continous time between the {file1} and previously processed file."
                )
                break

        last_date = ds["time"].values[-1]

        S_merge_.extend(ds["S"].values)
        T_merge_.extend(ds["T"].values)

        # Close the file
        ds.close()

        # Convert the lists to numpy arrays
        S_merge = np.array(S_merge_)
        T_merge = np.array(T_merge_)

        # Get the dimensions of the data
        time_steps, yc_dim_, xc_dim_ = S_merge.shape
        # Define the number of steps per period
        steps_per_period = 15 * 24
        valid_steps = (time_steps // steps_per_period) * steps_per_period
        if valid_steps != 0:
            # Reshape the data to have the 15-day periods
            S_merge = S_merge[:valid_steps, :, :]
            S_reshaped = S_merge.reshape(-1, steps_per_period, yc_dim_, xc_dim_)

            T_merge = T_merge[:valid_steps, :, :]
            T_reshaped = T_merge.reshape(-1, steps_per_period, yc_dim_, xc_dim_)

            # Mask the NaN values
            masked_S = np.ma.masked_array(S_reshaped, np.isnan(S_reshaped))
            masked_T = np.ma.masked_array(T_reshaped, np.isnan(T_reshaped))

            # Calculate the 15-day average and standard deviation and don't consider NaN values
            if S_15day_avg is None:  # First iteration
                S_15day_avg = np.mean(masked_S, axis=1).filled(np.nan)
                T_15day_avg = np.mean(masked_T, axis=1).filled(np.nan)

                S_15day_sd = np.std(masked_S, axis=1).filled(np.nan)
                T_15day_sd = np.std(masked_T, axis=1).filled(np.nan)
                # Calculate the ratio of dry measurements to all measurements over a 15-day periods
                dry_measurement_ratio = (
                    np.isnan(S_reshaped).sum(axis=1) / S_reshaped.shape[1] * 100
                )

            else:  # Subsequent iterations
                S_15day_avg = np.concatenate(
                    (S_15day_avg, np.mean(masked_S, axis=1).filled(np.nan)), axis=0
                )
                T_15day_avg = np.concatenate(
                    (T_15day_avg, np.mean(masked_T, axis=1).filled(np.nan)), axis=0
                )

                S_15day_sd = np.concatenate(
                    (S_15day_sd, np.std(masked_S, axis=1).filled(np.nan)), axis=0
                )
                T_15day_sd = np.concatenate(
                    (T_15day_sd, np.std(masked_T, axis=1).filled(np.nan)), axis=0
                )
                # Calculate the ratio of dry measurements to all measurements over a 15-day periods
                dry_measurement_ratio = np.concatenate(
                    (
                        dry_measurement_ratio,
                        np.isnan(S_reshaped).sum(axis=1) / S_reshaped.shape[1] * 100,
                    ),
                    axis=0,
                )

            # Remove the processed 15-day periods from the lists
            S_merge_ = S_merge_[valid_steps:]
            T_merge_ = T_merge_[valid_steps:]

    # Convert double to float
    S_15day_avg = S_15day_avg.astype(np.float32)
    T_15day_avg = T_15day_avg.astype(np.float32)
    S_15day_sd = S_15day_sd.astype(np.float32)
    T_15day_sd = T_15day_sd.astype(np.float32)

    # Check if there are NaN values in the Wadden Sea area
    # - which means that there were dry points for the 15-day period
    area_values = S_15day_avg
    has_nan = np.isnan(area_values).sum()
    print(f"Number of NaN values in the area: {has_nan}")

    # Define start date based on the first data file
    ds_st = xr.open_dataset(data_dir + "/" + files_in_dir1[0])
    start_date_np = ds_st["time"].values[0]
    np_datetime_str = str(start_date_np)
    datetime_str_truncated = np_datetime_str[:26]
    start_date = datetime.strptime(datetime_str_truncated, "%Y-%m-%dT%H:%M:%S.%f")
    print(f"Start date: {start_date}")
    # Define the 15-day interval
    interval = timedelta(days=15)

    # Create the middle points for each 15-day period
    n_periods = T_15day_avg.shape[0]
    middle_times = [start_date + interval * i + interval / 2 for i in range(n_periods)]

    # Read data from the first data file to get the time units,
    # Use Dataset from netCDF4 library (not xarray - there is a problem in reading time units)
    ds_d = xr.open_dataset(data_dir + "/" + files_in_dir1[0], decode_times=False)
    # Create the time values for the dataset
    time = np.array(
        date2num(
            middle_times,
            units=ds_d["time"].attrs["units"],
            calendar=ds_d["time"].attrs["calendar"],
        ),
        dtype=np.float64,
    )

    # Read longtitude and latitude
    lonc = ds_d["lonc"].values
    latc = ds_d["latc"].values

    # Create a new file with the processed data
    create_file_to_save_processed_data(
        ds_st["xc"],
        ds_st["yc"],
        ds_d,
        time,
        S_15day_avg,
        T_15day_avg,
        S_15day_sd,
        T_15day_sd,
        dry_measurement_ratio,
        processed_data_file,
        latc,
        lonc,
    )
    ds_d.close()
    ds_st.close()


def create_file_to_save_processed_data(
    xc: np.array,
    yc: np.array,
    ds_d: xr.Dataset,
    time: np.array,
    S_avg: np.array,
    T_avg: np.array,
    S_sd: np.array,
    T_sd: np.array,
    dry_measurement_ratio: np.array,
    processed_data_file: str,
    latc: np.array,
    lonc: np.array,
):
    """
    Creates a new .nc file with the given data.

    Parameters:
    xc : np.array
        The x coordinates of the data.
    yc : np.array
        The y coordinates of the data.
    ds_d : xr.Dataset
        The dataset from the first data file.
    time : np.array
        The time values of the new data.
    S_avg : np.array
        The sea average salinity data.
    T_avg : np.array
        The temperature data.
    S_sd : np.array
        The standard deviation of the salinity data.
    T_sd : np.array
        The standard deviation of the temperature data.
    dry_measurement_ratio : np.array
        The ratio of dry measurements to all measurements over a 15-day period.
    processed_data_file : str
        The path to the file where the new data will be stored.
    latc : np.array
        The latitude data.
    lonc : np.array
        The longitude data.
    """
    # Get the current time to add to the history
    current_time = datetime.now().strftime("%a %b %d %H:%M:%S %Y")

    # Create a new dataset
    ds = xr.Dataset(
        {
            "S_avg": (
                ("time", "yc", "xc"),
                S_avg,
                create_variable_metadata(
                    "g kg-1", "salinity: on_scale", "15-day averaged salinity"
                ),
            ),
            "T_avg": (
                ("time", "yc", "xc"),
                T_avg,
                create_variable_metadata(
                    "degC", "temperature: on_scale", "15-day averaged temperature"
                ),
            ),
            "S_sd": (
                ("time", "yc", "xc"),
                S_sd,
                create_variable_metadata(
                    "g kg-1",
                    "salinity: on_scale",
                    "15-day standard deviation salinity",
                ),
            ),
            "T_sd": (
                ("time", "yc", "xc"),
                T_sd,
                create_variable_metadata(
                    "degC",
                    "temperature: on_scale",
                    "15-day standard deviation temperature",
                ),
            ),
            "exp_pct": (
                ("time", "yc", "xc"),
                dry_measurement_ratio,
                create_variable_metadata(
                    "%",
                    "%: on_scale",
                    "exposure percentage over a 15-day period",
                ),
            ),
            "lonc": (
                ("yc", "xc"),
                lonc,
                {
                    "units": "degrees_east",
                    "long_name": "longitude",
                },
            ),
            "latc": (
                ("yc", "xc"),
                latc,
                {
                    "units": "degrees_north",
                    "long_name": "latitude",
                },
            ),
        },
        coords={
            "xc": xc,
            "yc": yc,
            "time": (
                ["time"],
                time,
                {
                    "long_name": ds_d["time"].attrs["long_name"]
                    + " middle points of each 15-day period",
                    "units": ds_d["time"].attrs["units"],
                    "calendar": ds_d["time"].attrs["calendar"],
                },
            ),
        },
        attrs={
            "title": "Layers of 15-day aggregations of hydrodynamic quantities from the Dutch Wadden Sea",
            "conventions": "CF-1.12",
            "institution": "www.tue.nl; www.nioz.nl;www.io-warnemuende.de",
            "email": "m.duran.matute@tue.nl; theo.gerkema@nioz.nl; ulf.graewe@io-warnemuende.de",
            "source": "GETM (www.getm.eu)",
            "comment": "This data is provided as part of the NWO/ENW project: LOCO-EX (OCENW.KLEIN.138). The numerical simulations were done thanks to the North-German Supercomputing Alliance (HLRN).",
            "history": f"Created {datetime.now().replace(tzinfo=timezone.utc).isoformat(timespec='minutes')} using aggregates_15days.py",
        },
    )

    # Define the encoding for the dataset
    encoding_format = {
        "zlib": True,
        "complevel": 4,
        "shuffle": True,
    }
    encoding = {
        "S_avg": encoding_format,
        "T_avg": encoding_format,
        "S_sd": encoding_format,
        "T_sd": encoding_format,
        "exp_pct": encoding_format,
        "lonc": encoding_format,
        "latc": encoding_format,
    }

    # Save the dataset to a new file
    ds.to_netcdf(processed_data_file, format="NETCDF4", encoding=encoding)
    print("File with processed data created successfully.")


def create_variable_metadata(units, units_metadata, long_name):
    return {
        "_FillValue": np.nan,
        "units": units,
        "units_metadata": units_metadata,
        "long_name": long_name,
        "coordinates": "xc yc",
        "cell_methods": "time: point",
    }


def main():
    """
    This function processes .nc files with data gathered from Wadden Sea
    It takes input parameters for sea area boundaries, a directory of data
    to process, and an output path to the file where to save processed files.
    It uses argparse to parse command-line arguments. Once the arguments are parsed,
    it calls the `process_files` function to handle the core file processing logic.

    """

    # Process the file
    process_files()


if __name__ == "__main__":
    main()

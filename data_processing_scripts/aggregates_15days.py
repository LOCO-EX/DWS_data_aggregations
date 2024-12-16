import argparse
import numpy as np
import xarray as xr
import os
from datetime import datetime, timedelta
from netCDF4 import date2num, Dataset


def create_file_to_save_processed_data(
    xc: np.array,
    yc: np.array,
    units: str,
    time: np.array,
    S_avg: np.array,
    T_avg: np.array,
    S_sd: np.array,
    T_sd: np.array,
    dry_measurement_ratio: np.array,
    processed_data_file: str,
):
    """
    Creates a new .nc file with the given data.

    Parameters:
    xc : np.array
        The x coordinates of the data.
    yc : np.array
        The y coordinates of the data.
    units : str
        The units of the time data.
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
    """
    # Get the current time to add to the history
    current_time = datetime.now().strftime("%a %b %d %H:%M:%S %Y")

    # Create a new dataset
    ds = xr.Dataset(
        {
            "S_avg": (
                ("time", "yc", "xc"),
                S_avg,
                {
                    "_FillValue": np.nan,
                    "units": "g kg-1",
                    "units_metadata": "salinity: on_scale",
                    "long_name": "15-days averaged salinity",
                    "coordinates": "xc yc",
                    "cell_methods": "time: point",
                },
            ),
            "T_avg": (
                ("time", "yc", "xc"),
                T_avg,
                {
                    "_FillValue": np.nan,
                    "units": "degC",
                    "units_metadata": "temperature: on_scale",
                    "long_name": "15-days averaged temperature",
                    "coordinates": "xc yc",
                    "cell_methods": "time: point",
                },
            ),
            "S_sd": (
                ("time", "yc", "xc"),
                S_sd,
                {
                    "_FillValue": np.nan,
                    "units": "g kg-1",
                    "units_metadata": "salinity: on_scale",
                    "long_name": "15-days standard deviation salinity",
                    "coordinates": "xc yc",
                    "cell_methods": "time: point",
                },
            ),
            "T_sd": (
                ("time", "yc", "xc"),
                T_sd,
                {
                    "_FillValue": np.nan,
                    "units": "degC",
                    "units_metadata": "temperature: on_scale",
                    "long_name": "15-days standard deviation temperature",
                    "coordinates": "xc yc",
                    "cell_methods": "time: point",
                },
            ),
            "DW_ratio": (
                ("time", "yc", "xc"),
                dry_measurement_ratio,
                {
                    "_FillValue": np.nan,
                    "units": "%",
                    "units_metadata": "%: on_scale",
                    "long_name": "ratio of dry measurements to all measurements over a 15-day period",
                    "coordinates": "xc yc",
                    "cell_methods": "time: point",
                },
            ),
        },
        coords={
            "xc": xc,
            "yc": yc,
            "time": time,
        },
        attrs={
            "title": "Volume of the Dutch Wadden Sea per hour.",
            "conventions": "CF-1.12",
            "institution": "www.tue.nl; www.nioz.nl;www.io-warnemuende.de",
            "email": "m.duran.matute@tue.nl; theo.gerkema@nioz.nl; ulf.graewe@io-warnemuende.de ",
            "source": "GETM (www.getm.eu)",
            "comment": "This data is provided as part of the NWO/ENW project: LOCO-EX (OCENW.KLEIN.138). The numerical simulations were done thanks to the North-German Supercomputing Alliance (HLRN). ",
            "history": f"Created [{current_time}] using [{os.path.basename(__file__)}]",
        },
    )

    # # Add attributes to the time coordinate
    ds["time"].attrs["long_name"] = "time"
    ds["time"].attrs["units"] = units
    ds["time"].attrs["calendar"] = "standard"

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
        "DW_ratio": encoding_format,
    }

    # Save the dataset to a new file
    ds.to_netcdf(processed_data_file, format="NETCDF4", encoding=encoding)
    print("File with processed data created successfully.")


def process_files(dws_boundaries_area: str, data_dir: str, processed_data_file: str):
    """
    Processes the files from the given directory.
    It calculates the 15-day average and standard deviation of the salinity and temperature.
    It also calculates the ratio of dry measurements to all measurements over a 15-day period.

    Parameters:
    dws_boundaries_area : str
        The path to the file where the boundaries of the Wadden Sea area are stored.
    data_dir : str
        The path to the directory where the data to process is stored.
    processed_data_file : str
        The path to the file where the processed files will be stored.

    """

    # Load the boundaries and create a mask for the area of interest
    boundaries = xr.open_dataset(dws_boundaries_area)
    mask = boundaries["mask_dws"].values

    # Get list of files from the directory and sort them alphabetically
    files_in_dir1 = os.listdir(data_dir)
    files_in_dir1_sorted = sorted(files_in_dir1)

    # Create empty lists to store the data
    S_merge = []
    T_merge = []

    for file1 in files_in_dir1_sorted:

        # Load the file
        ds = xr.open_dataset(data_dir + "/" + file1)

        # Here should be inlcuded the check of the end and start date from the previous and following files

        # Prepare the mask
        expanded_mask = np.expand_dims(mask, axis=0)
        expanded_mask = np.repeat(expanded_mask, ds["S"].values.shape[0], axis=0)

        # Use the mask to select the area of interest
        S_merge.extend(np.where(expanded_mask == 1, ds["S"].values, np.nan))
        T_merge.extend(np.where(expanded_mask == 1, ds["T"].values, np.nan))

        ds.close()

    # Convert the lists to numpy arrays
    S_merge = np.array(S_merge)
    T_merge = np.array(T_merge)

    # Get the dimensions of the data
    time_steps, yc_dim_, xc_dim_ = S_merge.shape
    # Define the number of steps per period
    steps_per_period = 15 * 24
    valid_steps = (time_steps // steps_per_period) * steps_per_period

    # Reshape the data to have the 15-day periods
    S_merge = S_merge[:valid_steps, :, :]
    S_reshaped = S_merge.reshape(-1, steps_per_period, yc_dim_, xc_dim_)

    T_merge = T_merge[:valid_steps, :, :]
    T_reshaped = T_merge.reshape(-1, steps_per_period, yc_dim_, xc_dim_)

    # Calculate the 15-day average and standard deviation and don't consider NaN values
    S_15day_avg = np.nanmean(S_reshaped, axis=1)
    T_15day_avg = np.nanmean(T_reshaped, axis=1)

    S_15day_sd = np.nanstd(S_reshaped, axis=1)
    T_15day_sd = np.nanstd(T_reshaped, axis=1)

    # Calculate the ratio of dry measurements to all measurements over a 15-day periods
    dry_measurement_ratio = np.isnan(S_reshaped).sum(axis=1) / S_reshaped.shape[1] * 100
    # Don't consider the ratio outside of the Wadden Sea area
    expanded_mask = np.expand_dims(mask, axis=0)
    expanded_mask = np.repeat(expanded_mask, S_15day_avg.shape[0], axis=0)
    dry_measurement_ratio = np.where(expanded_mask == 1, dry_measurement_ratio, np.nan)

    # Check if there are NaN values in the Wadden Sea area
    # - which means that there were dry points for the 15-day period
    area_values = S_15day_avg[expanded_mask == 1]
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
    n_periods = time_steps // steps_per_period
    middle_times = [start_date + interval * i + interval / 2 for i in range(n_periods)]

    # Read data from the first data file to get the time units,
    # Use Dataset from netCDF4 library (not xarray - there is a problem in reading time units)
    ds_d = Dataset(data_dir + "/" + files_in_dir1[0], "r")
    # Create the time values for the dataset
    time = date2num(middle_times, units=ds_d["time"].units, calendar="standard")

    # Create a new file with the processed data
    create_file_to_save_processed_data(
        ds_st["xc"],
        ds_st["yc"],
        ds_d["time"].units,
        time,
        S_15day_avg,
        T_15day_avg,
        S_15day_sd,
        T_15day_sd,
        dry_measurement_ratio,
        processed_data_file,
    )


def main():
    """
    This function processes .nc files with data gathered from Wadden Sea
    It takes input parameters for sea area boundaries, a directory of data
    to process, and an output path to the file where to save processed files.
    It uses argparse to parse command-line arguments. Once the arguments are parsed,
    it calls the `process_files` function to handle the core file processing logic.

    """

    # Arguments parsing
    parser = argparse.ArgumentParser(description="Process a .nc files.")
    parser.add_argument(
        "dws_boundaries_area",
        type=str,
        help="Path to the file where the boundaries are stored.",
    )
    parser.add_argument(
        "data_dir",
        type=str,
        help="Path to the directory where the files are stored.",
    )
    parser.add_argument(
        "processed_data_file",
        nargs="?",
        default="15_days_avg_std.nc",
        type=str,
        help="Path to the file where the data will be stored.",
    )
    args = parser.parse_args()

    # Process the file
    process_files(args.dws_boundaries_area, args.data_dir, args.processed_data_file)


if __name__ == "__main__":
    main()

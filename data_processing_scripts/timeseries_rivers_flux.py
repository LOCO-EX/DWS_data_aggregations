#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

PATH_ROOT = Path(__file__).parent


def timeseries_rivers_flux(
    save=False,
    path_file: str | Path = "",
    path_output: str | Path = "",
):
    """Creates a netCDF with the volume flux from rivers into the Dutch Wadden Sea five times an hour.

    Parameters
    ----------
    save : bool, optional
        Whether to save the file, by default False
    path_file : str | Path, optional
        The relative path to the rivers file, by default ""
    path_output : str | Path, optional
        The relative path to the output directory, by default ""
    """

    # Get Rivers netCDF file
    if not Path(PATH_ROOT / path_file).resolve().is_file():
        print(f"{(PATH_ROOT / path_file).resolve()} does not exist")
        quit()
    ds_rivers = xr.open_dataset(PATH_ROOT / path_file)
    da_time = ds_rivers["time"]
    da_time.attrs = {"standard_name": "time", "long_name": "time of measurement"}
    nr_timesteps = len(da_time)

    # Create Station, Lat and Lon List
    station_names = np.array(
        [
            "Cleveringsluizen",
            "Deschans",
            "Dijkmanshuizen",
            "Harlingen",
            "Helsdeur",
            "Krassekeet",
            "Miedema",
            "Oostoever",
            "Ropta",
            "Zandkes",
            "Denoever",
            "Kornwerderzand",
        ]
    )
    nr_stations = len(station_names)
    lats = np.array(
        [
            53.412,
            53.029,
            53.054,
            53.180,
            52.943,
            53.098,
            53.310,
            52.932,
            53.207,
            53.064,
            52.937,
            53.075,
        ]
    )
    lons = np.array(
        [
            6.190,
            4.829,
            4.872,
            5.417,
            4.793,
            4.899,
            5.627,
            4.795,
            5.434,
            4.876,
            5.045,
            5.327,
        ]
    )

    # Create volume_flux array
    np_volume_flux = np.zeros(shape=(nr_timesteps, nr_stations))
    for i, station in enumerate(station_names):
        np_volume_flux[:, i] = ds_rivers[station].values

    # Sort by size
    index_by_size = np.argsort(np.sum(np_volume_flux, 0))[::-1]
    station_names[:] = station_names[index_by_size]
    lats[:] = lats[index_by_size]
    lons[:] = lons[index_by_size]
    np_volume_flux[:, :] = np_volume_flux[:, index_by_size]

    # Create volume_flux meta data
    volume_flux_meta = dict(
        long_name="water volume flux",
        standard_name="water_volume_transport_into_sea_water_from_rivers",
        units="m3 s-1",
        coordinates="lat lon station_name",
    )

    # Create a new Dataset
    ds = xr.Dataset(
        {
            "volume_flux": (("time", "station"), np_volume_flux, volume_flux_meta),
            "lon": (
                ("station"),
                lons,
                {
                    "standard_name": "longitude",
                    "long_name": "station longitude",
                    "units": "degrees_east",
                },
            ),
            "lat": (
                ("station"),
                lats,
                {
                    "standard_name": "latitude",
                    "long_name": "station latitude",
                    "units": "degrees_north",
                },
            ),
            "station_name": (
                ("station"),
                station_names,
                {"long_name": "station name", "cf_role": "timeseries_id"},
            ),
        },
        coords={
            "time": da_time,
        },
        attrs={
            "title": "Volume flux from rivers into the Dutch Wadden Sea five times per hour.",
            "conventions": "CF-1.12",
            "featureType": "timeSeries",
            "institution": "www.tue.nl; www.nioz.nl; www.io-warnemuende.de",
            "email": "m.duran.matute@tue.nl; theo.gerkema@nioz.nl; ulf.graewe@io-warnemuende.de",
            "source": "GETM (www.getm.eu)",
            "comment": "This data is provided as part of the NWO/ENW project: LOCO-EX (OCENW.KLEIN.138). The numerical simulations were done thanks to the North-German Supercomputing Alliance (HLRN). ",
            "history": f"Created {datetime.now().replace(tzinfo=timezone.utc).isoformat(timespec="minutes")} using {Path(__file__).name}",
        },
    )

    # Define the encoding for the dataset
    encoding_format = {
        "zlib": True,
        "complevel": 4,
        "shuffle": True,
    }
    encoding = {
        "volume_flux": encoding_format,
        "lon": encoding_format,
        "lat": encoding_format,
    }

    # Save the dataset to a new file
    if save:
        ds.to_netcdf(
            Path(PATH_ROOT / path_output / "rivers_volume_flux.nc").resolve(),
            format="NETCDF4",
            encoding=encoding,
        )

    print("Executed succesfully")


if __name__ == "__main__":
    timeseries_rivers_flux(
        save=True,
        path_file=Path("../../../Data/GETM_RIVERS/rivers.v02.nc"),
        path_output=Path("../../data"),
    )

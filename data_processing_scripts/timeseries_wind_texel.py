#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from datetime import datetime, timezone

import xarray as xr

PATH_ROOT = Path(__file__).parent


def timeseries_wind_texel(
    path_file: str | Path = "",
    path_output: str | Path = "",
):
    """Adds metadata to a netCDF with the wind texel middle point of the Dutch Wadden Sea.

    Parameters
    ----------
    path_file : str | Path, optional
        The relative path to the wind texel file, by default ""
    path_output : str | Path, optional
        The relative path to the output directory, by default ""
    """

    # Get Rivers netCDF file
    if not Path(PATH_ROOT / path_file).resolve().is_file():
        print(f"{(PATH_ROOT / path_file).resolve()} does not exist")
        quit()
    ds_wind = xr.open_dataset(PATH_ROOT / path_file, decode_times=False)

    # Add attributes to the ds_wind dataset
    ds_wind.attrs["title"] = "Wind texel midpoint"
    ds_wind.attrs["conventions"] = "CF-1.12"
    ds_wind.attrs["institution"] = "www.tue.nl; www.nioz.nl;www.io-warnemuende.de"
    ds_wind.attrs["email"] = (
        "m.duran.matute@tue.nl; theo.gerkema@nioz.nl; ulf.graewe@io-warnemuende.de"
    )
    ds_wind.attrs["source"] = "GETM (www.getm.eu)"
    ds_wind.attrs["comment"] = (
        "This data is provided as part of the NWO/ENW project: LOCO-EX (OCENW.KLEIN.138). The numerical simulations were done thanks to the North-German Supercomputing Alliance (HLRN)."
    )
    ds_wind.attrs["history"] = (
        f"Created {datetime.now().replace(tzinfo=timezone.utc).isoformat(timespec='minutes')}."
    )

    # Define the encoding for the dataset
    encoding_format = {
        "zlib": True,
        "complevel": 4,
        "shuffle": True,
    }
    encoding = {
        "u10": encoding_format,
        "v10": encoding_format,
    }

    # Save the dataset to a new file
    ds_wind.to_netcdf(
        Path(PATH_ROOT / path_output / "DWS.wind_texel_mid_point.nc").resolve(),
        format="NETCDF4",
        encoding=encoding,
    )

    print("Saved succesfully")


if __name__ == "__main__":
    timeseries_wind_texel(
        path_file=Path("wind_texel/wind_texel_mid_point.nc"),
        path_output=Path("wind_texel"),
    )

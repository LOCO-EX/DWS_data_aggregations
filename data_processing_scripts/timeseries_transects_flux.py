#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

PATH_ROOT = Path(__file__).parent


def create_volume_flux_da(transect: str, file: Path):
    ds = xr.open_dataset(file)
    vars = list(ds.keys())
    vars_to_keep = ["hn", "uu", "vv"]
    vars_to_drop = list(set(vars) - set(vars_to_keep))
    ds_core = ds.drop_vars(vars_to_drop)

    xc = ds["xc"].values
    yc = ds["yc"].values
    dx = xc[-1] - xc[0]
    dy = yc[-1] - yc[0]
    normal_vector = (dy, -dx) / np.sqrt(dy**2 + dx**2)  # Vector normal to the transect line

    # Layer height (m) * Point width (m) * Velocity (m/s) = Volume flux m3s-1
    # td Point width WRONG; not all 200
    if "uu" in vars and "vv" in vars:
        ds_core = ds_core.assign(
            volume_flux=ds_core["hn"]
            * 200
            * (normal_vector[0] * ds_core["uu"] + normal_vector[1] * ds_core["vv"])
        )
        ds_sum = ds_core.sum(["level", "nbdyp"])
    elif "uu" in vars:
        ds_core = ds_core.assign(volume_flux=ds_core["hn"] * 200 * (normal_vector[0] * ds_core["uu"]))
        ds_sum = ds_core.sum(["level", "xc", "yc"])
    elif "vv" in vars:
        ds_core = ds_core.assign(volume_flux=ds_core["hn"] * 200 * (normal_vector[1] * ds_core["vv"]))
        ds_sum = ds_core.sum(["level", "xc", "yc"])
    else:
        print("Error with velocity")
        quit()

    da = ds_sum["volume_flux"].assign_coords({"transect": transect}).expand_dims(dim={"transect": 1})
    da_n = xr.DataArray(
        data=normal_vector,
        dims=["unit_vector"],
        coords={"unit_vector": ["x", "y"]},
        attrs={"long_name": "Normal vector"},
        name="normal_vector",
    )
    da_n = da_n.assign_coords({"transect": transect}).expand_dims(dim={"transect": 1})

    return da, da_n


def timeseries_transects_flux(
    save=False,
    path_files: str | Path = "",
    path_output: str | Path = "",
):
    """Creates a netCDF file with the volume flux [m3s-1] trough predetermined transects in the Dutch Wadden Sea

    Parameters
    ----------
    save : bool, optional
        Whether to save the file, by default False
    path_files : str | Path, optional
        The path to the folder with the folders with the transect files, by default ""
    path_output : str | Path, optional
        The path where the file will be saved, by default ""
    """

    transects = ["Eierlandsgat", "MardiepTahlweg", "Marsdiep", "VlieInlet", "Watershed3", "Watershed5"]

    das_volume_flux = []
    das_normal_vector = []
    for transect in transects:
        for file in os.listdir(Path(path_files) / transect):  # All the different time
            da, da_n = create_volume_flux_da(transect, Path(path_files) / transect / file)
            das_volume_flux.append(da)
            das_normal_vector.append(da_n)

    ds_merge = xr.merge(das_volume_flux)
    ds_merge["time"].attrs = {"standard_name": "time", "long_name": "time of measurement"}
    ds_merge_normal = xr.merge(das_normal_vector)

    # Create volume_flux meta data; td change
    volume_flux_meta = dict(
        long_name="water volume flux",
        # standard_name="water_volume_transport_into_sea_water_from_rivers",
        units="m3 s-1",
        # coordinates="lat lon station_name",
    )

    # Create a new Dataset
    ds = xr.Dataset(
        {
            "volume_flux": (("transect", "time"), ds_merge["volume_flux"].values, volume_flux_meta),
            "normal_vector": (("transect", "unit_vector"), ds_merge_normal["normal_vector"].values),
            "transect_name": (
                ("transect"),
                transects,
                {"long_name": "transect", "cf_role": "timeseries_id"},
            ),
        },
        coords={"time": ds_merge["time"], "unit_vector": ["x", "y"], "transect": np.arange(len(transects))},
        attrs={
            "title": "Volume flux from transects into the Dutch Wadden Sea three times per hour.",
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
        "normal_vector": encoding_format,
    }

    # Save the dataset to a new file
    if save:
        ds.to_netcdf(
            Path(PATH_ROOT / path_output / "TR.volume_flux.nc").resolve(),
            format="NETCDF4",
            encoding=encoding,
        )

    print("Executed succesfully")


if __name__ == "__main__":
    timeseries_transects_flux(
        save=True,
        path_files=Path(""),
    )

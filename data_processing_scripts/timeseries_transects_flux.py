#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

PATH_ROOT = Path(__file__).parent


def add_transect(da: xr.DataArray, transect: str):
    return da.assign_coords({"transect": transect}).expand_dims(dim={"transect": 1})


def create_volume_flux_das(transect: str, file: Path):
    ds = xr.open_dataset(file)
    vars = list(ds.keys())
    dims = list(ds.dims.keys())
    vars_to_keep = ["hn", "uu", "vv", "salt"]
    vars_to_drop = list(set(vars) - set(vars_to_keep))
    ds_core = ds.drop_vars(vars_to_drop)

    lat_1 = ds["latc"].values.flatten()[0]
    lon_1 = ds["lonc"].values.flatten()[0]
    lat_2 = ds["latc"].values.flatten()[-1]
    lon_2 = ds["lonc"].values.flatten()[-1]

    xc = ds["xc"].values
    yc = ds["yc"].values

    dx = xc[-1] - xc[0]
    dy = yc[-1] - yc[0]
    normal_vector = (dy, -dx) / np.sqrt(
        dy**2 + dx**2
    )  # Vector normal to the transect line

    if transect == "Eierlandsgat":
        point_width = 209.785
    elif transect == "VlieInlet":
        point_width = 210.950
    elif (
        transect == "Marsdiep"
        or transect == "Borndiep"
        or transect == "Pinkegat"
        or transect == "Watershed1"
        or transect == "Watershed2"
        or transect == "Watershed3"
        or transect == "Watershed4"
        or transect == "Watershed5"
    ):
        point_width = 200
    else:
        print(f"{transect} not a recognized transect")
        quit()

    # Layer height (m) * Point width (m) * Velocity (m/s) = Volume flux m3s-1
    ds_core = ds_core.assign(
        volume_flux=ds_core["hn"]
        * point_width
        * (normal_vector[0] * ds_core["uu"] + normal_vector[1] * ds_core["vv"])
    )

    vars_to_sum = ["level", "nbdyp"] if "nbdyp" in dims else ["level", "xc", "yc"]

    ds_core = ds_core.assign(salinity_flux=lambda x: x.volume_flux * x.salt)
    ds_sum = ds_core.sum(vars_to_sum)

    da_vol_flux = add_transect(ds_sum["volume_flux"], transect)
    da_salt_flux = add_transect(ds_sum["salinity_flux"], transect)
    da_n = xr.DataArray(
        data=normal_vector,
        dims=["unit_vector"],
        attrs={"long_name": "Normal vector"},
        name="normal_vector",
    )
    da_n = add_transect(da_n, transect)
    da_coords_1 = xr.DataArray(
        data=[lat_1, lon_1],
        dims=["coordinates"],
        attrs={"description": "The coordinates of the first point in the transect"},
        name="transect_point_1",
    )
    da_coords_1 = add_transect(da_coords_1, transect)
    da_coords_2 = xr.DataArray(
        data=[lat_2, lon_2],
        dims=["coordinates"],
        attrs={"description": "The coordinates of the last point in the transect"},
        name="transect_point_2",
    )
    da_coords_2 = add_transect(da_coords_2, transect)

    ds.close()

    return da_vol_flux, da_salt_flux, da_n, da_coords_1, da_coords_2


def create_merged_das(transects: list[str], path_files: str | Path):
    das_volume_flux = []
    das_salinity_flux = []
    das_normal_vector = []
    das_coords_1 = []
    das_coords_2 = []
    for transect in transects:
        for file in os.listdir(Path(path_files) / transect):  # All the different time
            da_vol_flux, da_salt_flux, da_n, da_coords_1, da_coords_2 = (
                create_volume_flux_das(transect, Path(path_files) / transect / file)
            )
            das_volume_flux.append(da_vol_flux)
            das_salinity_flux.append(da_salt_flux)
            das_normal_vector.append(da_n)
            das_coords_1.append(da_coords_1)
            das_coords_2.append(da_coords_2)

    ds_merge_vol_flux = xr.merge(das_volume_flux)
    ds_merge_vol_flux["time"].attrs = {
        "standard_name": "time",
        "long_name": "time of measurement",
    }
    da_time = ds_merge_vol_flux["time"]

    np_merge_vol_flux = ds_merge_vol_flux["volume_flux"].values.transpose()
    np_merge_salt_flux = xr.merge(das_salinity_flux)["salinity_flux"].values.transpose()
    np_merge_normal = xr.merge(das_normal_vector)["normal_vector"].values
    np_merge_coord_1 = xr.merge(das_coords_1)["transect_point_1"].values
    np_merge_coord_2 = xr.merge(das_coords_2)["transect_point_2"].values

    return (
        da_time,
        np_merge_vol_flux,
        np_merge_salt_flux,
        np_merge_normal,
        np_merge_coord_1,
        np_merge_coord_2,
    )


def timeseries_transects_flux(
    save=False,
    path_files: str | Path = "",
    path_output: str | Path = "",
):
    """Creates a netCDF file with the volume flux [m3s-1] and the salinity flux [1e-3m3s-1] trough predetermined transects in the Dutch Wadden Sea

    Parameters
    ----------
    save : bool, optional
        Whether to save the file, by default False
    path_files : str | Path, optional
        The path to the folder with the folders with the transect files, by default ""
    path_output : str | Path, optional
        The path where the file will be saved, by default ""
    """

    # transects = ["Eierlandsgat", "MardiepTahlweg", "Marsdiep", "VlieInlet", "Watershed3", "Watershed5", "Borndiep", "C0", "C1", "C2", "C5", "C6", "Pinkegat", "Watershed1", "Watershed2", "Watershed4"]
    transects = [
        "Eierlandsgat",
        "Marsdiep",
        "VlieInlet",
        "Watershed3",
        "Watershed5",
        "Borndiep",
        "Pinkegat",
        "Watershed1",
        "Watershed2",
        "Watershed4",
    ]

    (da_time, np_vol_flux, np_salt_flux, np_normal, np_point_1, np_point_2) = (
        create_merged_das(transects, path_files)
    )

    # Create volume and salinity flux meta data
    volume_flux_meta = dict(
        long_name="water volume flux",
        units="m3 s-1",
        description="Positive into the DWS, for boundary transects",
    )
    salinity_flux_meta = dict(
        long_name="water salinity flux",
        units="1e-3 m3 s-1",
        description="Positive into the DWS, for boundary transects",
    )

    # Create a new Dataset
    ds = xr.Dataset(
        data_vars={
            "volume_flux": (
                ("time", "transect"),
                np_vol_flux,
                volume_flux_meta,
            ),
            "salinity_flux": (
                ("time", "transect"),
                np_salt_flux,
                salinity_flux_meta,
            ),
            "normal_vector": (("transect", "unit_vector"), np_normal),
            "transect_point_1": (("transect", "coordinates"), np_point_1),
            "transect_point_2": (("transect", "coordinates"), np_point_2),
            "transect_name": (
                ("transect"),
                transects,
                {"long_name": "transect", "cf_role": "timeseries_id"},
            ),
        },
        coords={
            "time": da_time,
            "unit_vector": ["x", "y"],
            "coordinates": ["latitude", "longitude"],
            "transect": np.arange(len(transects)),
        },
        attrs={
            "title": "Volume and salinity flux from transects into the Dutch Wadden Sea three times per hour.",
            "conventions": "CF-1.12",
            "featureType": "timeSeries",
            "institution": "www.tue.nl; www.nioz.nl; www.io-warnemuende.de",
            "email": "m.duran.matute@tue.nl; theo.gerkema@nioz.nl; ulf.graewe@io-warnemuende.de",
            "source": "GETM (www.getm.eu)",
            "comment": "This data is provided as part of the NWO/ENW project: LOCO-EX (OCENW.KLEIN.138). The numerical simulations were done thanks to the North-German Supercomputing Alliance (HLRN).",
            "history": f"Created {datetime.now().replace(tzinfo=timezone.utc).isoformat(timespec="minutes")}.",
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
        "salinity_flux": encoding_format,
        "normal_vector": encoding_format,
        "transect_point_1": encoding_format,
        "transect_point_2": encoding_format,
    }

    # Save the dataset to a new file
    if save:
        ds.to_netcdf(
            Path(PATH_ROOT / path_output / "TR.volume_salt_flux.nc").resolve(),
            format="NETCDF4",
            encoding=encoding,
        )

    print("Executed succesfully")


if __name__ == "__main__":
    timeseries_transects_flux(
        save=True,
        path_files=Path(""),
    )

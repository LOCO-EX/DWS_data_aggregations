#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
import xarray as xr
from skimage.morphology import flood_fill

PATH_ROOT = Path(__file__).parent


def create_boundary(path_to_tracer: str | Path = "", path_to_output: str | Path = ""):
    """Create a netCDF DWS mask

    This function downloads the file containing the boundary inddices.
    It will add two points to the mask to ensure the boundary encloses an area.
    The area is then determined using a floodfill algorithm.
    This area is converted to a mask and saved as an netCDF file.

    Parameters
    ----------
    path_to_tracer : str | Path
        The relative path to the file with the tracer, by default ""
    path_to_output : str | Path
        The relative path to the folder where the file will be output, by default ""
    """

    # Download boundary nc from https://data.4tu.nl/datasets/6929d89f-e8cb-463d-b490-3265132841f5
    file_response = requests.get(
        "https://data.4tu.nl/file/6929d89f-e8cb-463d-b490-3265132841f5/3b205b0a-eb9b-4765-9fce-5336538efec3"
    )
    path_temp_boundary_file = PATH_ROOT / path_to_output / "boundary.nc"
    open(path_temp_boundary_file, "wb").write(file_response.content)

    # Open DS boundary file
    ds_boundary = xr.open_dataset(path_temp_boundary_file)
    # Get boundary indices
    bounds = ds_boundary["bdr_dws"].values.copy()
    # Delete boundary file
    ds_boundary.close()
    path_temp_boundary_file.unlink(missing_ok=False)

    # Open DS trace
    ds_trace = xr.open_dataset(PATH_ROOT / path_to_tracer)
    # Get x and y values
    da_xs = ds_trace["xc"]
    da_ys = ds_trace["yc"]
    xs = da_xs.values
    ys = da_ys.values
    # Close DS
    ds_trace.close()

    # Check if xy resolution is homogeneous and isotropic
    if (xs[1] - xs[0]) == (ys[1] - ys[0]) == (xs[-1] - xs[-2]) == (ys[-1] - ys[-2]):
        xy_resolution = xs[1] - xs[0]
    else:
        print("Error")
        quit()

    x_min = xs[0]
    y_min = ys[0]

    # Intitialize mask; 1 at boundary, 0 elsewhere
    mask_boundary = np.zeros((xs.size, ys.size), dtype=int)
    for x, y in bounds:
        mask_boundary[int((x - x_min) // xy_resolution), int((y - y_min) // xy_resolution)] = 1

    # Interpolate 2 points to close the boundary
    mask_boundary[745, 198] = 1
    mask_boundary[184, 116] = 1

    # Fill the inside area with 2
    bounds_as_2d_fill = flood_fill(mask_boundary, (400, 200), 2, connectivity=1)  # 2 as inner value

    # Make mask where it is 2
    mask_outside_area = bounds_as_2d_fill == 2
    # Transpose for (yc, xc) convention
    mask_outside_area = np.transpose(mask_outside_area)

    ds = xr.Dataset(
        data_vars={
            "mask_dws": (
                ["yc", "xc"],
                mask_outside_area,
                {"coordinates": "xc yc", "long_name": "mask of dws"},
            ),
        },
        coords=dict(xc=ds_trace["xc"], yc=ds_trace["yc"]),
        attrs=dict(
            title="Dutch Wadden Sea - 200 m resolution: area mask",
            Conventions="CF-1.12",
            institution="www.tue.nl; www.nioz.nl; www.io-warnemuende.de",
            email="m.duran.matute@tue.nl; theo.gerkema@nioz.nl; ulf.graewe@io-warnemuende.de",
            source="GETM (www.getm.eu)",
            comment="This data is provided as part of the NWO/ENW project: LOCO-EX (OCENW.KLEIN.138). The numerical simulations were done thanks to the North-German Supercomputing Alliance (HLRN).",
            history=f"Created {datetime.now().replace(tzinfo=timezone.utc).isoformat(timespec="minutes")} using create_boundary.py",
        ),
    )

    ds.to_netcdf(PATH_ROOT / path_to_output / "dws_boundary_area_200x200m.nc", "w", format="NETCDF4")

    print("File created")


if __name__ == "__main__":
    create_boundary()

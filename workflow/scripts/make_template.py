"""Make the analysis grid: one raster that every covariate is reprojected onto.

Cells are res_m wide, with edges on multiples of res_m. Cells inside the
boundary are 1, cells outside are missing, so the template also works as a
land mask.
"""

import geopandas as gpd
import numpy as np
import rioxarray  # noqa: F401  (adds the .rio methods used below)
import xarray as xr


def run(boundary_path: str, out_path: str, res_m: float) -> None:
    boundary = gpd.read_file(boundary_path)

    xmin, ymin, xmax, ymax = boundary.total_bounds
    xmin, ymin = np.floor(xmin / res_m) * res_m, np.floor(ymin / res_m) * res_m
    xmax, ymax = np.ceil(xmax / res_m) * res_m, np.ceil(ymax / res_m) * res_m
    xs = np.arange(xmin + res_m / 2, xmax, res_m)  # cell centres
    ys = np.arange(ymax - res_m / 2, ymin, -res_m)  # top to bottom

    grid = xr.DataArray(
        np.ones((len(ys), len(xs)), dtype="float32"),
        coords={"y": ys, "x": xs},
        dims=("y", "x"),
    )
    grid = grid.rio.write_crs(boundary.crs).rio.write_nodata(np.nan)
    grid = grid.rio.clip(boundary.geometry, boundary.crs, all_touched = True, drop=False)
    grid.rio.to_raster(out_path)


if __name__ == "__main__":
    run(
        boundary_path=snakemake.input[0],
        out_path=snakemake.output[0],
        res_m=snakemake.config["general"]["cell_res_km"] * 1000,
    )
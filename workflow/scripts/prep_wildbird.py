"""Combine species' eBird relative abundance layers into one waterbird layer.

For each species:
  1. crop the global weekly raster to (a little beyond) Australia,
  2. reproject onto the common analysis grid (bilinear),
  3. optionally average weeks into months.
The layers are then summed across species, and cells outside the
boundary are set to missing.
"""

from pathlib import Path

import pandas as pd
import rasterio
import rioxarray
import xarray as xr
from rasterio.enums import Resampling
from rasterio.warp import transform_bounds


def load_species(path: Path, template: xr.DataArray) -> xr.DataArray:
    """Read one species' weekly raster and put it on the analysis grid."""
    da = rioxarray.open_rasterio(path, masked=True)  # missing cells become NaN

    # each band is one week; its description holds the week's date
    with rasterio.open(path) as src:
        dates = pd.to_datetime(list(src.descriptions))
    da = da.assign_coords(band=dates)

    # crop to the grid's extent plus a 100 km margin *before* reprojecting: much faster
    # than reprojecting the whole globe, and the margin avoids edge effects
    # from the bilinear interpolation
    xmin, ymin, xmax, ymax = template.rio.bounds()
    margin = 100_000  # 100 km
    bounds = transform_bounds(
        template.rio.crs, # source crs
        da.rio.crs, # destination crs
        xmin - margin, ymin - margin, xmax + margin, ymax + margin,
        densify_pts=21,
    )
    da = da.rio.clip_box(*bounds)

    return da.rio.reproject_match(template, resampling=Resampling.bilinear)


def run(tif_dir: str, template_path: str, out_path: str, temporal_res: str) -> None:
    template = rioxarray.open_rasterio(template_path).squeeze("band", drop=True) # turns the template's single band into a plain 2-D grid, so that .where() lines up with the y and x dimensions.

    layers = []
    for path in sorted(Path(tif_dir).glob("*.tif")):
        da = load_species(path, template)

        if temporal_res == "month":
            da = da.groupby("band.month").mean()  # average the weeks in each month
            da = da.rename(month="band") # rename layers to 'band'

        # da = da / da.max()  # scale to 0-1 by the species' maximum (NaN ignored)
        layers.append(da)
        print(f"Processed {path.name}")

    # sum across species; cells where every species is missing become 0,
    combined = xr.concat(layers, dim="species").sum("species")

    # set cells outside the boundary to missing, using the template as a mask
    combined = combined.where(template.notnull().values)
    combined = combined.rio.write_crs(template.rio.crs)

    if temporal_res == "month":
        names = [pd.Timestamp(2000, m, 1).strftime("%b") for m in combined["band"].values]
    else:
        names = [pd.Timestamp(d).strftime("%Y-%m-%d") for d in combined["band"].values]
    
    out = combined.astype("float32")
    out.attrs["long_name"] = tuple(names)  # one name per band, written as band descriptions
    out.rio.to_raster(out_path)

if __name__ == "__main__":
    run(
        tif_dir=snakemake.input.tifs,
        template_path=snakemake.input.template,
        out_path=snakemake.output[0],
        temporal_res=snakemake.params.temporal_res
    )
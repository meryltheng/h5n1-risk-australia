"""Put GLW4 chicken density (animals per km²) onto the analysis grid."""

import rioxarray
from rasterio.enums import Resampling
from rasterio.warp import transform_bounds


def run(tif_path: str, template_path: str, out_path: str) -> None:
    template = rioxarray.open_rasterio(template_path).squeeze("band", drop=True)
    da = rioxarray.open_rasterio(tif_path, masked=True)  # missing cells become NaN

    # crop to the grid's extent plus a 100 km margin before reprojecting
    xmin, ymin, xmax, ymax = template.rio.bounds()
    margin = 100_000  # metres
    bounds = transform_bounds(
        template.rio.crs, da.rio.crs,
        xmin - margin, ymin - margin, xmax + margin, ymax + margin,
        densify_pts=21,
    )
    da = da.rio.clip_box(*bounds)

    # each 25 km cell gets the mean density of the ~10 km source cells it covers
    da = da.rio.reproject_match(template, resampling=Resampling.average)

    # set cells outside the boundary to missing, using the template as a mask
    da = da.where(template.notnull().values)

    da.astype("float32").rio.to_raster(out_path)


if __name__ == "__main__":
    run(
        tif_path=snakemake.input.tif,
        template_path=snakemake.input.template,
        out_path=snakemake.output[0],
    )
"""Make the Australia boundary used to crop and mask every layer.

Clips the ABS Australia outline to an extent that excludes the external
territories, merges it into one shape and reprojects it to the analysis CRS.
"""

import geopandas as gpd
from shapely.geometry import box


def run(in_path: str, out_path: str, extent: list[float], crs: str) -> None:
    aus = gpd.read_file(in_path)
    aus = aus[aus.geometry.notna()]  # ABS includes some records with no shape
    aus = aus.clip(box(*extent))  # extent is [xmin, ymin, xmax, ymax] in lon/lat
    aus = aus.dissolve().to_crs(crs)
    aus[["geometry"]].to_file(out_path)


if __name__ == "__main__":
    run(
        in_path=snakemake.input[0],
        out_path=snakemake.output[0],
        extent=snakemake.config["general"]["extent_lonlat"],
        crs=snakemake.config["general"]["proj_crs"]
    )
"""Give each DAFF record a location: a random point inside its ABS locality.

DAFF reports locations only by locality name and state. Each record is
matched to an ABS Suburbs and Localities (SAL) polygon by name and state,
then given one random point inside that polygon. Records that cannot be
matched keep empty coordinates and are listed in the printed output.
"""

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

# DAFF uses abbreviations; the ABS file uses full state names
STATE_NAMES = {
    "NSW": "New South Wales",
    "VIC": "Victoria",
    "QLD": "Queensland",
    "SA": "South Australia",
    "WA": "Western Australia",
    "TAS": "Tasmania",
    "NT": "Northern Territory",
    "ACT": "Australian Capital Territory",
}


def match_localities(daff: pd.DataFrame, sal: gpd.GeoDataFrame) -> pd.Series:
    """For each DAFF record, the row label of its SAL polygon (missing if unmatched).

    ABS adds a state suffix to names that occur in more than one state, e.g.
    "Richmond (Tas.)". The suffix is removed, and records are matched on name
    (ignoring upper/lower case) and state together.
    """
    keys = pd.DataFrame({
        "name": sal["SAL_NAME21"].str.strip().str.replace(r"\s*\(.*\)$", "", regex=True).str.lower(),
        "state": sal["STE_NAME21"],
        "sal_row": sal.index,
    })
    # a name + state pair should be unique; drop any that are not, so they stay unmatched
    keys = keys.drop_duplicates(["name", "state"], keep=False)

    records = pd.DataFrame({
        "name": daff["locality"].str.strip().str.lower(),
        "state": daff["jurisdiction"].str.strip().str.upper().map(STATE_NAMES),
    })
    matched = records.merge(keys, on=["name", "state"], how="left")
    return pd.Series(matched["sal_row"].values, index=daff.index)


def random_point_in(polygon, rng: np.random.Generator) -> tuple[float, float]:
    """A random point inside a polygon: draw points in its bounding box until one falls inside."""
    minx, miny, maxx, maxy = polygon.bounds
    while True:
        x, y = rng.uniform(minx, maxx), rng.uniform(miny, maxy)
        if polygon.contains(Point(x, y)):
            return x, y


def run(daff_path: str, sal_path: str, out_path: str, seed: int) -> None:
    daff = pd.read_csv(daff_path)
    sal = gpd.read_file(sal_path)
    sal = sal[sal.geometry.notna() & ~sal.geometry.is_empty]  # some ABS records have no shape

    sal_row = match_localities(daff, sal)

    unmatched = daff.loc[sal_row.isna(), ["locality", "jurisdiction"]].drop_duplicates()
    print(f"{sal_row.isna().sum()} of {len(daff)} DAFF records not matched to a locality")
    if len(unmatched):
        print(unmatched.to_string(index=False))

    rng = np.random.default_rng(seed)
    daff["sal_code"] = pd.NA
    daff["longitude"] = np.nan
    daff["latitude"] = np.nan
    for i, row in sal_row.dropna().items():
        polygon = sal.loc[row]
        daff.loc[i, "sal_code"] = polygon["SAL_CODE21"]
        daff.loc[i, ["longitude", "latitude"]] = random_point_in(polygon.geometry, rng)

    daff.to_csv(out_path, index=False)


if __name__ == "__main__":
    run(
        daff_path=snakemake.input.daff,
        sal_path=snakemake.input.sal,
        out_path=snakemake.output[0],
        seed=snakemake.config["general"]["seed"],
    )
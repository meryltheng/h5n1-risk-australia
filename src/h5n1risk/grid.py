"""Link point records to the analysis grid (the template raster)."""

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform as transform_coords
from rasterio.transform import rowcol, xy


def cells_from_template(template_path: str) -> pd.DataFrame:
    """One row per grid cell inside the boundary: cell_id, row, col, and centre x, y in km."""
    with rasterio.open(template_path) as src:
        data = src.read(1)
        transform = src.transform
    rows, cols = np.nonzero(np.isfinite(data))         # cells outside the boundary are NaN
    xs, ys = xy(transform, rows, cols)                 # cell centres, metres
    return pd.DataFrame({
        "cell_id": np.arange(len(rows)),
        "row": rows,
        "col": cols,
        "x": np.asarray(xs) / 1000,
        "y": np.asarray(ys) / 1000,
    })


def assign_cells(lon, lat, template_path: str, cells: pd.DataFrame) -> np.ndarray:
    """cell_id for each longitude/latitude point; -1 if it falls outside the grid."""
    with rasterio.open(template_path) as src:
        transform, crs, shape = src.transform, src.crs, src.shape

    px, py = transform_coords("EPSG:4326", crs, list(lon), list(lat))
    r, c = (np.asarray(v) for v in rowcol(transform, px, py))

    lookup = np.full(shape, -1)                        # (row, col) -> cell_id
    lookup[cells["row"], cells["col"]] = cells["cell_id"]

    ids = np.full(len(r), -1)
    inside = (r >= 0) & (r < shape[0]) & (c >= 0) & (c < shape[1])
    ids[inside] = lookup[r[inside], c[inside]]
    return ids
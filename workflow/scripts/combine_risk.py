"""Combine the kernel layer with the waterbird and poultry layers into relative risk maps.

All inputs are already on the template grid, so cells line up directly: no
extraction to a separate grid and no joining by row order.

    wildlife risk = rescale( kernel x waterbird )
    poultry risk  = additive:        rescale( kernel x (waterbird + poultry) )   [as in R]
                    multiplicative:  wildlife risk x poultry exposure

The results rank cells; they are not probabilities of infection.
"""

import matplotlib.pyplot as plt
import numpy as np
import rasterio

from h5n1risk import maps


def rescale01(a: np.ndarray) -> np.ndarray:
    """Rescale to 0-1 by the minimum and maximum, ignoring missing cells."""
    lo, hi = np.nanmin(a), np.nanmax(a)
    return (a - lo) / (hi - lo)


def compute_risk(kernel, waterbird, poultry, spillover="additive", half_sat=None):
    """Relative risk surfaces from the three layers (2-D arrays on the same grid).

    spillover : "additive" reproduces the R version; "multiplicative" requires
        both virus arrival and poultry, so high poultry density cannot make up
        for low wildlife risk.
    half_sat : for "multiplicative", the poultry density (birds per km²) at
        which exposure reaches half its maximum. None = median of cells with
        poultry.
    """
    z_wb = rescale01(waterbird)
    wildlife = rescale01(kernel * z_wb)

    if spillover == "additive":
        z_pou = rescale01(poultry)
        poultry_risk = rescale01(kernel * (z_wb + z_pou))
    elif spillover == "multiplicative":
        if half_sat is None:
            # density above which half of all birds live (cells are equal-area,
            # so weighting each cell by its density weights it by its birds)
            positive = np.sort(poultry[poultry > 0])
            half_sat = positive[np.searchsorted(np.cumsum(positive) / positive.sum(), 0.5)]
        print(f"Poultry half-saturation density: {half_sat:.0f} birds per km²")
        z_pou = poultry / (poultry + half_sat)          # 0 with no poultry, towards 1 when dense
        poultry_risk = rescale01(wildlife * z_pou)
    else:
        raise ValueError(f"spillover must be 'additive' or 'multiplicative', not {spillover!r}")

    return {"waterbird_scaled": z_wb, "poultry_scaled": z_pou,
            "rel_risk_wildlife": wildlife, "rel_risk_poultry": poultry_risk}


def read_band(path, band=1):
    with rasterio.open(path) as src:
        return src.read(band).astype(float), src.profile, src.descriptions


def run(kernel_path, waterbird_path, poultry_path, out_path, map_path,
        risk_month, spillover, half_sat):
    kernel, profile, _ = read_band(kernel_path, 1)
    detected, _, _ = read_band(kernel_path, 2)
    waterbird, wb_profile, wb_names = read_band(waterbird_path, risk_month)
    poultry, pou_profile, _ = read_band(poultry_path, 1)
    poultry[np.isnan(poultry) & ~np.isnan(detected)] = 0.0 # fill nans within aus with 0

    # all layers must be on the same grid, or cells would not line up
    for name, p in [("waterbird", wb_profile), ("poultry", pou_profile)]:
        if (p["transform"], p["width"], p["height"]) != (profile["transform"], profile["width"], profile["height"]):
            raise ValueError(f"The {name} layer is not on the kernel layer's grid.")
    print(f"Waterbird layer: band {risk_month} ({wb_names[risk_month - 1]})")

    layers = compute_risk(kernel, waterbird, poultry, spillover, half_sat)
    layers["kernel"] = kernel
    layers["detected"] = detected

    profile.update(count=len(layers), dtype="float32", nodata=np.nan)
    with rasterio.open(out_path, "w", **profile) as dst:
        for i, (name, values) in enumerate(layers.items(), start=1):
            dst.write(values.astype("float32"), i)
            dst.set_band_description(i, name)

    # maps: wildlife and poultry risk side by side
    with rasterio.open(out_path) as src:
        b = src.bounds
    extent = [b.left / 1000, b.right / 1000, b.bottom / 1000, b.top / 1000]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    maps.draw_risk_map(axes[0], layers["rel_risk_wildlife"], detected, extent, "H5N1 in wild birds", show_detected=False)
    im = maps.draw_risk_map(axes[1], layers["rel_risk_poultry"], detected, extent,
                            f"Spillover to poultry ({spillover})", show_detected=False)
    fig.colorbar(im, ax=axes, label="Relative risk", shrink=0.7)
    fig.savefig(map_path, dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    run(
        kernel_path=snakemake.input.kernel,
        waterbird_path=snakemake.input.waterbird,
        poultry_path=snakemake.input.poultry,
        out_path=snakemake.output.risk,
        map_path=snakemake.output.map,
        risk_month=snakemake.params.risk["risk_month"],
        spillover=snakemake.params.risk["spillover"],
        half_sat=snakemake.params.risk["poultry_half_sat"],
    )
"""Map the kernel-only relative risk from one fitted kernel."""

import json

import matplotlib.pyplot as plt
import pandas as pd
import rasterio

from h5n1risk import grid, kernels, maps


def run(fits_path, table_path, risk_set_path, template_path, layer_path, map_path, kernel=None):
    with open(fits_path) as f:
        fits = {k: kernels.KernelFit.from_dict(v) for k, v in json.load(f).items()}
    name = kernel or pd.read_csv(table_path)["kernel"].iloc[0]  # best by AIC if not set
    print(f"Mapping kernel layer from: {name}")

    risk_set = pd.read_csv(risk_set_path, parse_dates=["detect_date"])
    activated = risk_set["detect_date"].notna()
    z, is_detected = kernels.predict_kernel_layer(
        risk_set[["x", "y"]].to_numpy(),
        risk_set.loc[activated, ["x", "y"]].to_numpy(),
        fits[name],
    )
    z[is_detected] = 1.0   # detected cells: maximum kernel value instead of masked
    
    grid.cells_to_raster(
        {"relative_risk": z, "detected": is_detected.astype(float)},
        risk_set, template_path, layer_path,
    )

    with rasterio.open(layer_path) as src:
        risk, detected = src.read(1), src.read(2)
        b = src.bounds
    extent = [b.left / 1000, b.right / 1000, b.bottom / 1000, b.top / 1000]  # km

    fig, ax = plt.subplots(figsize=(8, 7))
    im = maps.draw_risk_map(ax, risk, detected, extent,
                            f"Kernel-only relative risk ({kernels.LABELS[name]})")
    fig.colorbar(im, ax=ax, label="Relative risk (kernel only)", shrink=0.7)
    fig.savefig(map_path, dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    run(
        fits_path=snakemake.input.fits,
        table_path=snakemake.input.table,
        risk_set_path=snakemake.input.risk_set,
        template_path=snakemake.input.template,
        layer_path=snakemake.output.layer,
        map_path=snakemake.output.map,
        kernel=snakemake.params.kernel,
    )
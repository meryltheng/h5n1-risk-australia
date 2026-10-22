"""Fit distance-transmission kernels to the combined outbreak records."""

import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from h5n1risk import grid, kernels


def run(outbreaks_path, template_path, risk_set_path, fits_path, table_path, plot_path, cfg, 
        drop_duplicates=False):
    outbreaks = pd.read_csv(outbreaks_path, parse_dates=["date"])
    outbreaks = outbreaks.dropna(subset=["date", "longitude", "latitude"])

    if drop_duplicates:
        # one record per locality + date + species; prefer WAHIS (recorded coordinates)
        n_before = len(outbreaks)
        outbreaks = (
            outbreaks.assign(_not_wahis=outbreaks["source"] != "WAHIS")
            .sort_values("_not_wahis", kind="stable")         # WAHIS rows first
            .drop_duplicates(["loc", "date", "species"], keep="first")
            .drop(columns="_not_wahis")
        )
        print(f"Dropped {n_before - len(outbreaks)} duplicate records")

    cells = grid.cells_from_template(template_path)
    outbreaks["cell_id"] = grid.assign_cells(
        outbreaks["longitude"], outbreaks["latitude"], template_path, cells
    )
    print(f"{(outbreaks['cell_id'] < 0).sum()} of {len(outbreaks)} records fall outside the grid and are dropped")

    risk_set = kernels.build_risk_set(cells, outbreaks[outbreaks["cell_id"] >= 0])
    kd = kernels.build_kernel_data(
        risk_set,
        step_days=cfg["step_days"],
        inf_offset=cfg["inf_offset_days"],
        inf_period=cfg["inf_period_days"],
        bin_km=cfg["bin_km"],
        max_dist_km=cfg["max_dist_km"] if cfg["max_dist_km"] is not None else np.inf,
    )

    fits = kernels.fit_kernels(kd, kernels=cfg["kernels"], mu=cfg["mu_km"])
    table = kernels.compare_kernels(fits)
    print(table.to_string(index=False))

    risk_set.to_csv(risk_set_path, index=False)
    table.to_csv(table_path, index=False)
    with open(fits_path, "w") as f:
        json.dump({name: fit.to_dict() for name, fit in fits.items()}, f, indent=2)

    # best = fits[table["kernel"].iloc[0]]
    # fig = kernels.plot_kernel(best, xlim=cfg["plot_xlim_km"], seed=cfg["seed"])
    # fig.savefig(plot_path, dpi=150)
    # one panel per kernel, best first, from the shortest distance observed at activation
    x = np.linspace(kd.exact_dist.min(), cfg["plot_xlim_km"], 400)
    quantiles = {
        name: np.nanquantile(kernels.kernel_curves(fits[name], x, seed=cfg["seed"]),
                             [0.025, 0.25, 0.5, 0.75, 0.975], axis=0)
        for name in table["kernel"]
    }
    # y-axis set by the median curves, so a very wide band can't flatten them;
    # bands taller than this are cut off at the top of the panel
    y_top = 1.5 * max(np.nanmax(q[2]) for q in quantiles.values())

    n_rows = int(np.ceil(len(table) / 2))
    fig, axes = plt.subplots(n_rows, 2, figsize=(10, 4 * n_rows), sharex=True, sharey=True)
    axes = axes.ravel()  # grid of panels -> flat list, filled row by row
    for ax, (_, row) in zip(axes, table.iterrows()):
        q = quantiles[row["kernel"]]
        ax.fill_between(x, q[0], q[4], color="salmon", alpha=0.4, label="95%")
        ax.fill_between(x, q[1], q[3], color="salmon", alpha=0.7, label="50%")
        ax.plot(x, q[2], color="firebrick", label="median")
        ax.set(title=f"{row['label']} (ΔAIC {row['dAIC']:.1f})")
        if np.nanmax(q[4]) > y_top:
            ax.text(0.98, 0.95, "95% band extends off the plot", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8, style="italic")
    ax.set_ylim(0, y_top)  # applies to every panel, as the y-axis is shared
    for ax in axes[len(table):]:  # hide spare panels if there are fewer than 4 kernels
        ax.set_visible(False)
    for ax in axes[-2:]:
        ax.set_xlabel("Distance (km)")
    for ax in axes[::2]:
        ax.set_ylabel("Transmission rate (per cell per step)")
    axes[0].legend(loc="upper center")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)


if __name__ == "__main__":
    run(
        outbreaks_path=snakemake.input.outbreaks,
        template_path=snakemake.input.template,
        risk_set_path=snakemake.output.risk_set,
        fits_path=snakemake.output.fits,
        table_path=snakemake.output.table,
        plot_path=snakemake.output.plot,
        cfg={**snakemake.config["kernel"], "seed": snakemake.config["general"]["seed"]},
    )

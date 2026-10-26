"""Shared map drawing for gridded risk layers."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

# stepped colour scale (as in the R maps): most cells are near 0, so even
# steps would show only the cells immediately around detections
BREAKS = [0, 0.01, 0.05, 0.1, 0.25, 0.4, 0.6, 0.8, 1]


def draw_risk_map(ax, risk: np.ndarray, detected: np.ndarray, extent, title: str,
                  show_detected: bool = True):
    """Draw a 0-1 risk grid on ax, optionally with detected cells in blue. Returns the image (for a colour bar).

    extent : [left, right, bottom, top] of the grid, in the units to label the axes with
    show_detected : cover detected cells in blue; if False they keep their risk colour
    """
    cmap = plt.get_cmap("magma", len(BREAKS) - 1)
    im = ax.imshow(risk, extent=extent, cmap=cmap, norm=BoundaryNorm(BREAKS, cmap.N))
    if show_detected:
        # detected cells in a colour outside the risk scale (magma's low end is black)
        ax.imshow(np.where(detected == 1, 1.0, np.nan), extent=extent,
                  cmap=ListedColormap(["deepskyblue"]))
        ax.legend(handles=[Patch(color="deepskyblue", label="Detected")], loc="lower left")
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.set_axis_off()
    return im
"""Distance-transmission kernels fitted to presence-only outbreak records.

Python port of fit_kernels_wahis.R. Same kernel shapes, likelihood and
optional long-range component; the grid itself is built elsewhere (grid.py)
so that the kernel fit and the covariates share one set of cells.

Model, for a cell i at time step t:

    lambda_i(t) = h_i^a * sum over infectious cells j of
                  [ beta0 * K(d_ij) + beta2 * G(d_ij; mu) ]

    P(cell i activates at step t) = 1 - exp(-lambda_i(t))

Activated cells contribute log(1 - exp(-lambda)) at their activation step;
every cell contributes -lambda for each earlier step it stayed at risk.
"""

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import minimize

# ---------------------------------------------------------------------------
# 1. Kernel definitions
# ---------------------------------------------------------------------------
# Each gives a relative transmission rate at distance d (km). beta0 scales
# them, so K only carries the shape.

KERNELS = {
    "ftail": lambda d, rho0, psi: 1 / (1 + (d / rho0) ** psi),
    "pow1": lambda d, rho0, psi: 1 / (1 + d**psi),
    "pow2": lambda d, rho0, psi: (1 + d / rho0) ** (-psi),
    "cauchy": lambda d, rho0, psi: psi / (psi**2 + d**2),
    "exp1": lambda d, rho0, psi: np.exp(-psi * d),
    "exp2": lambda d, rho0, psi: 1 - np.exp(-((d / rho0) ** (-psi))),
    "exp3": lambda d, rho0, psi: np.exp(-((d / rho0) ** psi)),
}

# pow1, cauchy and exp1 have no rho0, so it must not be estimated for them
USES_RHO0 = {"ftail": True, "pow1": False, "pow2": True, "cauchy": False,
             "exp1": False, "exp2": True, "exp3": True}

# two-parameter kernels: the ones with a length scale of their own
KERNELS_2P = ["ftail", "pow2", "exp2", "exp3"]

LABELS = {"ftail": "Fat tailed", "pow1": "Power law type 1",
          "pow2": "Power law type 2", "cauchy": "Cauchy",
          "exp1": "Exponential type 1", "exp2": "Exponential type 2",
          "exp3": "Exponential type 3"}


def g_long(d, mu):
    """Long-range component: fixed exponential shape, mu (km) fixed, not estimated."""
    return np.exp(-d / mu)


PARAM_NAMES = ["log_beta0", "log_rho0", "log_psi", "a", "log_beta2"]
DEFAULT_START = {"log_beta0": np.log(0.01), "log_rho0": np.log(20),
                 "log_psi": np.log(2), "a": 0.0, "log_beta2": np.log(1e-4)}


# ---------------------------------------------------------------------------
# 2. Risk set
# ---------------------------------------------------------------------------

def build_risk_set(cells: pd.DataFrame, detections: pd.DataFrame) -> pd.DataFrame:
    """Give each grid cell the date of its earliest detection (missing if none).

    cells      : one row per cell, with cell_id (0, 1, 2, ...) and x, y in km
    detections : one row per record, with cell_id and date

    The extent of the grid directly affects the fitted kernel: it sets how
    many never-activated cells enter the survival term, and a larger area
    pushes the estimate towards steeper decay. Choose it on ecological
    grounds and report the choice.
    """
    first = detections.groupby("cell_id")["date"].min()
    risk_set = cells.copy()
    risk_set["detect_date"] = risk_set["cell_id"].map(first)
    n_act = risk_set["detect_date"].notna().sum()
    print(f"Risk set: {len(risk_set)} cells, {n_act} activated ({100 * n_act / len(risk_set):.1f}%)")
    return risk_set


# ---------------------------------------------------------------------------
# 3. Likelihood data
# ---------------------------------------------------------------------------

@dataclass
class KernelData:
    """Everything the likelihood needs, precomputed once."""
    act_mat: sparse.csr_matrix   # [cells x distance bins] counts of at-risk step-pairs
    bins: np.ndarray             # midpoint distance (km) of each bin
    exact_cell: np.ndarray       # target cell of each infectious pair at activation
    exact_dist: np.ndarray       # its exact distance (km)
    seeds: np.ndarray            # cells activated with no infectious source
    n_cells: int
    n_steps: int


def build_kernel_data(
    risk_set: pd.DataFrame,
    step_days: int = 1,
    inf_offset: float = -7,
    inf_period: float = 14,
    bin_km: float = 1.0,
    max_dist_km: float = np.inf,
) -> KernelData:
    """Turn activation dates into the pair counts the likelihood uses.

    Each cell's detection date sets its ACTIVATION step (when it stops being
    a target). Its window as a SOURCE is set separately:

        infectious from  activation + inf_offset
        infectious until that + inf_period - 1

    inf_offset is normally negative: birds shed before they die and carcasses
    are found later, so a cell can be a source before it is detected. Neither
    value has data behind it for wild birds, so sweep both.
    """
    x = risk_set["x"].to_numpy(float)
    y = risk_set["y"].to_numpy(float)
    dates = risk_set["detect_date"]
    n = len(risk_set)

    d0, d1 = dates.min(), dates.max()
    act = np.floor((dates - d0).dt.days.to_numpy(float) / step_days) + 1  # NaN if never activated
    n_steps = int((d1 - d0).days // step_days + 1)

    off_steps = int(round(inf_offset / step_days))
    inf_steps = max(1, int(round(inf_period / step_days)))
    if inf_offset != 0 and off_steps == 0:
        warnings.warn(f"inf_offset ({inf_offset} d) rounds to 0 steps at step_days = {step_days}")
    if inf_period < step_days:
        warnings.warn(f"inf_period ({inf_period} d) is shorter than step_days; rounded up to one step")

    s = act + off_steps                                # first infectious step
    e = s + inf_steps - 1                              # last infectious step
    risk_end = np.where(np.isnan(act), n_steps, act - 1)  # last step contributing survival

    src = np.flatnonzero(~np.isnan(s))                 # cells that ever transmit
    if len(src) == 0:
        raise ValueError("No activated cells: nothing to fit.")

    # distances (km) from every cell to every ever-infectious cell
    dist = np.hypot(x[:, None] - x[src], y[:, None] - y[src])

    # survival pair counts: steps that source j is infectious while cell i is
    # still at risk. Steps before the study period (s < 1) are dropped.
    lo = np.maximum(s[src], 1)[None, :]
    hi = np.minimum(e[src][None, :], risk_end[:, None])
    count = np.maximum(hi - lo + 1, 0)
    count[src, np.arange(len(src))] = 0               # no self-transmission

    keep = (count > 0) & (dist <= max_dist_km)
    if not keep.any():
        raise ValueError("No infectious/at-risk pairs: check dates and max_dist_km.")
    rows, cols = np.nonzero(keep)
    bin_idx = np.floor(dist[rows, cols] / bin_km).astype(int)
    n_bins = bin_idx.max() + 1
    bins = (np.arange(n_bins) + 0.5) * bin_km
    # repeated (row, bin) entries are summed
    act_mat = sparse.csr_matrix((count[rows, cols], (rows, bin_idx)), shape=(n, n_bins))

    # exact distances to every source infectious at each cell's activation step
    exact_cell, exact_dist, seeds = [], [], []
    for i in np.flatnonzero(~np.isnan(act)):
        t = act[i]
        j = np.flatnonzero((s[src] <= t) & (t <= e[src]) & (src != i) & (dist[i] <= max_dist_km))
        if len(j) == 0:
            seeds.append(i)
            continue
        exact_cell.append(np.full(len(j), i))
        exact_dist.append(dist[i, j])

    if seeds:
        print(f"{len(seeds)} index case(s) with no infectious source: conditioned on, not fitted.")
    kd = KernelData(
        act_mat=act_mat,
        bins=bins,
        exact_cell=np.concatenate(exact_cell) if exact_cell else np.array([], dtype=int),
        exact_dist=np.concatenate(exact_dist) if exact_dist else np.array([]),
        seeds=np.array(seeds, dtype=int),
        n_cells=n,
        n_steps=n_steps,
    )
    print(f"Likelihood data: {n} cells x {n_bins} distance bins, {len(kd.exact_dist)} exact pairs.")
    return kd


# ---------------------------------------------------------------------------
# 4. Likelihood
# ---------------------------------------------------------------------------

def make_nll(kd: KernelData, kernel: str, hdens=None, estimate_a: bool = True, mu=None):
    """Negative log-likelihood for one kernel, as a function of a dict of parameters.

    Rates are a SUM of two components, each with its own amplitude (beta0,
    beta2); mu = None gives the single-kernel model (beta2 = 0).
    hdens is an optional per-cell multiplier (e.g. human density), scaled to
    mean 1 so beta0 stays the rate for an average cell.
    """
    K = KERNELS[kernel]
    uses_rho0 = USES_RHO0[kernel]
    use_long = mu is not None

    if hdens is None:
        h = np.ones(kd.n_cells)
        estimate_a = False
    else:
        h = np.asarray(hdens, float) / np.nanmean(hdens)
        floor = np.nanmin(h[h > 0]) / 2
        h[~np.isfinite(h) | (h <= 0)] = floor          # avoid 0 ** a

    ex_cells, ex_group = np.unique(kd.exact_cell, return_inverse=True)
    g_bins = g_long(kd.bins, mu) if use_long else 0.0   # mu is fixed: precompute once
    g_exact = g_long(kd.exact_dist, mu) if use_long else 0.0

    def nll(p: dict) -> float:
        beta0 = np.exp(p["log_beta0"])
        beta2 = np.exp(p["log_beta2"]) if use_long else 0.0
        rho0 = np.exp(p["log_rho0"]) if uses_rho0 else np.nan
        psi = np.exp(p["log_psi"])
        hh = h ** p["a"] if estimate_a else h

        with np.errstate(all="ignore"):
            k_bins = K(kd.bins, rho0, psi)
            if not np.all(np.isfinite(k_bins)) or np.any(k_bins < 0):
                return 1e10
            rate_bins = beta0 * k_bins + beta2 * g_bins

            # survival: sum of lambda over every at-risk cell-step
            log_p_surv = np.sum(hh * (kd.act_mat @ rate_bins))

            # activation: lambda at each cell's activation step
            k_exact = K(kd.exact_dist, rho0, psi)
            if not np.all(np.isfinite(k_exact)):
                return 1e10
            rate_exact = beta0 * k_exact + beta2 * g_exact
            lam = hh[ex_cells] * np.bincount(ex_group, weights=rate_exact, minlength=len(ex_cells))
            lam[~np.isfinite(lam)] = 0
            log_p_act = np.sum(np.log(-np.expm1(-np.maximum(lam, 1e-12))))  # log(1 - exp(-lam))

        value = -(log_p_act - log_p_surv)
        return float(value) if np.isfinite(value) else 1e10

    return nll


# ---------------------------------------------------------------------------
# 5. Fitting
# ---------------------------------------------------------------------------

@dataclass
class KernelFit:
    kernel: str
    params: dict          # all parameters, including fixed ones
    free: list            # names of the estimated parameters
    nll: float
    cov: np.ndarray | None  # covariance of the free parameters (None if unavailable)
    converged: bool
    mu: float | None

    @property
    def n_par(self) -> int:
        return len(self.free)

    @property
    def aic(self) -> float:
        return 2 * self.nll + 2 * self.n_par

    def to_dict(self) -> dict:
        """Plain values, for saving as JSON."""
        return {"kernel": self.kernel, "params": self.params, "free": self.free,
                "nll": self.nll, "aic": self.aic, "converged": self.converged, "mu": self.mu,
                "cov": None if self.cov is None else self.cov.tolist()}

    @classmethod
    def from_dict(cls, d: dict) -> "KernelFit":
        cov = None if d["cov"] is None else np.array(d["cov"])
        return cls(d["kernel"], d["params"], d["free"], d["nll"], cov, d["converged"], d["mu"])


def numerical_hessian(f, x: np.ndarray, rel_step: float = 1e-4) -> np.ndarray:
    """Second derivatives of f at x by central differences."""
    n = len(x)
    steps = rel_step * np.maximum(1.0, np.abs(x))
    hess = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            ei = np.zeros(n)
            ej = np.zeros(n)
            ei[i] = steps[i]
            ej[j] = steps[j]
            hess[i, j] = hess[j, i] = (
                f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)
            ) / (4 * steps[i] * steps[j])
    return hess


def fit_kernel(kd: KernelData, kernel: str, hdens=None, estimate_a=None, mu=None,
               start: dict | None = None) -> KernelFit:
    """Maximum-likelihood fit of one kernel by Nelder-Mead on the log-parameters."""
    if estimate_a is None:
        estimate_a = hdens is not None
    start = {**DEFAULT_START, **(start or {})}

    fixed = {}
    if not USES_RHO0[kernel]:
        fixed["log_rho0"] = 0.0
    if not estimate_a:
        fixed["a"] = 0.0
    if mu is None:
        fixed["log_beta2"] = 0.0           # ignored by the likelihood when mu is None
    free = [p for p in PARAM_NAMES if p not in fixed]

    nll = make_nll(kd, kernel, hdens=hdens, estimate_a=estimate_a, mu=mu)

    def objective(theta):
        return nll({**fixed, **dict(zip(free, theta))})

    options = {"maxiter": 5000, "maxfev": 10000, "xatol": 1e-8, "fatol": 1e-10}

    def nelder_mead(x0):
        # R's optim() builds its first simplex with a step of 0.1 * max(|x0|) in
        # every parameter. scipy's default step is only 5% of each value, which
        # explores much less and can slide along ridges R escapes, so copy R.
        size = max(0.1 * np.max(np.abs(x0)), 0.1)
        simplex = np.vstack([x0, x0 + size * np.eye(len(x0))])
        return minimize(objective, x0, method="Nelder-Mead",
                        options={**options, "initial_simplex": simplex})

    # several starting distance scales: the likelihood can have a ridge towards
    # rho0 -> 0 (a pure power law), and a single start can end up on it
    rho0_starts = [5, 20, 80] if USES_RHO0[kernel] else [None]
    result = None
    for rho0 in rho0_starts:
        s = dict(start) if rho0 is None else {**start, "log_rho0": np.log(rho0)}
        r = nelder_mead(np.array([s[p] for p in free]))
        if result is None or r.fun < result.fun:
            result = r
    result = nelder_mead(result.x)          # restart from the best: Nelder-Mead can stall
    x = result.x

    # covariance from the curvature of the likelihood at the optimum
    cov = None
    try:
        cov = np.linalg.inv(numerical_hessian(objective, x))
        if not np.all(np.isfinite(cov)) or np.any(np.diag(cov) <= 0):
            cov = None
    except np.linalg.LinAlgError:
        pass
    if cov is None:
        warnings.warn(f"{kernel}: covariance unavailable (likelihood flat or not at a minimum)")

    params = {**fixed, **dict(zip(free, map(float, x)))}
    return KernelFit(kernel, params, free, float(result.fun), cov, bool(result.success), mu)


def fit_kernels(kd: KernelData, kernels=KERNELS_2P, **kwargs) -> dict[str, KernelFit]:
    """Fit several kernel shapes; kwargs are passed on to fit_kernel()."""
    fits = {}
    for kernel in kernels:
        try:
            fits[kernel] = fit_kernel(kd, kernel, **kwargs)
        except Exception as err:  # one failed shape should not stop the others
            warnings.warn(f"{kernel} failed: {err}")
            continue
        print(f"  {kernel:<7}  nll = {fits[kernel].nll:10.3f}")
    return fits


def compare_kernels(fits: dict[str, KernelFit]) -> pd.DataFrame:
    """AIC table across kernels, best first."""
    rows = []
    for name, f in fits.items():
        p = f.params
        rows.append({
            "kernel": name,
            "label": LABELS[name],
            "npar": f.n_par,
            "nll": f.nll,
            "AIC": f.aic,
            "beta0": np.exp(p["log_beta0"]),
            "rho0": np.exp(p["log_rho0"]) if USES_RHO0[name] else np.nan,
            "psi": np.exp(p["log_psi"]),
            "beta2": np.exp(p["log_beta2"]) if f.mu is not None else np.nan,
            "mu": f.mu if f.mu is not None else np.nan,
            "a": p["a"],
            "converged": f.converged,
        })
    table = pd.DataFrame(rows).sort_values("AIC", ignore_index=True)
    table["dAIC"] = table["AIC"] - table["AIC"].min()
    return table


# ---------------------------------------------------------------------------
# 6. Fitted kernel with uncertainty
# ---------------------------------------------------------------------------

def kernel_curves(fit: KernelFit, x: np.ndarray, n_sim: int = 5000, seed: int = 0) -> np.ndarray:
    """Transmission rate at distances x for parameter sets drawn from the fit's
    approximate sampling distribution (one row per draw). Uses the point
    estimate only if the covariance is unavailable."""
    est = np.array([fit.params[p] for p in fit.free])
    if fit.cov is None:
        draws = est[None, :]
    else:
        draws = np.random.default_rng(seed).multivariate_normal(est, fit.cov, size=n_sim)

    K = KERNELS[fit.kernel]
    curves = np.empty((len(draws), len(x)))
    for i, theta in enumerate(draws):
        p = {**fit.params, **dict(zip(fit.free, theta))}
        rho0 = np.exp(p["log_rho0"]) if USES_RHO0[fit.kernel] else np.nan
        with np.errstate(all="ignore"):
            curves[i] = np.exp(p["log_beta0"]) * K(x, rho0, np.exp(p["log_psi"]))
            if fit.mu is not None:  # plot the combined rate, not K alone
                curves[i] += np.exp(p["log_beta2"]) * g_long(x, fit.mu)
    return curves


def plot_kernel(fit: KernelFit, xlim: float = 200, n_sim: int = 5000, seed: int = 0):
    """Median fitted rate with 50% and 95% bands."""
    import matplotlib.pyplot as plt

    x = np.linspace(0, xlim, 400)
    curves = kernel_curves(fit, x, n_sim=n_sim, seed=seed)
    q = np.nanquantile(curves, [0.025, 0.25, 0.5, 0.75, 0.975], axis=0)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.fill_between(x, q[0], q[4], color="salmon", alpha=0.4, label="95%")
    ax.fill_between(x, q[1], q[3], color="salmon", alpha=0.7, label="50%")
    ax.plot(x, q[2], color="firebrick", label="median")
    ax.set(xlabel="Distance (km)", ylabel="Transmission rate (per cell per step)",
           title=LABELS[fit.kernel])
    ax.legend()
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 6b. National kernel layer
# ---------------------------------------------------------------------------

def predict_kernel_layer(cells_xy: np.ndarray, activated_xy: np.ndarray, fit: KernelFit,
                         normalise: bool = True, include_long: bool = False
                         ) -> tuple[np.ndarray, np.ndarray]:
    """Kernel-only relative risk for every cell: the kernel summed over detected cells.

    cells_xy, activated_xy : two-column arrays of x, y (km)
    include_long : add the long-range component G. Default False: G stands in
        for waterbird movement, which tracks rain and ephemeral water rather
        than distance, so the covariate layers should carry that halo instead.

    Cells containing a detection are masked (returned as NaN) and flagged in
    is_detected, to be mapped separately as "detected". Their value includes
    the kernel at distance 0, K(0), which the data never inform (the shortest
    observed distance is one cell width), so it would otherwise set the scale
    for every other cell. Normalisation therefore uses undetected cells only.

    Returns (z, is_detected).
    """
    if len(activated_xy) == 0:               # no sources: no pressure anywhere
        warnings.warn("No detected cells given: the kernel layer is zero everywhere.")
        return np.zeros(len(cells_xy)), np.zeros(len(cells_xy), dtype=bool)

    p = fit.params
    rho0 = np.exp(p["log_rho0"]) if USES_RHO0[fit.kernel] else np.nan
    psi = np.exp(p["log_psi"])
    d = np.hypot(cells_xy[:, None, 0] - activated_xy[None, :, 0],
                 cells_xy[:, None, 1] - activated_xy[None, :, 1])
    with np.errstate(all="ignore"):
        z = KERNELS[fit.kernel](d, rho0, psi).sum(axis=1)
        if include_long and fit.mu is not None:
            z = np.exp(p["log_beta0"]) * z + np.exp(p["log_beta2"]) * g_long(d, fit.mu).sum(axis=1)

    is_detected = d.min(axis=1) < 1e-9          # cells containing a detection
    if normalise:
        z = z / np.nanmax(z[~is_detected])      # scale set by undetected cells
    z[is_detected] = np.nan
    return z, is_detected


# ---------------------------------------------------------------------------
# 7. Synthetic check: simulate from known parameters, then fit
# ---------------------------------------------------------------------------

def demo(seed: int = 42, grid_n: int = 30, cell_km: float = 25, n_steps: int = 60,
         beta0: float = 0.02, rho0: float = 40, psi: float = 2.5):
    """Simulate spread on a square grid with the fat-tailed kernel, then fit.

    The fitted ftail parameters should come out close to the true ones.

    The simulation uses inf_offset = +1 step (a cell becomes infectious the
    step after it activates), so the simulated process is exactly the one the
    likelihood assumes. With a negative offset the simulator cannot give a
    cell infectiousness before its own activation, but the likelihood assumes
    it does, so the truth would not be recoverable.
    """
    rng = np.random.default_rng(seed)
    gx, gy = np.meshgrid(np.arange(grid_n) * cell_km, np.arange(grid_n) * cell_km)
    x, y = gx.ravel(), gy.ravel()
    n = len(x)
    dist = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])

    inf_offset, inf_period = 1, 14
    act = np.full(n, np.nan)
    act[rng.choice(n, size=2, replace=False)] = 1          # two seeds

    for t in range(2, n_steps + 1):
        s = act + inf_offset
        e = s + inf_period - 1
        infectious = np.flatnonzero((s <= t) & (t <= e))   # comparisons with NaN are False
        if len(infectious) == 0:
            continue
        at_risk = np.flatnonzero(np.isnan(act))
        lam = beta0 * KERNELS["ftail"](dist[np.ix_(at_risk, infectious)], rho0, psi).sum(axis=1)
        act[at_risk[rng.random(len(at_risk)) < 1 - np.exp(-lam)]] = t
    print(f"Simulated: {np.sum(~np.isnan(act))} of {n} cells activated.")

    cells = pd.DataFrame({
        "cell_id": np.arange(n), "x": x, "y": y,
        "detect_date": pd.Timestamp("2026-06-01") + pd.to_timedelta(act - 1, unit="D"),
    })
    kd = build_kernel_data(cells, inf_offset=inf_offset, inf_period=inf_period, bin_km=2)
    fits = fit_kernels(kd)
    print(f"\nTruth: beta0 = {beta0}, rho0 = {rho0}, psi = {psi}\n")
    print(compare_kernels(fits).to_string())
    return kd, fits
"""Multi-run comparison of WRF experiments (a `runs_config` run set) against
GOES / MRMS / HRRR: panel frames, difference maps, and time series.

Three products, all built from the WRF runs' *cheap* fields (no CRTM):

  * `panels_tb` / `panels_refl` -- top row: observations (GOES ch.13 + HRRR
    OLR-Tb, or MRMS + HRRR `refc`) and the reference run; bottom row: the
    other runs.
  * `pair_<field>` -- for every pair `(a, b)` in the config's `pairs`, frames
    of run a | run b | a - b (no observations) for each field: `tb`, `refl`,
    `precip` (15-min rate) and the cloud-water, rain, cloud-ice, snow and
    graupel paths `clwp`/`rwp`/`ciwp`/`swp`/`gwp` (see
    `wrf.read_hydrometeor_paths`). `pair_timeseries` adds domain-mean
    precipitation / path time series with the a - b differences, per pair.
  * `timeseries` -- domain statistics vs. time for every run (and, where
    available, the observations restricted to the WRF domain's footprint): mean/min Tb,
    cold-cloud fraction, reflectivity coverage > 30 / 40 dBZ, and inter-run
    spread. Cached to ``<output-root>/compare/<case>/timeseries_stats.npz``.

**Methodology note (same as the d2 plots, unlike the d1 `tb4`):** the
model Tb here is the OLR-derived one (Yang & Slingo 2001, `OLR = sigma*Tf**4`)
for *both* WRF and HRRR -- no CRTM. GOES ch.13 is a true narrow-band Tb, so
the OLR-derived model Tb runs systematically warm against it. Compare runs
against each other on this scale; use the GOES panel for pattern, not
absolute value.

    python -m bnf_wrf_eval.run_compare --config configs/da_20250520.json \\
        --which panels_tb,panels_refl,pair_tb,pair_refl --cadence 60
    python -m bnf_wrf_eval.run_compare --config ... --which timeseries
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import sys
import time
from pathlib import Path

import numpy as np

from . import goes as goes_reader
from . import hrrr as hrrr_reader
from . import mrms as mrms_reader
from . import plotting as P
from . import wrf as wrf_reader
from .crtm_precompute import _parse_wrfout_time, list_wrfout_times
from .runs_config import RunSet, load_runset

# hydrometeor path kinds -> wrfout mixing-ratio variable
SPECIES = {"clwp": "QCLOUD", "rwp": "QRAIN", "ciwp": "QICE", "swp": "QSNOW", "gwp": "QGRAUP"}
PAIR_FIELDS = ["tb", "refl", "precip"] + list(SPECIES)
PAIR_TYPES = [f"pair_{k}" for k in PAIR_FIELDS] + ["pair_timeseries"]
TYPES = ["panels_tb", "panels_refl", "timeseries"] + PAIR_TYPES
# frame cadence (minutes) per type, else the --cadence option. Every pair field is
# 15-min (for movies); the multi-run panels (with hourly HRRR) are hourly. The
# time-series types read only WRF, so they use every 15-min output.
TYPE_CADENCE = {f"pair_{k}": 15 for k in PAIR_FIELDS}
TYPE_CADENCE.update({"timeseries": 15, "pair_timeseries": 15})
PRECIP_WINDOW_MIN = 15  # precipitation rate = accumulation over this window, scaled to mm/h
COLD_TB_K = 235.0
REFL_THRESHOLDS = (30.0, 40.0)

_RUN_COLORS = ["#1b6ca8", "#d95f02", "#4d9221", "#7b3294", "#666666"]
# one distinct colour per run pair in the inter-run difference panels (up to 6 pairs for 4 runs)
_PAIR_COLORS = ["#1b6ca8", "#d95f02", "#4d9221", "#7b3294", "#a6761d", "#e7298a", "#666666", "#1b9e77"]


# --------------------------------------------------------------- field access
def _wrf_field(kind: str, run_dir: Path, when: dt.datetime, domain: str):
    f = P._find_wrf_file(run_dir, domain, when)
    if kind == "tb":
        return wrf_reader.read_brightness_temperature_from_olr(f)
    return wrf_reader.read_column_max_reflectivity(f)


def _obs_fields(kind: str, when: dt.datetime, opts):
    """[(lon, lat, values, title)] for the observation row."""
    hrrr = P._find_hrrr_file(opts.hrrr_dir, when)
    if kind == "tb":
        gf = P._resolve_goes_file(when, opts.goes_dir, 13, goes_reader.DEFAULT_SATELLITE, False)
        lo, la, v, t = goes_reader.read_brightness_temperature(gf)
        g = (lo, la, v, f"GOES ch.13 Tb (K)  {t:%H:%M} UTC scan")
        lo, la, v, _ = hrrr_reader.read_brightness_temperature_from_olr(hrrr)
        return [g, (lo, la, v, "HRRR OLR-derived Tb (K)")]
    mf = P._resolve_mrms_file(when, opts.mrms_dir, False)
    lo, la, v, t = mrms_reader.read_composite_reflectivity(mf)
    m = (lo, la, v, f"MRMS composite refl. (dBZ)  {t:%H:%M} UTC")
    lo, la, v, _ = hrrr_reader.read_composite_reflectivity(hrrr)
    return [m, (lo, la, v, "HRRR composite refl. (dBZ)")]


def _style(kind: str):
    if kind == "tb":
        return P.DEFAULT_CTT_LEVELS, P.DEFAULT_TB_CMAP, "both", "Brightness Temp. (K)"
    return P.DEFAULT_REFL_LEVELS, "turbo", "max", "Comp. Refl (dBZ)"


def _subplots(nrows: int, ncols: int, proj, panel_w: float, panel_h: float):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(panel_w * ncols, panel_h * nrows),
        subplot_kw={"projection": proj}, constrained_layout=True, squeeze=False,
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)
    return fig, axes


def _domain(rs: RunSet, domain: str, pad_deg: float):
    ref = P._first_wrf_file(rs.ref_dir, domain)
    proj = wrf_reader.get_lambert_projection(ref)
    lon, lat, _, _ = wrf_reader.read_field(ref, "HGT")
    return proj, lon, lat, P._domain_extent(lon, lat, pad_deg), P._domain_outline(lon, lat)


# -------------------------------------------------------------------- panels
def plot_panels(kind: str, when: dt.datetime, rs: RunSet, opts, out_file, pad_deg=0.5):
    import matplotlib.pyplot as plt

    proj, dlon, dlat, extent, outline = _domain(rs, opts.domain, pad_deg)
    levels, cmap, extend, cbar_label = _style(kind)
    obs = _obs_fields(kind, when, opts)
    top = [rs.ref_run]                                   # reference run joins the obs row
    bottom = [n for n in rs.runs if n != rs.ref_run]
    ncols = max(len(obs) + len(top), len(bottom))
    fig, axes = _subplots(2, ncols, proj, 5.2, 3.7)
    for ax in axes.flat:
        ax.set_visible(False)

    mesh = None
    for j, (lo, la, v, title) in enumerate(obs):
        lo, la, v = P._crop_to_extent(lo, la, v, extent)
        axes[0, j].set_visible(True)
        mesh = P._plot_panel(axes[0, j], lo, la, v, levels, cmap, extent, title,
                             extend=extend, domain_outline=outline)
    slots = [(0, len(obs) + j, n) for j, n in enumerate(top)] + [(1, j, n) for j, n in enumerate(bottom)]
    for r, c, name in slots:
        lo, la, v, _ = _wrf_field(kind, rs.runs[name], when, opts.domain)
        axes[r, c].set_visible(True)
        mesh = P._plot_panel(axes[r, c], lo, la, v, levels, cmap, extent, rs.labels[name],
                             extend=extend, domain_outline=outline)
    fig.colorbar(mesh, ax=list(axes.flat), label=cbar_label, shrink=0.6, pad=0.01)
    what = "OLR-derived Tb" if kind == "tb" else "column-max reflectivity"
    fig.suptitle(f"{rs.case[:8]}  {what}  --  {when:%Y-%m-%d %H:%M} UTC", fontsize=15, fontweight="bold")
    fig.savefig(out_file, dpi=110, bbox_inches="tight")
    plt.close(fig)


def _floor0(v):
    """Reflectivity below 0 dBZ (WRF's -30 dBZ 'no echo' fill) -> 0 so
    differences measure echo, not the fill value."""
    return np.maximum(v, 0.0)


# ------------------------------------------------------ precip / water paths
def _trunc_cmap(name: str, lo: float = 0.2, hi: float = 1.0):
    """Sequential colormap without its near-white low end, so the lowest
    contour level is still visible on a white map background."""
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.pyplot as plt

    return LinearSegmentedColormap.from_list(f"{name}_t", plt.get_cmap(name)(np.linspace(lo, hi, 256)))


# kind -> (title, colorbar label, levels, cmap, diff levels, min value shown,
# colorbar extend, stat shown next to the mean). Linear (uniform-width) levels;
# values below `min shown` are left blank so the lowest bin does not paint every
# trace-amount cell. Difference levels are uniform too, with a central bin
# straddling zero. Each species has its own scale (about the domain's 99th
# percentile at the top; snow >> graupel > rain > cloud water >> cloud ice).
_G = "g m$^{-2}$"
_PAIR_KINDS = {
    "tb": ("OLR-derived brightness temperature", "Brightness Temp. (K)",
           P.DEFAULT_CTT_LEVELS, P.DEFAULT_TB_CMAP, np.arange(-32.5, 33, 5), -np.inf, "both", "min"),
    "refl": ("Column-max reflectivity (floored at 0 dBZ)", "Comp. Refl (dBZ)",
             P.DEFAULT_REFL_LEVELS, "turbo", np.arange(-35, 36, 10), -np.inf, "max", "max"),
    "precip": ("Surface precipitation rate (15-min accumulation x 4)", "mm / h",
               np.arange(0, 63, 3), _trunc_cmap("YlGnBu"), np.arange(-27, 28, 6), 0.1, "max", "max"),
    "clwp": ("Cloud liquid water path (CLWP)", _G,
             np.arange(0, 2200, 200), _trunc_cmap("Blues"), np.arange(-900, 901, 200), 10.0, "max", "max"),
    "rwp": ("Rain water path (RWP)", _G,
            np.arange(0, 2750, 250), _trunc_cmap("YlGnBu"), np.arange(-1125, 1126, 250), 10.0, "max", "max"),
    "ciwp": ("Cloud ice water path (CIWP)", _G,
             np.arange(0, 110, 10), _trunc_cmap("BuPu"), np.arange(-45, 46, 10), 1.0, "max", "max"),
    "swp": ("Snow water path (SWP)", _G,
            np.arange(0, 11000, 1000), _trunc_cmap("PuRd"), np.arange(-4500, 4501, 1000), 10.0, "max", "max"),
    "gwp": ("Graupel water path (GWP)", _G,
            np.arange(0, 2200, 200), _trunc_cmap("Oranges"), np.arange(-900, 901, 200), 10.0, "max", "max"),
}


def _pair_field(kind: str, run_dir: Path, when: dt.datetime, domain: str):
    """(lon, lat, values) for one run; units per `_PAIR_KINDS`."""
    f = P._find_wrf_file(run_dir, domain, when)
    if kind == "tb":
        lo, la, tb, _ = wrf_reader.read_brightness_temperature_from_olr(f)
        return lo, la, tb
    if kind == "refl":
        lo, la, refl, _ = wrf_reader.read_column_max_reflectivity(f)
        return lo, la, _floor0(refl)  # WRF's -30 dBZ "no echo" fill -> 0
    if kind == "precip":
        lo, la, acc, _ = wrf_reader.read_accumulated_precipitation(f)
        prev = P._find_wrf_file(run_dir, domain, when - dt.timedelta(minutes=PRECIP_WINDOW_MIN))
        _, _, acc0, _ = wrf_reader.read_accumulated_precipitation(prev)
        # accumulation over the window, scaled to mm/h (15 min -> x4)
        return lo, la, np.maximum(acc - acc0, 0.0) * (60.0 / PRECIP_WINDOW_MIN)
    var = SPECIES[kind]
    lo, la, paths, _ = wrf_reader.read_hydrometeor_paths(f, (var,))
    return lo, la, 1000.0 * paths[var]  # kg/m2 -> g/m2


def plot_pair(kind: str, when: dt.datetime, rs: RunSet, opts, out_file, pair, pad_deg=0.5):
    """Run a | run b | a - b for one field (`PAIR_FIELDS`)."""
    import matplotlib.pyplot as plt

    a, b = pair
    proj, dlon, dlat, extent, outline = _domain(rs, opts.domain, pad_deg)
    title, unit, levels, cmap, dlevels, vmin, extend, stat = _PAIR_KINDS[kind]
    fig, axes = _subplots(1, 3, proj, 5.2, 4.4)
    lo, la, va = _pair_field(kind, rs.runs[a], when, opts.domain)
    _, _, vb = _pair_field(kind, rs.runs[b], when, opts.domain)
    show = lambda v: np.where(v >= vmin, v, np.nan)
    # mean and max (min for Tb) over the whole domain (unmasked field)
    extreme = np.nanmin if stat == "min" else np.nanmax
    label = lambda n, v: f"{rs.labels[n]}\nmean {np.nanmean(v):,.2f}, {stat} {extreme(v):,.2f}"
    m1 = P._plot_panel(axes[0, 0], lo, la, show(va), levels, cmap, extent, label(a, va),
                       extend=extend, domain_outline=outline)
    P._plot_panel(axes[0, 1], lo, la, show(vb), levels, cmap, extent, label(b, vb),
                  extend=extend, domain_outline=outline)
    d = va - vb
    # pcolormesh, not contourf: a difference field is often fine-grained speckle of both
    # signs, and contourf's thousands of tiny polygons then render with the last-drawn
    # (positive) levels visually covering the negative ones -- a one-signed-looking map of
    # a field whose positive and negative parts actually cancel.
    m2 = P._plot_panel(axes[0, 2], lo, la, d, dlevels, "BrBG" if kind not in ("tb", "refl") else "RdBu_r",
                       extent, f"{a} - {b}\nmean {np.nanmean(d):+,.2f}, RMS {np.sqrt(np.nanmean(d**2)):,.2f}",
                       extend="both", domain_outline=outline, plot_type="pcolormesh")
    dunit = "K" if kind == "tb" else ("dBZ" if kind == "refl" else unit)
    fig.colorbar(m1, ax=list(axes[0, :2]), label=unit, shrink=0.8, pad=0.01)
    fig.colorbar(m2, ax=[axes[0, 2]], label=f"difference ({dunit})", shrink=0.8, pad=0.01)
    fig.suptitle(f"{rs.case[:8]}  {title}  --  {when:%Y-%m-%d %H:%M} UTC",
                 fontsize=14, fontweight="bold")
    fig.savefig(out_file, dpi=110, bbox_inches="tight")
    plt.close(fig)


def compute_pair_timeseries(rs: RunSet, opts, times, pair, log=print) -> dict:
    """Domain-mean accumulated precip (mm) and per-species paths (g/m^2) per
    run at every time; columns: `precip`, then `SPECIES` in order."""
    names = list(pair)
    out = {"times": np.array([np.datetime64(t) for t in times])}
    for n in names:
        out[f"{n}/vals"] = np.full((len(times), 1 + len(SPECIES)), np.nan)
    t0 = time.time()
    for i, when in enumerate(times):
        for n in names:
            f = P._find_wrf_file(rs.runs[n], opts.domain, when)
            _, _, acc, _ = wrf_reader.read_accumulated_precipitation(f)
            _, _, paths, _ = wrf_reader.read_hydrometeor_paths(f, tuple(SPECIES.values()))
            out[f"{n}/vals"][i] = [np.nanmean(acc)] + [1000 * np.nanmean(paths[v]) for v in SPECIES.values()]
        if i % 20 == 0:
            log(f"[{i + 1}/{len(times)}] {when:%m-%d %H:%M}  {time.time() - t0:.0f}s")
    return out


def plot_pair_timeseries(stats: dict, rs: RunSet, out_file, pair):
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    a, b = pair
    times = stats["times"].astype("datetime64[s]").astype(object)
    va, vb = stats[f"{a}/vals"], stats[f"{b}/vals"]
    dt_h = np.diff(stats["times"].astype("datetime64[s]").astype(float)) / 3600.0
    rate = lambda v: np.diff(v[:, 0]) / dt_h  # mm/h, between consecutive outputs
    fig, axes = plt.subplots(4, 2, figsize=(14, 15), sharex=True, constrained_layout=True)
    ax = axes.flat
    cols = {a: _RUN_COLORS[0], b: _RUN_COLORS[1]}
    ra, rb = rate(va), rate(vb)
    for n, r in ((a, ra), (b, rb)):
        ax[0].plot(times[1:], r, color=cols[n], lw=1.2, label=rs.labels[n])
    ax[1].plot(times[1:], ra - rb, color="k", lw=1.2)
    ax[0].set_title("Domain-mean precipitation rate (mm/h)", fontsize=11)
    ax[1].set_title(f"Precipitation rate difference, {a} - {b} (mm/h)", fontsize=11)
    ax[0].legend(fontsize=8)
    ax[1].axhline(0, color="0.5", lw=0.8)
    diff_colors = _RUN_COLORS + ["#a6761d"]
    for j, kind in enumerate(SPECIES, start=1):
        title = _PAIR_KINDS[kind][0].split("(")[-1].rstrip(")")
        for n, v in ((a, va), (b, vb)):
            ax[1 + j].plot(times, v[:, j], color=cols[n], lw=1.2)
        ax[1 + j].set_title(f"Domain-mean {title} (g/m$^2$)", fontsize=11)
        ax[7].plot(times, va[:, j] - vb[:, j], color=diff_colors[j - 1], lw=1.2, label=title)
    ax[7].set_title(f"Path differences, {a} - {b} (g/m$^2$)", fontsize=11)
    ax[7].axhline(0, color="0.5", lw=0.8)
    ax[7].legend(fontsize=8, ncol=2)
    for a_ in ax:
        a_.grid(alpha=0.3)
    for a_ in axes[-1]:
        a_.xaxis.set_major_locator(mdates.HourLocator(byhour=[0, 12]))
        a_.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:00"))
        a_.tick_params(axis="x", labelrotation=30)
    fig.suptitle(f"{rs.case[:8]}  {a} vs {b}: precipitation and water paths (domain means)",
                 fontsize=14, fontweight="bold")
    fig.savefig(out_file, dpi=110, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------- time series
def _inside(lo, la, v, proj, box, extent):
    """Values of an obs/HRRR field whose points fall inside the WRF domain's
    Lambert x/y rectangle (the WRF grid is regular in that projection). Area
    statistics are then taken on the obs' own grid -- no regridding, and the
    grids are near-uniform in area (GOES ~2 km, MRMS ~1 km, HRRR 3 km)."""
    lo, la, v = P._crop_to_extent(lo, la, v, extent, pad_cells=0)
    xy = proj.transform_points(P.ccrs.PlateCarree(), lo, la)
    x, y = xy[..., 0], xy[..., 1]
    m = (x >= box[0]) & (x <= box[1]) & (y >= box[2]) & (y <= box[3])
    return v[m]


def _stats(tb, refl):
    return [
        np.nanmean(tb), np.nanmin(tb), np.nanmean(tb < COLD_TB_K),
        np.nanmean(refl > REFL_THRESHOLDS[0]), np.nanmean(refl > REFL_THRESHOLDS[1]),
    ]


STAT_NAMES = ["mean_tb", "min_tb", "cold_frac", "cov30", "cov40"]


def compute_timeseries(rs: RunSet, opts, times, log=print) -> dict:
    """Per-time stats for each run (+ pairwise RMS diffs, + obs inside the WRF domain box)."""
    proj, dlon, dlat, extent, _ = _domain(rs, opts.domain, 0.5)
    xy = proj.transform_points(P.ccrs.PlateCarree(), dlon, dlat)
    box = (xy[..., 0].min(), xy[..., 0].max(), xy[..., 1].min(), xy[..., 1].max())
    names = list(rs.runs)
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]
    out = {"times": np.array([np.datetime64(t) for t in times])}
    for n in names:
        out[f"{n}/stats"] = np.full((len(times), 5), np.nan)
    for a, b in pairs:
        out[f"rms_tb/{a}|{b}"] = np.full(len(times), np.nan)
        out[f"rms_refl/{a}|{b}"] = np.full(len(times), np.nan)
    for o in ("goes", "mrms", "hrrr"):
        out[f"{o}/stats"] = np.full((len(times), 5), np.nan)

    t0 = time.time()
    for i, when in enumerate(times):
        tb, rf = {}, {}
        for n in names:
            try:
                lo, la, tb[n], _ = _wrf_field("tb", rs.runs[n], when, opts.domain)
                _, _, rf[n], _ = _wrf_field("refl", rs.runs[n], when, opts.domain)
            except FileNotFoundError:
                tb.pop(n, None)
                continue
            rf[n] = _floor0(rf[n])
            out[f"{n}/stats"][i] = _stats(tb[n], rf[n])
        for a, b in pairs:
            if a in tb and b in tb:
                out[f"rms_tb/{a}|{b}"][i] = np.sqrt(np.nanmean((tb[a] - tb[b]) ** 2))
                out[f"rms_refl/{a}|{b}"][i] = np.sqrt(np.nanmean((rf[a] - rf[b]) ** 2))
        # obs on the WRF grid (each independent; missing files -> NaN)
        for o, fn in (("goes", lambda w: (P._resolve_goes_file(w, opts.goes_dir, 13, goes_reader.DEFAULT_SATELLITE, False), goes_reader.read_brightness_temperature)),
                      ("mrms", lambda w: (P._resolve_mrms_file(w, opts.mrms_dir, False), mrms_reader.read_composite_reflectivity))):
            try:
                f, reader = fn(when)
                lo2, la2, v2, _ = reader(f)
                g = _inside(lo2, la2, v2, proj, box, extent)
                s = np.full(5, np.nan)
                if o == "goes":
                    s[:3] = [np.nanmean(g), np.nanmin(g), np.nanmean(g < COLD_TB_K)]
                else:
                    g = np.nan_to_num(_floor0(g))
                    s[3:] = [np.mean(g > REFL_THRESHOLDS[0]), np.mean(g > REFL_THRESHOLDS[1])]
                out[f"{o}/stats"][i] = s
            except (FileNotFoundError, KeyError, ValueError):
                pass
        if when.minute == 0:
            try:
                hf = P._find_hrrr_file(opts.hrrr_dir, when)
                lo2, la2, t2, _ = hrrr_reader.read_brightness_temperature_from_olr(hf)
                _, _, r2, _ = hrrr_reader.read_composite_reflectivity(hf)
                out["hrrr/stats"][i] = _stats(_inside(lo2, la2, t2, proj, box, extent),
                                              _floor0(_inside(lo2, la2, r2, proj, box, extent)))
            except FileNotFoundError:
                pass
        if i % 20 == 0:
            log(f"[{i + 1}/{len(times)}] {when:%m-%d %H:%M}  {time.time() - t0:.0f}s")
    return out


def plot_timeseries(stats: dict, rs: RunSet, out_file):
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    # The first output time is the initial-condition snapshot: no radiation has
    # run yet (OLR = 0 -> Tb ~ 0 K), which would wreck every y-axis. Drop it.
    stats = {k: (v[1:] if np.ndim(v) else v) for k, v in stats.items()}
    times = stats["times"].astype("datetime64[s]").astype(object)
    names = list(rs.runs)
    color = {n: _RUN_COLORS[i % len(_RUN_COLORS)] for i, n in enumerate(names)}
    titles = ["Domain-mean OLR-derived Tb (K)", "Domain-min OLR-derived Tb (K)",
              f"Cold-cloud fraction (Tb < {COLD_TB_K:.0f} K)",
              f"Reflectivity coverage > {REFL_THRESHOLDS[0]:.0f} dBZ",
              f"Reflectivity coverage > {REFL_THRESHOLDS[1]:.0f} dBZ"]
    fig, axes = plt.subplots(4, 2, figsize=(15, 13), sharex=True, constrained_layout=True)
    ax = axes.flat
    for k in range(5):
        for j, n in enumerate(names):
            # the reference run (drawn first) is thick + translucent so a
            # run that reproduces it exactly stays visible on top of it
            ax[k].plot(times, stats[f"{n}/stats"][:, k], color=color[n],
                       lw=4.5 if j == 0 else 1.3, alpha=0.4 if j == 0 else 1.0, label=rs.labels[n])
        # obs: GOES only has Tb stats, MRMS only refl; HRRR both (hourly)
        for o, sty, lab in (("goes", ("k", "--"), "GOES (Tb stats: true ch.13, runs warm)"),
                            ("mrms", ("k", "--"), "MRMS"), ("hrrr", ("0.4", ":"), "HRRR (hourly)")):
            y = stats[f"{o}/stats"][:, k]
            if np.isfinite(y).any():
                m = np.isfinite(y)
                ax[k].plot(np.array(times)[m], y[m], color=sty[0], ls=sty[1], lw=1.1,
                           marker="o" if o == "hrrr" else None, ms=3, label=lab)
        ax[k].set_title(titles[k], fontsize=11)
        ax[k].grid(alpha=0.3)
    ax[0].legend(fontsize=7.5, loc="best")
    pair_keys = [k.split("/", 1)[1] for k in stats if k.startswith("rms_tb/")]
    for j, pk in enumerate(pair_keys):
        ax[5].plot(times, stats[f"rms_tb/{pk}"], color=_PAIR_COLORS[j % len(_PAIR_COLORS)], lw=1.2, label=pk.replace("|", " vs "))
        ax[6].plot(times, stats[f"rms_refl/{pk}"], color=_PAIR_COLORS[j % len(_PAIR_COLORS)], lw=1.2, label=pk.replace("|", " vs "))
    ax[5].set_title("Inter-run RMS diff, Tb (K)", fontsize=11)
    ax[6].set_title("Inter-run RMS diff, reflectivity (dBZ, floored at 0)", fontsize=11)
    for a in (ax[5], ax[6]):
        a.set_yscale("symlog", linthresh=0.01)
        a.grid(alpha=0.3)
    ax[5].legend(fontsize=7.5)
    # ensemble range of the cold-cloud fraction across runs
    cf = np.column_stack([stats[f"{n}/stats"][:, 2] for n in names])
    ax[7].plot(times, np.nanmax(cf, 1) - np.nanmin(cf, 1), color="k", lw=1.3)
    ax[7].set_title("Across-run range (max - min) of cold-cloud fraction", fontsize=11)
    ax[7].grid(alpha=0.3)
    for a in axes[-1]:
        a.xaxis.set_major_locator(mdates.HourLocator(byhour=[0, 12]))
        a.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H:00", tz=None))
        a.tick_params(axis="x", labelrotation=30)
    fig.suptitle(f"{rs.case[:8]} run comparison -- domain statistics", fontsize=15, fontweight="bold")
    fig.savefig(out_file, dpi=110, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------- cli
def _times(rs: RunSet, domain: str, cadence: int, start, end) -> list[dt.datetime]:
    """Times present in *every* run (a run still in progress shortens the list)."""
    sets = []
    for d in rs.runs.values():
        ts = {_parse_wrfout_time(f, domain) for f in
              list_wrfout_times(d, domain, every_minutes=None if cadence <= 15 else cadence)}
        sets.append(ts)
    common = sorted(set.intersection(*sets))
    return [t for t in common if (not start or t >= start) and (not end or t <= end)]


def _build_parser():
    p = argparse.ArgumentParser(prog="python -m bnf_wrf_eval.run_compare", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True)
    p.add_argument("--which", default=",".join(TYPES))
    p.add_argument("--cadence", type=int, default=60, help="minutes, for panels/diffs (side/timeseries are always 15)")
    p.add_argument("--domain", default="d01")
    p.add_argument("--output-root", default="outputs")
    p.add_argument("--hrrr-dir", default="satoshi_forcing_data/hrrr/hrrrnat_data")
    p.add_argument("--goes-dir", default="goes_data")
    p.add_argument("--mrms-dir", default="mrms_data")
    p.add_argument("--start", type=dt.datetime.fromisoformat)
    p.add_argument("--end", type=dt.datetime.fromisoformat)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--task-index", type=int)
    p.add_argument("--task-count", type=int)
    p.add_argument("--limit", type=int)
    p.add_argument("--list", action="store_true")
    return p


def main(argv=None) -> int:
    import matplotlib

    matplotlib.use("Agg")
    a = _build_parser().parse_args(argv)
    rs = load_runset(a.config, a.domain)
    types = [t.strip() for t in a.which.split(",") if t.strip()]
    bad = set(types) - set(TYPES)
    if bad:
        raise SystemExit(f"unknown type(s) {sorted(bad)}; choose from {TYPES}")
    root = Path(a.output_root) / "compare" / rs.case
    times_by_cadence: dict[int, list] = {}

    def times_for(t):
        c = TYPE_CADENCE.get(t, a.cadence)
        if c not in times_by_cadence:
            times_by_cadence[c] = _times(rs, a.domain, c, a.start, a.end)
        return times_by_cadence[c]

    times = times_for("timeseries")

    work = []
    for t in types:
        if t in PAIR_TYPES and not rs.pairs:
            continue
        # panels / timeseries: one item each; pair types: one series per configured pair
        for pair in (rs.pairs if t in PAIR_TYPES else [None]):
            tag = f"_{pair[0]}_vs_{pair[1]}" if pair else ""
            if t in ("timeseries", "pair_timeseries"):
                work.append((t, None, root / f"{t}{tag}.png", pair))
                continue
            ts = times_for(t)
            for when in ts:
                # the precip rate needs the output PRECIP_WINDOW_MIN earlier (absent for the
                # first time, and after any gap in a run's output)
                if t == "pair_precip" and when - dt.timedelta(minutes=PRECIP_WINDOW_MIN) not in ts:
                    continue
                out = root / "frames" / t / f"{t}{tag}_{rs.case[:8]}_{when:%Y%m%d_%H%MZ}.png"
                work.append((t, when, out, pair))
    idx = a.task_index if a.task_index is not None else (int(os.environ["SLURM_ARRAY_TASK_ID"]) if "SLURM_ARRAY_TASK_ID" in os.environ else None)
    cnt = a.task_count if a.task_count is not None else (int(os.environ["SLURM_ARRAY_TASK_COUNT"]) if "SLURM_ARRAY_TASK_COUNT" in os.environ else None)
    if idx is not None:
        work = work[idx::cnt]
    if a.limit:
        work = work[: a.limit]
    print(f"run_compare: runs={list(rs.runs)} pairs={rs.pairs} frames={len(work)}", flush=True)
    if a.list:
        for w in work:
            print(w[0], w[1], w[2])
        return 0

    failed = 0
    for i, (t, when, out, extra) in enumerate(work, 1):
        out.parent.mkdir(parents=True, exist_ok=True)
        if t not in ("timeseries", "pair_timeseries") and out.exists() and not a.overwrite:
            print(f"[{i}/{len(work)}] SKIP {out.name}")
            continue
        t0 = time.time()
        try:
            if t == "timeseries":
                cache = root / "timeseries_stats.npz"
                if cache.exists() and not a.overwrite:
                    stats = dict(np.load(cache))
                else:
                    stats = compute_timeseries(rs, a, times)
                    np.savez(cache, **stats)
                plot_timeseries(stats, rs, out)
            elif t == "pair_timeseries":
                cache = out.with_name(out.stem + "_stats.npz")
                if cache.exists() and not a.overwrite:
                    stats = dict(np.load(cache))
                else:
                    stats = compute_pair_timeseries(rs, a, times, extra)
                    np.savez(cache, **stats)
                plot_pair_timeseries(stats, rs, out, extra)
            elif t in PAIR_TYPES:
                plot_pair(t.split("_", 1)[1], when, rs, a, out, extra)
            else:
                plot_panels(t.split("_")[1], when, rs, a, out)
            print(f"[{i}/{len(work)}] OK {out.name} {time.time() - t0:.1f}s", flush=True)
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"[{i}/{len(work)}] FAIL {out.name}: {exc!r}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

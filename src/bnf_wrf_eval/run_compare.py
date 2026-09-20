"""Multi-run comparison of WRF experiments (a `runs_config` run set) against
GOES / MRMS / HRRR: panel frames, difference maps, and time series.

Three products, all built from the WRF runs' *cheap* fields (no CRTM):

  * `panels_tb` / `panels_refl` -- top row: observations (GOES ch.13 + HRRR
    OLR-Tb, or MRMS + HRRR `refc`); bottom row: every run.
  * `diffs_tb` / `diffs_refl` -- one difference map per `diffs` pair in the
    config (``a - b``), one row.
  * `side_tb` / `side_refl` -- two runs (config `side`) side by side, 15-min
    frames, for movies.
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
        --which panels_tb,panels_refl,diffs_tb,diffs_refl --cadence 60
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

TYPES = ["panels_tb", "panels_refl", "diffs_tb", "diffs_refl", "side_tb", "side_refl", "timeseries"]
# frame cadence (minutes) per type; None -> the --cadence option. The side-by-side
# and time-series types read only WRF, so they use every 15-min output.
TYPE_CADENCE = {"side_tb": 15, "side_refl": 15, "timeseries": 15}
TB_DIFF_LEVELS = np.arange(-30, 31, 5)
REFL_DIFF_LEVELS = np.arange(-40, 41, 10)
COLD_TB_K = 235.0
REFL_THRESHOLDS = (30.0, 40.0)

_RUN_COLORS = ["#1b6ca8", "#d95f02", "#4d9221", "#7b3294", "#666666"]


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
    names = list(rs.runs)
    ncols = max(len(names), 2)
    fig, axes = _subplots(2, ncols, proj, 5.2, 3.7)
    for ax in axes.flat:
        ax.set_visible(False)

    mesh = None
    for j, (lo, la, v, title) in enumerate(_obs_fields(kind, when, opts)):
        lo, la, v = P._crop_to_extent(lo, la, v, extent)
        axes[0, j].set_visible(True)
        mesh = P._plot_panel(axes[0, j], lo, la, v, levels, cmap, extent, title,
                             extend=extend, domain_outline=outline)
    for j, name in enumerate(names):
        lo, la, v, _ = _wrf_field(kind, rs.runs[name], when, opts.domain)
        axes[1, j].set_visible(True)
        mesh = P._plot_panel(axes[1, j], lo, la, v, levels, cmap, extent, rs.labels[name],
                             extend=extend, domain_outline=outline)
    fig.colorbar(mesh, ax=list(axes.flat), label=cbar_label, shrink=0.6, pad=0.01)
    what = "OLR-derived Tb" if kind == "tb" else "column-max reflectivity"
    fig.suptitle(f"{rs.case[:8]}  {what}  --  {when:%Y-%m-%d %H:%M} UTC", fontsize=15, fontweight="bold")
    fig.savefig(out_file, dpi=110, bbox_inches="tight")
    plt.close(fig)


def plot_side(kind: str, when: dt.datetime, rs: RunSet, opts, out_file, pad_deg=0.5):
    import matplotlib.pyplot as plt

    proj, dlon, dlat, extent, outline = _domain(rs, opts.domain, pad_deg)
    levels, cmap, extend, cbar_label = _style(kind)
    fig, axes = _subplots(1, 2, proj, 5.6, 4.6)
    mesh = None
    for ax, name in zip(axes.flat, rs.side):
        lo, la, v, _ = _wrf_field(kind, rs.runs[name], when, opts.domain)
        mesh = P._plot_panel(ax, lo, la, v, levels, cmap, extent,
                             rs.side_labels.get(name, rs.labels[name]),
                             extend=extend, domain_outline=outline)
    fig.colorbar(mesh, ax=list(axes.flat), label=cbar_label, shrink=0.8, pad=0.01)
    what = "OLR-derived Tb" if kind == "tb" else "column-max reflectivity"
    fig.suptitle(f"{rs.case[:8]}  {what}  --  {when:%Y-%m-%d %H:%M} UTC", fontsize=14, fontweight="bold")
    fig.savefig(out_file, dpi=110, bbox_inches="tight")
    plt.close(fig)


def _floor0(v):
    """Reflectivity below 0 dBZ (WRF's -30 dBZ 'no echo' fill) -> 0 so
    differences measure echo, not the fill value."""
    return np.maximum(v, 0.0)


def plot_diffs(kind: str, when: dt.datetime, rs: RunSet, opts, out_file, pad_deg=0.5):
    import matplotlib.pyplot as plt

    proj, dlon, dlat, extent, outline = _domain(rs, opts.domain, pad_deg)
    levels = TB_DIFF_LEVELS if kind == "tb" else REFL_DIFF_LEVELS
    n = len(rs.diffs)
    fig, axes = _subplots(1, n, proj, 5.4, 4.6)
    mesh = None
    stats = []
    for ax, (a, b) in zip(axes.flat, rs.diffs):
        lo, la, va, _ = _wrf_field(kind, rs.runs[a], when, opts.domain)
        _, _, vb, _ = _wrf_field(kind, rs.runs[b], when, opts.domain)
        if kind == "refl":
            va, vb = _floor0(va), _floor0(vb)
        d = va - vb
        rms = float(np.sqrt(np.nanmean(d**2)))
        mesh = P._plot_panel(ax, lo, la, d, levels, "RdBu_r", extent,
                             f"{a} - {b}\nRMS {rms:.2f}, max|d| {np.nanmax(np.abs(d)):.1f}",
                             extend="both", domain_outline=outline)
    unit = "K" if kind == "tb" else "dBZ"
    fig.colorbar(mesh, ax=list(axes.flat), label=f"difference ({unit})", shrink=0.7, pad=0.01)
    what = "OLR-derived Tb" if kind == "tb" else "column-max reflectivity (floored at 0 dBZ)"
    fig.suptitle(f"{rs.case[:8]}  {what} differences  --  {when:%Y-%m-%d %H:%M} UTC",
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
        ax[5].plot(times, stats[f"rms_tb/{pk}"], color=_RUN_COLORS[j % 5], lw=1.2, label=pk.replace("|", " vs "))
        ax[6].plot(times, stats[f"rms_refl/{pk}"], color=_RUN_COLORS[j % 5], lw=1.2, label=pk.replace("|", " vs "))
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
        if t == "timeseries":
            work.append((t, None, root / "timeseries.png"))
            continue
        if (t.startswith("diffs") and not rs.diffs) or (t.startswith("side") and not rs.side):
            continue
        for when in times_for(t):
            work.append((t, when, root / "frames" / t / f"{t}_{rs.case[:8]}_{when:%Y%m%d_%H%MZ}.png"))
    idx = a.task_index if a.task_index is not None else (int(os.environ["SLURM_ARRAY_TASK_ID"]) if "SLURM_ARRAY_TASK_ID" in os.environ else None)
    cnt = a.task_count if a.task_count is not None else (int(os.environ["SLURM_ARRAY_TASK_COUNT"]) if "SLURM_ARRAY_TASK_COUNT" in os.environ else None)
    if idx is not None:
        work = work[idx::cnt]
    if a.limit:
        work = work[: a.limit]
    print(f"run_compare: runs={list(rs.runs)} diffs={rs.diffs} frames={len(work)}", flush=True)
    if a.list:
        for w in work:
            print(w[0], w[1], w[2])
        return 0

    failed = 0
    for i, (t, when, out) in enumerate(work, 1):
        out.parent.mkdir(parents=True, exist_ok=True)
        if t != "timeseries" and out.exists() and not a.overwrite:
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
            elif t.startswith("side"):
                plot_side(t.split("_")[1], when, rs, a, out)
            elif t.startswith("panels"):
                plot_panels(t.split("_")[1], when, rs, a, out)
            else:
                plot_diffs(t.split("_")[1], when, rs, a, out)
            print(f"[{i}/{len(work)}] OK {out.name} {time.time() - t0:.1f}s", flush=True)
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"[{i}/{len(work)}] FAIL {out.name}: {exc!r}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

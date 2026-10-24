"""Batch rendering of the comparison and single-panel plots.

Separates plotting (many independent, ~1-3 s cartopy renders) from movie
assembly (`bnf_wrf_eval.make_movies`), the same way `crtm_precompute` split
the CRTM forward model from plotting. Run this to render every frame of

  * the hourly 4-panel Tb comparison, per run, in two variants -- `tb4`
    (WRF ctt | WRF CRTM | HRRR | GOES) and `tb4_simple` (WRF OLR-fit Tb |
    WRF CRTM | HRRR | GOES) -- and the 3-panel reflectivity comparison
    (WRF | HRRR | MRMS), per run; and
  * the six single-panel series -- WRF CRTM Tb, WRF reflectivity (both
    15-min, per run), HRRR Tb, HRRR reflectivity (both hourly, per case),
    GOES Tb (5-min, per case), MRMS reflectivity (15-min, per case).

Runs as a module (like `crtm_precompute`, deliberately no console-script):

    python -m bnf_wrf_eval.batch_plot --which comparison --output-root outputs
    python -m bnf_wrf_eval.batch_plot --which single --output-root outputs

`--which` takes `comparison`, `single`, `all`, or any comma-separated list
of the individual type names (`tb4`, `tb4_simple`, `refl3`, `wrf_tb`,
`wrf_refl`, `hrrr_tb`, `hrrr_refl`, `goes_tb`, `mrms_refl`).

The flat work list is deterministically ordered, so `--task-index` /
`--task-count` (or the `SLURM_ARRAY_TASK_ID` / `SLURM_ARRAY_TASK_COUNT` a
job array sets) slice it `worklist[idx::cnt]` -- strided, so each task gets
an even mix of types / case days / times. Frames already on disk are
skipped (restartable); per-frame exceptions are retried then logged and
skipped so one bad frame does not abort the batch.

Frames land in `<output-root>/frames/<type>/`. The observation files GOES
Tb and MRMS reflectivity need must already be local (run
`scripts/download_obs.py` first); pass `--auto-download` to fetch missing
ones inline instead (races if many tasks want the same missing file, so
prefer the pre-download).

STATUS (2026-08-30): new. The per-frame plotting functions it drives live
in `bnf_wrf_eval.plotting` and are the same ones the notebooks use.
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import sys
import time
from pathlib import Path

from .crtm_precompute import list_wrfout_times
from .plotting import _case_output_path, _output_path
from .runs_config import load_runset

DEFAULT_DOMAIN = "d01"
DEFAULT_RUNS_ROOT = "satoshi_testruns"
DEFAULT_OUTPUT_ROOT = "outputs"
DEFAULT_HRRR_DIR = "satoshi_forcing_data/hrrr/hrrrnat_data"
DEFAULT_GOES_DIR = "goes_data"
DEFAULT_MRMS_DIR = "mrms_data"
DEFAULT_CACHE_DIR = "crtm_cache"

# The nine d1 HRRR3 runs: three case days x three microphysics schemes.
CASES = [
    "20250502lassobnfwrfhrrr3",
    "20250520lassobnfwrfhrrr3",
    "20250917lassobnfwrfhrrr3",
]
SCHEMES = {
    "rund1": "2.5-km (D1) Thompson",
    "rund1-mp52": "2.5-km (D1) P3 (mp=52)",
    "rund1-mp53": "2.5-km (D1) P3 (mp=53)",
    "rund1-dynlit": "2.5-km (D1) Thompson + dyn-lightning",
}
REF_SCHEME = "rund1"  # subdir used for the (time-independent) domain geometry

# The 500-m ndown runs, for the "d2run"-scope 4-panel comparisons. Each d2
# run's d1 (2.5-km parent) panel is resolved by
# `plotting._D2_TO_D1_RUN` (rund2 -> rund1, rund2-dynlit -> rund1-dynlit).
D2_SCHEMES = {
    "rund2": "500-m (D2) Thompson",
    "rund2-dynlit": "500-m (D2) Thompson + dyn-lightning",
}

# The 500-m + 100-m (online one-way nest) runs, for the "d3run"-scope
# 4-panel comparisons. WRF numbers domains within the run, so the 100-m
# "d3" is `wrfout_d02_*` and its 500-m parent "d2" is `wrfout_d01_*`, both
# in the same run directory (see `plotting.plot_run_d3_*`). Frame times are
# the hourly d3 times that also have a d2 file (a few d2 hours are missing
# at restart boundaries).
D3_SCHEMES = {
    "rund2d3": "100-m (D3) Thompson",
}
D3_DOMAIN = "d02"
D3_PARENT_DOMAIN = "d01"

# The 500-m runs driven directly by hourly HRRRv4 analyses, for the
# "d2hrrr"-scope 4-panel comparisons against the free-run-driven (ndown from
# d1) `*hrrr3/rund2` on the same grid -- resolved by
# `plotting._HRRR_DRIVEN_TO_FREE_CASE` (hrrr4 -> hrrr3). Frame times are the
# hourly times both runs have.
D2HRRR_CASES = [
    "20250502lassobnfwrfhrrr4",
    "20250520lassobnfwrfhrrr4",
    "20250917lassobnfwrfhrrr4",
]
D2HRRR_SCHEMES = {
    "rund2": "500-m (D2) Thompson: free-run- vs HRRR-driven",
}

# Only `rund1-dynlit` has `dyn_lightning_option=1` (LPOS/LNEG/LNEU); the
# other three `SCHEMES` entries would raise KeyError from
# `wrf.read_dyn_lightning_flash_counts` immediately. A separate one-entry
# schemes dict (matching the `D2_SCHEMES` pattern above) keeps the
# "lightning_*" types from wastefully attempting -- and failing on --
# rund1/rund1-mp52/rund1-mp53.
DYNLIT_SCHEMES = {"rund1-dynlit": SCHEMES["rund1-dynlit"]}

# Types needing a *pair* of consecutive-hour wrfout times (the difference
# of LPOS/LNEG/LNEU between them), not a single time -- so the run's very
# last hourly timestamp (no successor) must be excluded from the
# worklist, unlike every other "run"-scope type.
PAIRED_TYPES = {"lightning_cg", "lightning_total"}

# type -> (scope, cadence-minutes, output-filename prefix). scope "run"
# iterates every case x d1 scheme; "case" iterates case days only; "d2run"
# iterates every case x d2 scheme (D2_SCHEMES); "d3run" every case x d3
# scheme (D3_SCHEMES); "d2hrrr" every HRRR-driven case (D2HRRR_CASES) x
# D2HRRR_SCHEMES.
TYPES: dict[str, tuple[str, int, str]] = {
    "tb4": ("run", 60, "wrf_ctt_crtm_hrrr_goes_tb_comparison"),
    "tb4_simple": ("run", 60, "wrf_simple_tb_crtm_hrrr_goes_tb_comparison"),
    "refl3": ("run", 60, "wrf_hrrr_mrms_refl_comparison"),
    "wrf_tb": ("run", 15, "wrf_crtm_tb"),
    "wrf_refl": ("run", 15, "wrf_refl"),
    "wrf_olr_tb": ("run", 15, "wrf_olr_tb"),
    "hrrr_tb": ("case", 60, "hrrr_tb"),
    "hrrr_refl": ("case", 60, "hrrr_refl"),
    "goes_tb": ("case", 5, "goes_tb"),
    "mrms_refl": ("case", 15, "mrms_refl"),
    # d2 (500-m) 4-panel comparisons: d2 | d1 | HRRR | obs, hourly, all
    # Tb panels OLR-derived (no CRTM -- unlike the d1 `tb4`/`wrf_tb`).
    "d2_tb4": ("d2run", 60, "d2_olr_tb_4panel"),
    "d2_refl4": ("d2run", 60, "d2_refl_4panel"),
    # d3 (100-m) 4-panel comparisons: d3 | d2 | HRRR | obs, hourly, Tb
    # OLR-derived for every model panel like d2.
    "d3_tb4": ("d3run", 60, "d3_olr_tb_4panel"),
    "d3_refl4": ("d3run", 60, "d3_refl_4panel"),
    # d2 (500-m) free-run-driven | HRRR-driven | HRRR | obs, hourly, Tb
    # OLR-derived for every model panel like d2.
    "d2hrrr_tb4": ("d2hrrr", 60, "d2_free_vs_hrrr_olr_tb_4panel"),
    "d2hrrr_refl4": ("d2hrrr", 60, "d2_free_vs_hrrr_refl_4panel"),
    # WRF (rund1-dynlit only) vs. observed lightning, hourly-accumulated
    # flash counts. scope "dynlit" -> DYNLIT_SCHEMES (just rund1-dynlit);
    # `when` is each window's *start* time (see PAIRED_TYPES above).
    "lightning_cg": ("dynlit", 60, "lightning_cg_wrf_mrms"),
    "lightning_total": ("dynlit", 60, "lightning_total_wrf_glm"),
}
GROUPS = {
    "comparison": ["tb4", "tb4_simple", "refl3"],
    "single": ["wrf_tb", "wrf_refl", "wrf_olr_tb", "hrrr_tb", "hrrr_refl", "goes_tb", "mrms_refl"],
    "d2": ["d2_tb4", "d2_refl4"],
    "d3": ["d3_tb4", "d3_refl4"],
    "d2hrrr": ["d2hrrr_tb4", "d2hrrr_refl4"],
    "lightning": ["lightning_cg", "lightning_total"],
    "all": list(TYPES),
}


# Filled by `apply_runs_config` (`--runs-config`): explicit run dirs that
# replace the built-in `<runs-root>/<case>/<scheme>` layout.
RUN_DIRS: dict[tuple[str, str], Path] = {}
REF_DIRS: dict[str, Path] = {}
# Types that read the CRTM cache (would compute inline for unknown runs).
CRTM_TYPES = {"tb4", "wrf_tb"}


def _run_dir(runs_root: Path, case: str, scheme: str) -> Path:
    return RUN_DIRS.get((case, scheme), runs_root / case / scheme)


def _ref_dir(runs_root: Path, case: str) -> Path:
    return REF_DIRS.get(case, runs_root / case / REF_SCHEME)


def apply_runs_config(path: str | Path, domain: str = DEFAULT_DOMAIN) -> None:
    """Replace CASES/SCHEMES with the runs of a `runs_config` JSON file.

    Scheme keys are the run dirs' basenames (what `_render` looks labels up
    by), so they must be unique. The d2 / d3 / dyn-lightning scopes are
    emptied -- they assume the built-in layout and dyn-lightning output."""
    global CASES, SCHEMES, D2_SCHEMES, D3_SCHEMES, D2HRRR_SCHEMES, DYNLIT_SCHEMES
    rs = load_runset(path, domain)
    CASES = [rs.case]
    SCHEMES, D2_SCHEMES, D3_SCHEMES, D2HRRR_SCHEMES, DYNLIT_SCHEMES = {}, {}, {}, {}, {}
    RUN_DIRS.clear()
    for name, d in rs.runs.items():
        if d.name in SCHEMES:
            raise SystemExit(f"run dirs share basename {d.name!r}; batch_plot keys labels on it")
        SCHEMES[d.name] = rs.labels[name]
        RUN_DIRS[(rs.case, d.name)] = d
    REF_DIRS.clear()
    REF_DIRS[rs.case] = rs.ref_dir


def _resolve_types(which: str) -> list[str]:
    names: list[str] = []
    for token in which.split(","):
        token = token.strip()
        if not token:
            continue
        if token in GROUPS:
            names.extend(GROUPS[token])
        elif token in TYPES:
            names.append(token)
        else:
            raise SystemExit(f"unknown --which token {token!r}")
    # de-dupe, keep order
    seen: set[str] = set()
    return [n for n in names if not (n in seen or seen.add(n))]


def _run_span(run_dir: Path, domain: str) -> tuple[dt.datetime, dt.datetime]:
    files = list_wrfout_times(run_dir, domain)
    if not files:
        raise FileNotFoundError(f"no wrfout_{domain}_* in {run_dir}")
    from .crtm_precompute import _parse_wrfout_time

    times = [_parse_wrfout_time(f, domain) for f in files]
    return min(times), max(times)


def _time_grid(start: dt.datetime, end: dt.datetime, step_minutes: int) -> list[dt.datetime]:
    out, t = [], start
    step = dt.timedelta(minutes=step_minutes)
    while t <= end:
        out.append(t)
        t += step
    return out


def _wrf_times(run_dir: Path, domain: str, cadence_minutes: int) -> list[dt.datetime]:
    from .crtm_precompute import _parse_wrfout_time

    every = None if cadence_minutes <= 15 else cadence_minutes
    files = list_wrfout_times(run_dir, domain, every_minutes=every)
    return [_parse_wrfout_time(f, domain) for f in files]


def build_worklist(
    types: list[str],
    *,
    runs_root: str | Path = DEFAULT_RUNS_ROOT,
    output_root: str | Path = DEFAULT_OUTPUT_ROOT,
    domain: str = DEFAULT_DOMAIN,
    start: dt.datetime | None = None,
    end: dt.datetime | None = None,
) -> list[tuple]:
    """Flat, deterministically ordered list of frames to render.

    Each item is `(ptype, key, when, out_path)` where `key` is the run
    directory (scope "run") or the case label (scope "case")."""
    runs_root = Path(runs_root)
    frames_root = Path(output_root) / "frames"
    items: list[tuple] = []

    for ptype in types:
        scope, cadence, prefix = TYPES[ptype]
        out_dir = frames_root / ptype

        if scope in ("run", "d2run", "d3run", "d2hrrr", "dynlit"):
            schemes = {
                "run": SCHEMES, "d2run": D2_SCHEMES, "d3run": D3_SCHEMES,
                "d2hrrr": D2HRRR_SCHEMES, "dynlit": DYNLIT_SCHEMES,
            }[scope]
            cases = D2HRRR_CASES if scope == "d2hrrr" else CASES
            for case in cases:
                for scheme in schemes:
                    run_dir = _run_dir(runs_root, case, scheme)
                    if scope == "d3run":
                        parent = set(_wrf_times(run_dir, D3_PARENT_DOMAIN, cadence))
                        times = [t for t in _wrf_times(run_dir, D3_DOMAIN, cadence) if t in parent]
                    elif scope == "d2hrrr":
                        from .plotting import _resolve_free_d2_run_dir

                        free = set(_wrf_times(_resolve_free_d2_run_dir(run_dir, None), domain, cadence))
                        times = [t for t in _wrf_times(run_dir, domain, cadence) if t in free]
                    else:
                        times = _wrf_times(run_dir, domain, cadence)
                    if ptype in PAIRED_TYPES and times:
                        times = times[:-1]  # drop the last hour: no successor to pair with
                    for when in times:
                        if start and when < start:
                            continue
                        if end and when > end:
                            continue
                        out_path = _output_path(out_dir, prefix, run_dir, when)
                        items.append((ptype, run_dir, when, out_path))
        else:  # scope == "case"
            for case in CASES:
                ref_dir = _ref_dir(runs_root, case)
                if ptype == "goes_tb":
                    # 5-min regular grid, snapped to the nearest ABI scan
                    # per frame -- there is no wrfout at :05/:10/... .
                    span_start, span_end = _run_span(ref_dir, domain)
                    grid = _time_grid(span_start, span_end, cadence)
                else:
                    # HRRR (hourly) and MRMS (15-min): reuse the run's own
                    # wrfout times so frame counts line up with the WRF
                    # single-panel series and every time has a match.
                    grid = _wrf_times(ref_dir, domain, cadence)
                for when in grid:
                    if start and when < start:
                        continue
                    if end and when > end:
                        continue
                    out_path = _case_output_path(out_dir, prefix, case, when)
                    items.append((ptype, case, when, out_path))

    items.sort(key=lambda it: (it[0], str(it[1]), it[2]))
    return items


def _render(item: tuple, opts: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from . import plotting

    ptype, key, when, out_path = item
    runs_root = Path(opts.runs_root)
    try:
        if ptype in ("tb4", "tb4_simple"):
            plotting.plot_run_tb_comparison_4panel(
                when, SCHEMES[Path(key).name], key,
                hrrr_base_dir=opts.hrrr_dir, goes_dir=opts.goes_dir,
                domain=opts.domain, out_file=out_path,
                auto_download_goes=opts.auto_download,
                wrf_panel="olr" if ptype == "tb4_simple" else "ctt",
                wrf_tb_kwargs={"cache_dir": opts.cache_dir, "require_cache": True},
            )
        elif ptype == "refl3":
            plotting.plot_run_refl_comparison(
                when, SCHEMES[Path(key).name], key,
                hrrr_base_dir=opts.hrrr_dir, mrms_dir=opts.mrms_dir,
                domain=opts.domain, out_file=out_path,
                auto_download_mrms=opts.auto_download,
            )
        elif ptype == "wrf_tb":
            plotting.plot_run_wrf_tb_single(
                when, SCHEMES[Path(key).name], key, domain=opts.domain,
                out_file=out_path,
                wrf_tb_kwargs={"cache_dir": opts.cache_dir, "require_cache": True},
            )
        elif ptype == "wrf_refl":
            plotting.plot_run_wrf_refl_single(
                when, SCHEMES[Path(key).name], key, domain=opts.domain,
                out_file=out_path,
            )
        elif ptype == "wrf_olr_tb":
            plotting.plot_run_wrf_olr_tb_single(
                when, SCHEMES[Path(key).name], key, domain=opts.domain,
                out_file=out_path,
            )
        elif ptype == "d2_tb4":
            plotting.plot_run_d2_tb_comparison_4panel(
                when, D2_SCHEMES[Path(key).name], key,
                hrrr_base_dir=opts.hrrr_dir, goes_dir=opts.goes_dir,
                domain=opts.domain, out_file=out_path,
                auto_download_goes=opts.auto_download,
            )
        elif ptype == "d2_refl4":
            plotting.plot_run_d2_refl_comparison_4panel(
                when, D2_SCHEMES[Path(key).name], key,
                hrrr_base_dir=opts.hrrr_dir, mrms_dir=opts.mrms_dir,
                domain=opts.domain, out_file=out_path,
                auto_download_mrms=opts.auto_download,
            )
        elif ptype == "d3_tb4":
            plotting.plot_run_d3_tb_comparison_4panel(
                when, D3_SCHEMES[Path(key).name], key,
                hrrr_base_dir=opts.hrrr_dir, goes_dir=opts.goes_dir,
                d3_domain=D3_DOMAIN, d2_domain=D3_PARENT_DOMAIN, out_file=out_path,
                auto_download_goes=opts.auto_download,
            )
        elif ptype == "d3_refl4":
            plotting.plot_run_d3_refl_comparison_4panel(
                when, D3_SCHEMES[Path(key).name], key,
                hrrr_base_dir=opts.hrrr_dir, mrms_dir=opts.mrms_dir,
                d3_domain=D3_DOMAIN, d2_domain=D3_PARENT_DOMAIN, out_file=out_path,
                auto_download_mrms=opts.auto_download,
            )
        elif ptype == "d2hrrr_tb4":
            plotting.plot_run_d2_forcing_tb_comparison_4panel(
                when, D2HRRR_SCHEMES[Path(key).name], key,
                hrrr_base_dir=opts.hrrr_dir, goes_dir=opts.goes_dir,
                domain=opts.domain, out_file=out_path,
                auto_download_goes=opts.auto_download,
            )
        elif ptype == "d2hrrr_refl4":
            plotting.plot_run_d2_forcing_refl_comparison_4panel(
                when, D2HRRR_SCHEMES[Path(key).name], key,
                hrrr_base_dir=opts.hrrr_dir, mrms_dir=opts.mrms_dir,
                domain=opts.domain, out_file=out_path,
                auto_download_mrms=opts.auto_download,
            )
        elif ptype == "lightning_cg":
            plotting.plot_run_lightning_cg_comparison(
                when, DYNLIT_SCHEMES[Path(key).name], key,
                mrms_dir=opts.mrms_dir, domain=opts.domain, out_file=out_path,
                auto_download_mrms=opts.auto_download,
            )
        elif ptype == "lightning_total":
            plotting.plot_run_lightning_total_comparison(
                when, DYNLIT_SCHEMES[Path(key).name], key,
                goes_dir=opts.goes_dir, domain=opts.domain, out_file=out_path,
                auto_download_glm=opts.auto_download,
            )
        elif ptype == "hrrr_tb":
            plotting.plot_case_hrrr_tb_single(
                when, key, _ref_dir(runs_root, key),
                hrrr_base_dir=opts.hrrr_dir, domain=opts.domain, out_file=out_path,
            )
        elif ptype == "hrrr_refl":
            plotting.plot_case_hrrr_refl_single(
                when, key, _ref_dir(runs_root, key),
                hrrr_base_dir=opts.hrrr_dir, domain=opts.domain, out_file=out_path,
            )
        elif ptype == "goes_tb":
            plotting.plot_case_goes_tb_single(
                when, key, _ref_dir(runs_root, key),
                goes_dir=opts.goes_dir, domain=opts.domain, out_file=out_path,
                auto_download_goes=opts.auto_download,
            )
        elif ptype == "mrms_refl":
            plotting.plot_case_mrms_refl_single(
                when, key, _ref_dir(runs_root, key),
                mrms_dir=opts.mrms_dir, domain=opts.domain, out_file=out_path,
                auto_download_mrms=opts.auto_download,
            )
        else:  # pragma: no cover - guarded by _resolve_types
            raise ValueError(f"unknown plot type {ptype!r}")
    finally:
        plt.close("all")


def run(worklist: list[tuple], opts: argparse.Namespace, log=print) -> dict:
    summary = {"ok": 0, "skipped": 0, "failed": 0, "failures": [], "seconds": 0.0}
    t_start = time.time()
    total = len(worklist)

    for i, item in enumerate(worklist, 1):
        ptype, key, when, out_path = item
        tag = f"{ptype} {Path(key).name if isinstance(key, Path) else key} {when:%Y-%m-%d_%H:%M}"
        out_path = Path(out_path)
        if out_path.exists() and not opts.overwrite:
            summary["skipped"] += 1
            log(f"[{i}/{total}] SKIP  {tag}")
            continue

        for attempt in range(opts.retry + 1):
            t0 = time.time()
            try:
                out_path.parent.mkdir(parents=True, exist_ok=True)
                _render(item, opts)
                summary["ok"] += 1
                log(f"[{i}/{total}] OK    {tag}  {time.time() - t0:.1f}s  -> {out_path.name}")
                break
            except Exception as exc:  # noqa: BLE001 - one bad frame must not abort the batch
                if attempt < opts.retry:
                    log(f"[{i}/{total}] RETRY {tag}  attempt {attempt + 1}: {exc!r}")
                    continue
                summary["failed"] += 1
                summary["failures"].append(f"{tag}: {exc!r}")
                log(f"[{i}/{total}] FAIL  {tag}  after {attempt + 1} attempt(s): {exc!r}")

    summary["seconds"] = time.time() - t_start
    return summary


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m bnf_wrf_eval.batch_plot",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--which", default="all",
                   help="comparison | single | all | comma-separated type names")
    p.add_argument("--runs-root", default=DEFAULT_RUNS_ROOT)
    p.add_argument("--runs-config", default=None,
                   help="JSON run-set (configs/*.json) replacing CASES/SCHEMES; "
                        "see bnf_wrf_eval.runs_config. Not valid for CRTM-based types.")
    p.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    p.add_argument("--domain", default=DEFAULT_DOMAIN)
    p.add_argument("--hrrr-dir", default=DEFAULT_HRRR_DIR)
    p.add_argument("--goes-dir", default=DEFAULT_GOES_DIR)
    p.add_argument("--mrms-dir", default=DEFAULT_MRMS_DIR)
    p.add_argument("--cache-dir", default=DEFAULT_CACHE_DIR)
    p.add_argument("--start", type=dt.datetime.fromisoformat, default=None)
    p.add_argument("--end", type=dt.datetime.fromisoformat, default=None)
    p.add_argument("--auto-download", action="store_true",
                   help="fetch missing GOES/MRMS files inline (races under a job array)")
    p.add_argument("--overwrite", action="store_true",
                   help="re-render frames that already exist on disk")
    p.add_argument("--retry", type=int, default=1,
                   help="extra attempts per frame on an exception")
    p.add_argument("--task-index", type=int, default=None)
    p.add_argument("--task-count", type=int, default=None)
    p.add_argument("--limit", type=int, default=None,
                   help="render at most this many frames (after task slicing)")
    p.add_argument("--list", action="store_true",
                   help="print the resolved work list and exit")
    return p


def _resolve_task_slice(args: argparse.Namespace) -> tuple[int | None, int | None]:
    idx, cnt = args.task_index, args.task_count
    if idx is None and "SLURM_ARRAY_TASK_ID" in os.environ:
        idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
    if cnt is None and "SLURM_ARRAY_TASK_COUNT" in os.environ:
        cnt = int(os.environ["SLURM_ARRAY_TASK_COUNT"])
    if (idx is None) != (cnt is None):
        raise SystemExit("task-index and task-count must be given together")
    return idx, cnt


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    types = _resolve_types(args.which)
    if args.runs_config:
        dropped = [t for t in types if t in CRTM_TYPES]
        if dropped:
            print(f"--runs-config: dropping CRTM-based types {dropped} (no CRTM cache for these runs)")
            types = [t for t in types if t not in CRTM_TYPES]
        apply_runs_config(args.runs_config, args.domain)

    worklist = build_worklist(
        types,
        runs_root=args.runs_root,
        output_root=args.output_root,
        domain=args.domain,
        start=args.start,
        end=args.end,
    )
    idx, cnt = _resolve_task_slice(args)
    scope = ""
    if idx is not None:
        worklist = worklist[idx::cnt]
        scope = f" (task {idx}/{cnt})"
    if args.limit is not None:
        worklist = worklist[: args.limit]
        scope += f" (limit {args.limit})"

    if args.list:
        for j, (ptype, key, when, out_path) in enumerate(worklist):
            name = Path(key).name if isinstance(key, Path) else key
            print(f"{j:5d}  {ptype:10s}  {name:24s}  {when:%Y-%m-%d_%H:%M}  {out_path}")
        print(f"# {len(worklist)} frame(s){scope}  types={','.join(types)}")
        return 0

    print(
        f"batch_plot: {len(worklist)} frame(s){scope}  types={','.join(types)}  "
        f"output_root={args.output_root}  auto_download={args.auto_download}",
        flush=True,
    )
    summary = run(worklist, args)
    print(
        f"batch_plot done{scope}: ok={summary['ok']} skipped={summary['skipped']} "
        f"failed={summary['failed']}  {summary['seconds'] / 60:.1f} min",
        flush=True,
    )
    for failure in summary["failures"]:
        print(f"  FAILED: {failure}", flush=True)
    return 1 if summary["failed"] else 0


if __name__ == "__main__":
    sys.exit(main())

"""JSON run-set config: point the batch plotters at run directories that
are not the built-in `satoshi_testruns/<case>/<scheme>` layout (e.g. the
`bnf_wrf_da` experiments).

Schema (see `configs/da_20250520.json`)::

    {"case": "<case label for the shared HRRR/GOES/MRMS single-panel files>",
     "ref_run": "<run name whose wrfout gives domain geometry>",
     "runs": [{"name": ..., "label": ..., "dir": ...}, ...],
     "diffs": [["a", "b"], ...],          # optional; a - b difference maps
     "side": ["left", "right"],           # optional; side-by-side movie pair
     "side_labels": {"name": "title"}}    # optional; titles for the side pair

A run whose directory holds no `wrfout_<domain>_*` file yet (e.g. a run that
has not been made) is dropped with a warning rather than failing, so one
config can list runs that come online later. Relative `dir`s resolve against
the working directory. Run names must be unique.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class RunSet:
    case: str
    ref_run: str
    runs: dict[str, Path]          # name -> run dir (only runs with output)
    labels: dict[str, str]         # name -> plot label
    diffs: list[tuple[str, str]] = field(default_factory=list)
    side: tuple[str, str] | None = None   # (left, right) for side-by-side movies
    side_labels: dict[str, str] = field(default_factory=dict)

    @property
    def ref_dir(self) -> Path:
        return self.runs[self.ref_run]


def load_runset(path: str | Path, domain: str = "d01") -> RunSet:
    cfg = json.loads(Path(path).read_text())
    runs, labels = {}, {}
    for r in cfg["runs"]:
        name, d = r["name"], Path(r["dir"])
        if name in runs or name in labels:
            raise ValueError(f"duplicate run name {name!r} in {path}")
        labels[name] = r.get("label", name)
        if d.is_dir() and any(d.glob(f"wrfout_{domain}_*")):
            runs[name] = d
        else:
            warnings.warn(f"run {name!r}: no wrfout_{domain}_* in {d}; skipping", stacklevel=2)
    ref = cfg["ref_run"]
    if ref not in runs:
        raise ValueError(f"ref_run {ref!r} has no output")
    diffs = [tuple(p) for p in cfg.get("diffs", []) if p[0] in runs and p[1] in runs]
    side = tuple(cfg["side"]) if cfg.get("side") and all(n in runs for n in cfg["side"]) else None
    return RunSet(cfg["case"], ref, runs, labels, diffs, side, cfg.get("side_labels", {}))

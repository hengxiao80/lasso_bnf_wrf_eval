"""JSON run-set config: point the batch plotters at run directories that
are not the built-in `satoshi_testruns/<case>/<scheme>` layout (e.g. the
`bnf_wrf_da` experiments).

Schema (see `configs/da_20250520.json`)::

    {"case": "<case label for the shared HRRR/GOES/MRMS single-panel files>",
     "ref_run": "<run name whose wrfout gives domain geometry>",
     "runs": [{"name": ..., "label": ..., "dir": ...}, ...],
     "pairs": [["a", "b"], ...]}          # optional; each pair gets a | b | a - b frames

The reference run (`ref_run`) also goes in the top row of the multi-run panels.
Pairs naming a run without output are dropped. """

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
    pairs: list[tuple[str, str]] = field(default_factory=list)  # (a, b): frames show a | b | a - b

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
    pairs = [tuple(p) for p in cfg.get("pairs", []) if all(n in runs for n in p)]
    return RunSet(cfg["case"], ref, runs, labels, pairs)

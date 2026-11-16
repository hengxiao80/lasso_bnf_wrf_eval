"""Reading the ensemble-mean reflectivity diagnostics of the `bnf_wrf_da`
WRF-DART runs.

The DA run's evaluation files (written in `bnf_wrf_da`, e.g.
`cases/20250520/da_theta_m0_n50/diag_conv_refl_f2/thompson_refl_mean/`) hold
one netCDF-4 file per analysis time, `thompson_refl_mean_d01_<YYYY-MM-DD_HH_MM_SS>.nc`,
with eight 2-D column maxima (dBZ) of ensemble-mean reflectivity plus `Times`,
`XLAT`, `XLONG` and a member wrfout's global attributes (so
`wrf.get_lambert_projection` works on them). Every field is the column maximum
of a 3-D ensemble mean, never a mean of the members' column maxima:

- `REFL_COLMAX_<stage>_DART`: DART's mean of `REFL_10CM` (`postassim_mean.nc`
  for the analysis, `input_mean.nc` for the forecast) -- the mean, in dBZ, of
  DART's regression-updated reflectivity, not computed from the hydrometeors.
- `REFL_COLMAX_<stage>_MEMDBZ`: each member's reflectivity recomputed from its
  hydrometeors (WRF's Thompson `calc_refl10cm`, run offline), averaged over the
  members in dBZ (members without echo at WRF's -35 dBZ floor).
- `REFL_COLMAX_<stage>_MEMZ`: the same, averaged in linear Z.
- `REFL_COLMAX_<stage>_STATE`: the reflectivity of the mean of the members'
  fields.

`<stage>` is `ANL` (the analysis, the state WRF starts each 15-min segment
from) or `FCST` (the 15-min forecast valid at the same time). The outermost row
and column are 0 dBZ in `DART` (as in WRF's own `REFL_10CM`, which WRF doesn't
compute there) and -35 dBZ in the recomputed fields.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np

from . import wrf as wrf_reader

STAGES = ("ANL", "FCST")
MEANS = ("DART", "MEMDBZ", "MEMZ", "STATE")


def find_da_refl_mean_file(mean_dir: str | Path, time: dt.datetime, domain: str = "d01") -> Path:
    """The `thompson_refl_mean_<domain>_<time>.nc` file in `mean_dir` for an
    analysis time; `FileNotFoundError` if there is none."""
    f = Path(mean_dir) / f"thompson_refl_mean_{domain}_{time:%Y-%m-%d_%H_%M_%S}.nc"
    if not f.exists():
        raise FileNotFoundError(f"No DA mean-reflectivity file for {time} at {f}")
    return f


def read_da_column_max_reflectivity(
    da_file: str | Path, mean: str, stage: str = "ANL"
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Column-max ensemble-mean reflectivity (dBZ) from a DA mean-reflectivity
    file: `REFL_COLMAX_<stage>_<mean>`, `mean` one of `MEANS`, `stage` one of
    `STAGES` (see the module docstring for what each is).

    Returns (lon, lat, values, valid_time), like the wrfout readers.
    """
    mean, stage = mean.upper(), stage.upper()
    if mean not in MEANS or stage not in STAGES:
        raise ValueError(f"mean must be one of {MEANS} and stage one of {STAGES}; got {mean!r}, {stage!r}")
    return wrf_reader.read_field(da_file, f"REFL_COLMAX_{stage}_{mean}")

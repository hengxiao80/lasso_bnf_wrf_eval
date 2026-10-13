"""Reading fields and map projection info from LASSO-BNF wrfout files."""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import cartopy.crs as ccrs
import netCDF4
import numpy as np


def get_lambert_projection(wrf_file: str | Path) -> ccrs.LambertConformal:
    """Cartopy Lambert Conformal projection matching a WRF domain's grid.

    Mirrors the convention wrf-python uses: the cone is defined by the
    domain's TRUELAT1/TRUELAT2 standard parallels, centered on STAND_LON /
    MOAD_CEN_LAT, on the spherical earth (radius 6,370,000 m) that WRF
    itself assumes rather than cartopy's default WGS84 ellipsoid.
    """
    with netCDF4.Dataset(wrf_file) as nc:
        globe = ccrs.Globe(ellipse=None, semimajor_axis=6370000, semiminor_axis=6370000)
        return ccrs.LambertConformal(
            central_longitude=float(nc.STAND_LON),
            central_latitude=float(nc.MOAD_CEN_LAT),
            standard_parallels=(float(nc.TRUELAT1), float(nc.TRUELAT2)),
            globe=globe,
        )


def read_cell_edges_xy(
    wrf_file: str | Path, proj: ccrs.Projection
) -> tuple[np.ndarray, np.ndarray]:
    """WRF's exact mass-grid cell edges, in `proj`'s projected x/y
    coordinates -- derived from the staggered U/V-point lat/lon arrays
    (`XLONG_U`/`XLAT_U`, `XLONG_V`/`XLAT_V`) rather than approximated by
    extrapolating outward from the mass-point (`XLONG`/`XLAT`) centers.

    WRF's Arakawa-C grid staggers U points in the west_east direction and
    V points in the south_north direction, each one point wider than the
    mass grid in that direction (confirmed against a real wrfout: mass
    (389, 549), U (389, 550), V (390, 549)) -- i.e. U/V points sit exactly
    on the west_east/south_north mass-cell edges. Both `point_lon`/
    `point_lat` and these U/V edges get projected through the *same*
    `proj` (matching WRF's own internal Lambert formula -- see
    `get_lambert_projection`) before any binning happens, so this is never
    comparing raw lat/lon against a rotated grid -- everything lives in one
    shared Lambert x/y plane.

    Reduces the full 2D projected U/V arrays to one representative row
    (for `x_edges`) and one representative column (for `y_edges`), on the
    assumption that they're separable -- checked empirically against a
    real wrfout: projected U-point x varies by at most ~3 m between the
    first and last row (out of a 2500 m cell), and projected V-point y by
    at most ~2 m between the first and last column, consistent with
    `XLONG_U`/`XLAT_V` being stored as float32 (~1 m precision at this
    magnitude) rather than any genuine non-uniformity -- i.e. WRF's grid
    really is regular in its own native projection, as expected by
    construction, and this reprojection preserves that to well under a
    MRMS pixel's size.

    Returns (x_edges, y_edges): 1D arrays of length (west_east+1) and
    (south_north+1) respectively, e.g. for use in
    `lightning.bin_points_to_wrf_grid`.
    """
    with netCDF4.Dataset(wrf_file) as nc:
        lon_u = np.ma.filled(nc.variables["XLONG_U"][0, ...], np.nan)
        lat_u = np.ma.filled(nc.variables["XLAT_U"][0, ...], np.nan)
        lon_v = np.ma.filled(nc.variables["XLONG_V"][0, ...], np.nan)
        lat_v = np.ma.filled(nc.variables["XLAT_V"][0, ...], np.nan)

    u_xyz = proj.transform_points(ccrs.PlateCarree(), lon_u[0, :], lat_u[0, :])
    v_xyz = proj.transform_points(ccrs.PlateCarree(), lon_v[:, 0], lat_v[:, 0])
    x_edges = u_xyz[:, 0]
    y_edges = v_xyz[:, 1]
    return x_edges, y_edges


def read_field(
    wrf_file: str | Path, varname: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Read a single-time 2D field plus its lat/lon grid from a wrfout file.

    Returns (lon, lat, values, valid_time); lon/lat/values all share shape
    (south_north, west_east).
    """
    with netCDF4.Dataset(wrf_file) as nc:
        values = np.ma.filled(nc.variables[varname][0, ...], np.nan)
        lon = np.ma.filled(nc.variables["XLONG"][0, ...], np.nan)
        lat = np.ma.filled(nc.variables["XLAT"][0, ...], np.nan)
        time_str = b"".join(nc.variables["Times"][0, :]).decode("utf-8")
        valid_time = dt.datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")
    return lon, lat, values, valid_time


def _get_wrf_getvar():
    """Import wrf-python's `getvar` lazily, working around two issues from
    it being last released in 2021 (see pyproject.toml for the packaging
    issues just getting it installed):

    - It needs `numpy.float_`, removed in NumPy 2.0.
    - Its `is_multi_file()` helper decides whether an input is a single
      NetCDF file vs. a sequence of them via `isinstance(wrfin, Iterable)`.
      Modern `netCDF4.Dataset` now implements `__iter__` (dict-like, over
      variable names) -- something wrf-python never anticipated -- so a
      plain open Dataset now gets misdetected as a multi-file sequence,
      and wrf-python then tries to `copy()` it (to get a "resettable"
      iterator over the sequence), which `Dataset` doesn't support,
      crashing with "Dataset is not picklable". Patched to still treat a
      plain Dataset as a single file; a real multi-file sequence (a list
      of Datasets, e.g.) is unaffected since it isn't a Dataset itself.
    """
    if not hasattr(np, "float_"):
        np.float_ = np.float64

    from wrf import util as wrf_util

    if not isinstance(wrf_util.is_multi_file, _PatchedIsMultiFile):
        wrf_util.is_multi_file = _PatchedIsMultiFile(wrf_util.is_multi_file)

    from wrf import getvar

    return getvar


class _PatchedIsMultiFile:
    def __init__(self, orig):
        self._orig = orig

    def __call__(self, wrfin):
        if isinstance(wrfin, netCDF4.Dataset):
            return False
        return self._orig(wrfin)


def read_cloud_top_temperature(
    wrf_file: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Cloud-top temperature (K), via `wrf-python`'s `ctt` diagnostic.

    `ctt` derives a physically-based cloud-top temperature from the model
    state (an optical-depth-threshold method over the hydrometeor
    profiles), designed specifically as the model-side analog of what a
    clean IR window channel -- e.g. GOES ABI channel 13, 10.3 micron --
    observes, for comparison against satellite brightness temperature (see
    `goes.read_brightness_temperature`). Falls back to clear-sky
    surface/skin temperature wherever `ctt` judges there's no cloud, same
    as it always does.
    """
    getvar = _get_wrf_getvar()

    with netCDF4.Dataset(wrf_file) as nc:
        ctt = getvar(nc, "ctt", units="K")
        lon = np.asarray(ctt.coords["XLONG"].values)
        lat = np.asarray(ctt.coords["XLAT"].values)
        values = np.asarray(ctt.values)
        time_str = b"".join(nc.variables["Times"][0, :]).decode("utf-8")
        valid_time = dt.datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")
    return lon, lat, values, valid_time


def brightness_temperature_from_olr(olr: np.ndarray) -> np.ndarray:
    """IR-window brightness temperature (K) from broadband TOA OLR (W/m^2).

    Uses the widely-used empirical fit

        OLR = sigma * Tf**4,   Tf = Tb * (a + b*Tb)

    with a = 1.228, b = -1.106e-3 K^-1 and sigma the Stefan-Boltzmann
    constant, inverted for Tb by solving the quadratic
    ``b*Tb**2 + a*Tb - Tf = 0``. The relation is from Yang and Slingo
    (2001, MWR 129, 784-801; coefficients originally Ohring et al. 1984)
    and is the same conversion PyFLEXTRKR applies in
    ``ftfunctions.olr_to_tb``.

    This is a cheap, no-forward-model alternative to CRTM-derived Tb
    (`bnf_wrf_eval.crtm.read_simulated_brightness_temperature`) and to the
    `ctt` diagnostic (`read_cloud_top_temperature`). Because OLR is a
    broadband flux integrated over the whole column rather than a narrow
    window-channel radiance, the result is an "effective" temperature that
    runs warm relative to a clean IR window (e.g. GOES ABI channel 13) and
    smooths out the very coldest convective cloud tops -- but it needs only
    a field already present in every wrfout / HRRR file.
    """
    sigma = 5.67e-8  # Stefan-Boltzmann constant, W m^-2 K^-4
    a = 1.228
    b = -1.106e-3  # K^-1
    tf = (np.asarray(olr, dtype=float) / sigma) ** 0.25
    return (-a + np.sqrt(a**2 + 4 * b * tf)) / (2 * b)


def read_brightness_temperature_from_olr(
    wrf_file: str | Path, varname: str = "OLR"
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """IR-window brightness temperature (K), converted from WRF's TOA OLR.

    Reads wrfout's `OLR` (top-of-atmosphere outgoing long-wave radiation,
    W/m^2) and passes it through `brightness_temperature_from_olr` -- see
    that function for the formula, its provenance, and the caveats about
    how a broadband-flux effective temperature differs from a true
    narrow-channel Tb. Provided as a fast counterpart to
    `read_cloud_top_temperature` and the CRTM path, comparable to HRRR's
    `hrrr.read_brightness_temperature_from_olr`.
    """
    lon, lat, olr, valid_time = read_field(wrf_file, varname)
    return lon, lat, brightness_temperature_from_olr(olr), valid_time


def read_dyn_lightning_flash_counts(
    wrf_file_t1: str | Path, wrf_file_t2: str | Path
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, tuple[dt.datetime, dt.datetime]]:
    """Per-grid-cell lightning flash counts between two wrfout times, from
    WRF's Dynamic Lightning Scheme (`dyn_lightning_option=1`; new in WRF
    v4.8.0; Lynn, B. H., Y. Yair, C. Price, G. Kelman, and A. Clark, 2012:
    Predicting cloud-to-ground and intracloud lightning in weather forecast
    models. Wea. Forecasting, 27, 1470-1486).

    `LPOS`/`LNEG`/`LNEU` ("Positive/Negative Cloud to Ground Lightning
    Density", "Intra-Cloud Lightning Density") are running totals of stroke
    counts accumulated since the start of the simulation -- confirmed
    against `phys/module_ltng_strokes.F`'s `flash` subroutine
    (`lpos = lpos + int(e_p/j_pos)`, a state variable never reset) and
    empirically (their domain sums only ever increase across a run's
    `wrfout` files) -- not per-output-interval counts. So the number of new
    flashes in (t1, t2] at each grid cell is the *difference* between the
    two files' values, computed here.

    Returns (lon, lat, cg_pos, cg_neg, ic, (t1, t2)); `cg_pos`/`cg_neg`/`ic`
    share the wrfout grid shape. Callers wanting total CG use
    `cg_pos + cg_neg`; total lightning (CG + intracloud) is
    `cg_pos + cg_neg + ic`.

    Raises KeyError if `wrf_file_t1` has no `LPOS`/`LNEG`/`LNEU` (i.e. the
    run wasn't configured with `dyn_lightning_option=1`).
    """
    with netCDF4.Dataset(wrf_file_t1) as nc1, netCDF4.Dataset(wrf_file_t2) as nc2:
        for varname in ("LPOS", "LNEG", "LNEU"):
            if varname not in nc1.variables:
                raise KeyError(
                    f"{wrf_file_t1} has no {varname} -- was this run configured with "
                    "dyn_lightning_option=1?"
                )

        def _diff(varname: str) -> np.ndarray:
            v1 = np.ma.filled(nc1.variables[varname][0, ...], np.nan)
            v2 = np.ma.filled(nc2.variables[varname][0, ...], np.nan)
            return v2 - v1

        cg_pos = _diff("LPOS")
        cg_neg = _diff("LNEG")
        ic = _diff("LNEU")
        lon = np.ma.filled(nc1.variables["XLONG"][0, ...], np.nan)
        lat = np.ma.filled(nc1.variables["XLAT"][0, ...], np.nan)
        t1 = dt.datetime.strptime(
            b"".join(nc1.variables["Times"][0, :]).decode("utf-8"), "%Y-%m-%d_%H:%M:%S"
        )
        t2 = dt.datetime.strptime(
            b"".join(nc2.variables["Times"][0, :]).decode("utf-8"), "%Y-%m-%d_%H:%M:%S"
        )
    return lon, lat, cg_pos, cg_neg, ic, (t1, t2)


def read_column_max_reflectivity(
    wrf_file: str | Path, varname: str = "REFL_10CM", levels_per_read: int = 10
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Column-max radar reflectivity (dBZ), computed from the 3D field.

    wrfout also carries a precomputed `REFD_MAX` diagnostic, but at least
    in the 20250917lassobnfwrfhrrr3 test run it is identically zero at
    every time (likely never accumulated/reset correctly for this run),
    so this takes the max over height of the instantaneous 3D
    `REFL_10CM` field instead, which does hold real values.

    The max is taken over `levels_per_read` vertical levels at a time, so
    peak memory stays bounded for very large grids (a 100-m d3
    `wrfout_d02` is 149 x 2290 x 2800, ~3.8 GB as float32 before masked-
    array copies). Same result as a single whole-column read.
    """
    with netCDF4.Dataset(wrf_file) as nc:
        var = nc.variables[varname]
        values = None
        for k in range(0, var.shape[1], levels_per_read):
            chunk = np.nanmax(np.ma.filled(var[0, k:k + levels_per_read], np.nan), axis=0)
            values = chunk if values is None else np.fmax(values, chunk)
        lon = np.ma.filled(nc.variables["XLONG"][0, ...], np.nan)
        lat = np.ma.filled(nc.variables["XLAT"][0, ...], np.nan)
        time_str = b"".join(nc.variables["Times"][0, :]).decode("utf-8")
        valid_time = dt.datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")
    return lon, lat, values, valid_time


def read_accumulated_precipitation(
    wrf_file: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Surface precipitation accumulated since simulation start (mm, i.e.
    kg m^-2 of liquid equivalent): grid-scale `RAINNC`. None of this
    project's runs use a cumulus scheme, so `RAINC` is not read. Rain rates
    come from differencing two files' totals."""
    return read_field(wrf_file, "RAINNC")


def read_hydrometeor_paths(
    wrf_file: str | Path,
    species: tuple[str, ...] = ("QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP"),
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], dt.datetime]:
    """Column-integrated path (kg m^-2) of each requested hydrometeor species.

    For species mixing ratio `q` (kg per kg of *dry* air, WRF's convention),
    ``path = sum_k q_k * m_k`` with ``m_k = (c1h*MUT + c2h) * (-dnw) / g``,
    layer k's dry-air mass per unit area from WRF's hybrid-coordinate mass
    (`MUT = MU + MUB`; `C1H`, `C2H`, `DNW` from the file -- for a pure-sigma
    run `C1H = 1`, `C2H = 0`, for a hybrid run they blend to a pressure
    coordinate aloft). `(c1h*MUT + c2h)*(-dnw)` is the layer's dry-air
    pressure thickness, the same term WRF's hydrostatic integration uses;
    dividing by g gives mass. Verified: layers sum to `MUT` (to ~0.1 Pa), and
    against an independent `rho*dz` integral (agrees to ~0.4-2.4%).

    Returns ``(lon, lat, {species_name: path}, valid_time)``; a species
    missing from the file raises `KeyError`.
    """
    g = 9.81
    with netCDF4.Dataset(wrf_file) as nc:
        mut = np.ma.filled(nc.variables["MU"][0], np.nan) + np.ma.filled(nc.variables["MUB"][0], np.nan)
        c1h = np.ma.filled(nc.variables["C1H"][0], np.nan)
        c2h = np.ma.filled(nc.variables["C2H"][0], np.nan)
        dnw = np.ma.filled(nc.variables["DNW"][0], np.nan)
        # (nz, ny, nx) dry-air mass per unit area of each layer
        mass = (c1h[:, None, None] * mut[None] + c2h[:, None, None]) * (-dnw)[:, None, None] / g
        paths = {
            name: np.sum(np.ma.filled(nc.variables[name][0], 0.0) * mass, axis=0)
            for name in species
        }
        lon = np.ma.filled(nc.variables["XLONG"][0, ...], np.nan)
        lat = np.ma.filled(nc.variables["XLAT"][0, ...], np.nan)
        time_str = b"".join(nc.variables["Times"][0, :]).decode("utf-8")
        valid_time = dt.datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")
    return lon, lat, paths, valid_time


def read_water_paths(
    wrf_file: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    """Total liquid (`QCLOUD + QRAIN`) and ice (`QICE + QSNOW + QGRAUP`) water
    paths, kg m^-2, as sums of `read_hydrometeor_paths`. Returns
    ``(lon, lat, lwp, iwp, valid_time)``."""
    lon, lat, p, t = read_hydrometeor_paths(wrf_file)
    return lon, lat, p["QCLOUD"] + p["QRAIN"], p["QICE"] + p["QSNOW"] + p["QGRAUP"], t

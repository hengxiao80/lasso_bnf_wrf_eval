"""Side-by-side WRF vs. HRRR vs. observations comparison plots."""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Colormap, LinearSegmentedColormap

from . import crtm as crtm_reader
from . import goes as goes_reader
from . import hrrr as hrrr_reader
from . import lightning as lightning_reader
from . import mrms as mrms_reader
from . import nalma as nalma_reader
from . import wrf as wrf_reader

DEFAULT_REFL_LEVELS = np.arange(5, 76, 5)
# Flash counts are sparse, non-negative integers (mostly 0), not a smooth
# field -- geometric-ish spacing so both a quiet cell (1-2 flashes) and an
# active storm core (dozens) get distinguishable colors.
DEFAULT_LIGHTNING_COUNT_LEVELS = np.array([1, 2, 3, 5, 10, 20, 50, 100])
# CG-only counts (`DEFAULT_LIGHTNING_COUNT_LEVELS`) never actually reach
# 100 in practice, but total lightning (CG + intracloud) does -- checked
# against a real active-storm hour: WRF's true max there was 280 with the
# 99.9th percentile of nonzero cells already at 148, so capping at 100
# would flatten the most-intense storm-core cells (exactly the ones worth
# distinguishing) into one "off-scale" bucket. One extra level pushes that
# ceiling out to 200.
DEFAULT_LIGHTNING_TOTAL_COUNT_LEVELS = np.array([1, 2, 3, 5, 10, 20, 50, 100, 200])
# NALMA's raw VHF *sources* (not flash-clustered -- see `nalma.py`'s module
# docstring) are a completely different scale from flash counts: checked
# against a real active-storm hour, ~3 million sources domain-wide vs.
# WRF's 58,441 flashes over the same window (~51 sources/flash, consistent
# with LMA literature), with a single WRF grid cell reaching 9540 sources.
# A wide geometric spread is needed to show both quiet fringes and the
# most active storm cores.
DEFAULT_LIGHTNING_SOURCE_COUNT_LEVELS = np.array([10, 30, 100, 300, 1000, 3000, 10000])
# 180-315 K to match the enhancement curve below: grayscale above the 240K
# convective-cloud-top threshold, a rainbow enhancement below it.
DEFAULT_CTT_LEVELS = np.arange(180, 316, 5)


def _tb_enhanced_colormap() -> LinearSegmentedColormap:
    """Grayscale-above/rainbow-below IR brightness-temperature colormap, of
    the kind commonly used to make cold convective cloud tops stand out
    (e.g. RAMMB/CIRA's "Rainbow IR" enhancement): black (warm) fading to
    white at 240 K, then a sharp jump into a cyan-to-magenta rainbow for
    the coldest (highest, most vigorous convective) cloud tops down to
    180 K. Matches a reference image the user supplied
    (`notebooks/tb-colorbar.png`).

    Calibrated against `DEFAULT_CTT_LEVELS` (180-315 K): matplotlib maps a
    colormap's [0, 1] domain evenly across whatever `levels` a caller
    passes to `contourf`, so the fixed physical anchor points below (240 K
    as the gray/rainbow boundary, etc.) only land where intended if
    `ctt_levels` keeps that same 180-315 K range. Passing a very different
    `ctt_levels` still works mechanically, just without that physical
    correspondence.
    """
    t_min, t_max = 180.0, 315.0

    def frac(t: float) -> float:
        return (t - t_min) / (t_max - t_min)

    stops = [
        (frac(180), (0.98, 0.90, 0.93)),  # pale pink (coldest, <180K arrow)
        (frac(190), (0.75, 0.30, 0.60)),  # magenta
        (frac(198), (0.45, 0.05, 0.10)),  # dark maroon
        (frac(206), (0.85, 0.10, 0.10)),  # red
        (frac(214), (0.95, 0.55, 0.05)),  # orange
        (frac(222), (0.95, 0.92, 0.15)),  # yellow
        (frac(230), (0.15, 0.70, 0.20)),  # green
        (frac(236), (0.10, 0.35, 0.85)),  # blue
        (frac(240) - 1e-6, (0.40, 0.90, 0.95)),  # cyan, right at the boundary
        (frac(240), (0.97, 0.97, 0.97)),  # sharp jump to near-white
        (frac(255), (0.80, 0.80, 0.80)),
        (frac(270), (0.55, 0.55, 0.55)),
        (frac(285), (0.35, 0.35, 0.35)),
        (frac(300), (0.15, 0.15, 0.15)),
        (frac(315), (0.02, 0.02, 0.02)),  # near-black (warmest)
    ]
    return LinearSegmentedColormap.from_list("tb_enhanced", stops)


DEFAULT_TB_CMAP = _tb_enhanced_colormap()

# DOE ARM Bankhead National Forest (BNF) site.
BNF_SITE_LAT = 34.342481
BNF_SITE_LON = -87.338177


def _domain_extent(lon: np.ndarray, lat: np.ndarray, pad_deg: float) -> tuple[float, float, float, float]:
    return (
        float(np.nanmin(lon)) - pad_deg,
        float(np.nanmax(lon)) + pad_deg,
        float(np.nanmin(lat)) - pad_deg,
        float(np.nanmax(lat)) + pad_deg,
    )


def _crop_to_extent(
    lon: np.ndarray, lat: np.ndarray, values: np.ndarray, extent: tuple[float, float, float, float], pad_cells: int = 30
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Slice a curvilinear grid down to the rows/cols overlapping `extent`.

    Keeps pcolormesh fast and avoids projecting the full (e.g. CONUS-wide
    HRRR) grid when only a small sub-domain is actually shown.

    `pad_cells` needs to be generous (not just 1-2 cells): pcolormesh
    mis-renders the outermost row/col of a rotated curvilinear mesh (shows
    as a spurious blank diagonal wedge) when that edge sits close to the
    visible map extent, so the crop must extend well past what's actually
    shown to push that artifact outside the visible area.
    """
    lon_min, lon_max, lat_min, lat_max = extent
    mask = (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)
    if not mask.any():
        return lon, lat, values
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    r0, r1 = max(rows.min() - pad_cells, 0), min(rows.max() + pad_cells + 1, lon.shape[0])
    c0, c1 = max(cols.min() - pad_cells, 0), min(cols.max() + pad_cells + 1, lon.shape[1])
    sl = (slice(r0, r1), slice(c0, c1))
    return lon[sl], lat[sl], values[sl]


def _domain_outline(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Trace the outer edge of a 2D curvilinear grid as one closed loop
    (its four sides, corner to corner), for drawing a domain boundary
    line."""
    lon_b = np.concatenate([lon[0, :], lon[:, -1], lon[-1, ::-1], lon[::-1, 0]])
    lat_b = np.concatenate([lat[0, :], lat[:, -1], lat[-1, ::-1], lat[::-1, 0]])
    return lon_b, lat_b


def _find_wrf_file(run_dir: Path, domain: str, time: dt.datetime) -> Path:
    wrf_file = run_dir / f"wrfout_{domain}_{time:%Y-%m-%d_%H_%M_%S}"
    if not wrf_file.exists():
        raise FileNotFoundError(f"No wrfout file for {time} at {wrf_file}")
    return wrf_file


def _find_hrrr_file(hrrr_base_dir: str | Path, time: dt.datetime) -> Path:
    hrrr_file = Path(hrrr_base_dir) / f"{time:%Y%m%d}" / f"hrrr.t{time:%H}z.wrfnatf00.grib2"
    if not hrrr_file.exists():
        raise FileNotFoundError(f"No HRRR native-level analysis file for {time} at {hrrr_file}")
    return hrrr_file


def _resolve_mrms_file(time: dt.datetime, mrms_dir: str | Path, auto_download: bool) -> Path:
    """Locate the MRMS file for `time`: from local filenames first (no S3
    call -- works on internet-less compute nodes), falling back to the
    S3-backed lookup+download only when `auto_download` is set."""
    mrms_dir = Path(mrms_dir)
    try:
        return mrms_reader.find_local_mrms_file(time, mrms_dir)
    except FileNotFoundError:
        if not auto_download:
            raise
    bucket, key, _ = mrms_reader.find_mrms_file(time)
    return mrms_reader.download_mrms_file(bucket, key, mrms_dir)


def _resolve_goes_file(
    time: dt.datetime, goes_dir: str | Path, channel: int, satellite: str, auto_download: bool
) -> Path:
    """GOES counterpart of `_resolve_mrms_file`."""
    goes_dir = Path(goes_dir)
    try:
        return goes_reader.find_local_goes_file(
            time, goes_dir, channel=channel, satellite=satellite
        )
    except FileNotFoundError:
        if not auto_download:
            raise
    bucket, key, _ = goes_reader.find_goes_file(time, channel=channel, satellite=satellite)
    return goes_reader.download_goes_file(bucket, key, goes_dir)


def _output_path(output_base_dir: str | Path, prefix: str, run_dir: Path, time: dt.datetime) -> Path:
    output_base_dir = Path(output_base_dir)
    output_base_dir.mkdir(parents=True, exist_ok=True)
    run_label = f"{run_dir.parent.name}_{run_dir.name}"
    return output_base_dir / f"{prefix}_{run_label}_{time:%Y%m%d_%H%MZ}.png"


def _plot_panel(
    ax,
    lon: np.ndarray,
    lat: np.ndarray,
    values: np.ndarray,
    levels: np.ndarray,
    cmap: str | Colormap,
    extent: tuple[float, float, float, float],
    title: str,
    extend: str = "both",
    domain_outline: tuple[np.ndarray, np.ndarray] | None = None,
    site: tuple[float, float] | None = (BNF_SITE_LON, BNF_SITE_LAT),
    plot_type: str = "contourf",
    dark_background: bool = False,
    lon_grid_step: float | None = None,
):
    fg_color = "white" if dark_background else "black"
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    if dark_background:
        ax.set_facecolor("black")
    # Pre-project lon/lat to the axes' own X/Y and call contourf/pcolormesh
    # with no `transform=` kwarg, rather than letting cartopy do it:
    # cartopy's contourf reprojection drops cells (leaves visible holes)
    # for a curvilinear grid that's rotated relative to the axes projection
    # -- e.g. HRRR's own Lambert grid plotted on WRF's Lambert axes -- even
    # though the same grid renders fine with pcolormesh.
    xyz = ax.projection.transform_points(ccrs.PlateCarree(), lon, lat)
    x, y = xyz[..., 0], xyz[..., 1]
    if plot_type == "contourf":
        mesh = ax.contourf(x, y, values, levels=levels, cmap=cmap, extend=extend)
    elif plot_type == "pcolormesh":
        # For sparse, mostly-NaN data (e.g. per-cell lightning flash
        # counts), `contourf` is the wrong tool: it fills the area *between*
        # grid points above a level, so an isolated non-zero cell entirely
        # surrounded by NaN neighbors has no contiguous region to fill and
        # renders as nothing -- real data silently vanishing, not a display
        # quirk. `pcolormesh` colors each cell on its own regardless of its
        # neighbors, which is what sparse per-cell data needs.
        from matplotlib.colors import BoundaryNorm

        # `.copy()` before `set_over`: `plt.get_cmap` returns matplotlib's
        # shared registered instance, so mutating it in place would leak
        # into every other plot using this same cmap name.
        cmap_obj = (plt.get_cmap(cmap) if isinstance(cmap, str) else cmap).copy()
        if dark_background and extend in ("max", "both"):
            # On a black background, the top of most sequential colormaps
            # (e.g. YlOrRd's dark red) is a genuinely dark color -- exactly
            # where the actual maxima land -- so it doesn't stand out
            # against black the way the low end does. Force the
            # off-the-top-of-the-scale bin to pure white instead.
            cmap_obj.set_over("white")
        norm = BoundaryNorm(levels, cmap_obj.N, extend=extend)
        mesh = ax.pcolormesh(x, y, values, cmap=cmap_obj, norm=norm, shading="nearest")
    else:
        raise ValueError(f"Unknown plot_type: {plot_type!r}")
    ax.coastlines(resolution="50m", linewidth=0.8, color=fg_color)
    # Pin an explicit resolution on every feature (BORDERS defaults to
    # cartopy's AdaptiveScaler, which picks 110m/50m/10m based on each
    # plot's map extent). Left un-pinned, a domain with a different extent
    # than previously plotted can make cartopy pick a resolution that
    # hasn't been downloaded yet and try to fetch it from
    # naturalearthdata.com mid-plot -- if that network call stalls, it
    # hangs the whole kernel in a way SIGINT can't reliably interrupt.
    # Pinning means every plot uses the same, already-cached files.
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.8, edgecolor=fg_color)
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.5, edgecolor=fg_color)
    # Explicit dict form (rather than the top_labels/right_labels booleans)
    # so cartopy doesn't fall back to its per-label nearest-edge geometry
    # guess, which for a Lambert projection can misplace some longitude
    # labels on the side panels.
    gl = ax.gridlines(
        draw_labels={"bottom": "x", "left": "y"}, linestyle="--",
        color="white" if dark_background else "gray", alpha=0.5,
    )
    gl.x_inline = False
    gl.y_inline = False
    # Keep longitude labels horizontal and pushed out below the axis,
    # rather than cartopy's default of rotating them to follow the
    # (slightly slanted, for a Lambert projection) gridline -- which
    # otherwise lands them inside the plot near the bottom edge.
    gl.rotate_labels = False
    gl.xpadding = 8
    if lon_grid_step is not None:
        # cartopy's automatic choice drops to 0.5-degree longitude lines on a
        # small (~3-degree-wide) map, whose bottom labels then collide.
        import matplotlib.ticker as mticker

        gl.xlocator = mticker.MultipleLocator(lon_grid_step)
    if dark_background:
        gl.xlabel_style = {"color": "white"}
        gl.ylabel_style = {"color": "white"}
        # The axes' own frame/border (cartopy's "geo" spine) defaults to
        # black -- invisible once both the figure and axes backgrounds are
        # also black, making it impossible to tell where one panel ends
        # and the next begins. Newer cartopy exposes it as a normal
        # matplotlib spine; fall back to the older `outline_patch` API.
        try:
            ax.spines["geo"].set_edgecolor("white")
        except (KeyError, AttributeError):
            ax.outline_patch.set_edgecolor("white")
    if domain_outline is not None:
        outline_lon, outline_lat = domain_outline
        ax.plot(
            outline_lon, outline_lat, transform=ccrs.PlateCarree(),
            color="red", linewidth=0.8, zorder=10,
        )
    if site is not None:
        site_lon, site_lat = site
        ax.plot(
            site_lon, site_lat, transform=ccrs.PlateCarree(),
            marker="*", markersize=10, color="red",
            markeredgecolor=fg_color, markeredgewidth=0.5, linestyle="none", zorder=11,
        )
    ax.set_title(title, fontsize=10, color=fg_color)
    return mesh


def _plot_row(
    fig,
    axes_row,
    panels: list[tuple[np.ndarray, np.ndarray, np.ndarray, str]],
    levels: np.ndarray,
    cmap: str | Colormap,
    extent: tuple[float, float, float, float],
    domain_outline: tuple[np.ndarray, np.ndarray],
    extend: str,
    colorbar_label: str,
    plot_type: str = "contourf",
    dark_background: bool = False,
    lon_grid_step: float | None = None,
):
    """Draw one row of side-by-side comparison panels (e.g. WRF / HRRR /
    obs for the same variable, all on the same levels/cmap), sharing one
    colorbar. `panels` is a list of (lon, lat, values, title) tuples, one
    per axes in `axes_row`. `plot_type` picks the fill method -- see
    `_plot_panel`; use `"pcolormesh"` for sparse, mostly-NaN data.
    `dark_background` switches map features/labels/colorbar text to white,
    for use with a black figure/axes background (set by the caller)."""
    mesh = None
    for ax, (lon, lat, values, title) in zip(axes_row, panels):
        mesh = _plot_panel(
            ax, lon, lat, values, levels, cmap, extent, title,
            extend=extend, domain_outline=domain_outline, plot_type=plot_type,
            dark_background=dark_background, lon_grid_step=lon_grid_step,
        )
    cbar = fig.colorbar(mesh, ax=axes_row[:], label=colorbar_label, shrink=0.5, pad=0.02)
    if dark_background:
        cbar.set_label(colorbar_label, color="white")
        cbar.ax.yaxis.set_tick_params(color="white")
        plt.setp(plt.getp(cbar.ax, "yticklabels"), color="white")


def plot_refl_comparison(
    wrf_file: str | Path,
    hrrr_file: str | Path,
    mrms_file: str | Path,
    refl_levels: np.ndarray = DEFAULT_REFL_LEVELS,
    refl_cmap: str = "turbo",
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (16, 5.5),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """Compare column-max radar reflectivity across WRF, HRRR, and MRMS
    (observations), at the single time each of the given files holds.

    Unlike the OLR-vs-brightness-temperature situation (see
    `plot_tb_comparison`), MRMS composite reflectivity is a direct,
    apples-to-apples match for WRF's simulated reflectivity: same physical
    quantity (dBZ), same basis as HRRR's `refc` -- see `mrms.py`.

    All three panels are drawn on the same Lambert Conformal projection
    (taken from the WRF domain) and the same map extent -- derived from
    the WRF domain's footprint, padded by `domain_pad_deg` -- so they're
    directly comparable. Because the WRF domain is much smaller than
    HRRR's/MRMS's, the padded extent extends past the WRF grid on all
    sides, so the WRF panel shows empty corners; that's expected. The WRF
    domain's actual footprint is outlined in thin red, and the BNF site is
    marked with a red star, on every panel.

    `refl_levels` sets the contourf contour intervals (bin edges) and can
    be overridden per call.

    Parameters
    ----------
    wrf_file : path to a wrfout_d0X_* file (single time per file).
    hrrr_file : path to a HRRR native-level analysis GRIB2 file
        (hrrr.tHHz.wrfnatf00.grib2).
    mrms_file : path to a downloaded MRMS composite reflectivity GRIB2.gz
        file (see `mrms.py`) for a scan time close to `wrf_file`'s.
    suptitle : if given, drawn as a big figure-level title above all three
        panels; see `plot_run_refl_comparison` for a wrapper that fills
        this in automatically.
    out_file : if given, the figure is saved there; otherwise it's left
        for the caller to show/save.

    Returns
    -------
    matplotlib.figure.Figure
    """
    wrf_lon, wrf_lat, wrf_refl, wrf_time = wrf_reader.read_column_max_reflectivity(wrf_file)
    proj = wrf_reader.get_lambert_projection(wrf_file)

    hrrr_lon, hrrr_lat, hrrr_refl, hrrr_time = hrrr_reader.read_composite_reflectivity(hrrr_file)
    mrms_lon, mrms_lat, mrms_refl, mrms_time = mrms_reader.read_composite_reflectivity(mrms_file)

    extent = _domain_extent(wrf_lon, wrf_lat, domain_pad_deg)
    hrrr_lon_c, hrrr_lat_c, hrrr_refl_c = _crop_to_extent(hrrr_lon, hrrr_lat, hrrr_refl, extent)
    mrms_lon_c, mrms_lat_c, mrms_refl_c = _crop_to_extent(mrms_lon, mrms_lat, mrms_refl, extent)
    domain_outline = _domain_outline(wrf_lon, wrf_lat)

    fig, axes = plt.subplots(
        1, 3, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)

    _plot_row(
        fig, axes,
        [
            (wrf_lon, wrf_lat, wrf_refl, "WRF max reflectivity (dBZ)"),
            (hrrr_lon_c, hrrr_lat_c, hrrr_refl_c, "HRRR composite reflectivity (dBZ)"),
            (
                mrms_lon_c, mrms_lat_c, mrms_refl_c,
                f"MRMS composite reflectivity (dBZ)\n{mrms_time:%Y-%m-%d %H:%M} UTC",
            ),
        ],
        refl_levels, refl_cmap, extent, domain_outline, extend="max", colorbar_label="dBZ",
    )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight="bold")

    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")

    return fig


def plot_run_refl_comparison(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    mrms_dir: str | Path = "mrms_data",
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_mrms: bool = False,
    **kwargs,
):
    """Wrapper around `plot_refl_comparison` that takes a time plus a WRF
    run's name and directory, locates the matching wrfout/HRRR/MRMS files
    itself, and labels the figure with a big suptitle showing `run_name`
    and `time`.

    Parameters
    ----------
    time : the output/analysis time to plot, e.g.
        `datetime.datetime(2025, 9, 16, 1)`, or an ISO string like
        `"2025-09-16 01:00"`.
    run_name : label for the run, used in the figure's suptitle (e.g.
        `"2.5-km (D1) Thompson"`).
    run_dir : directory holding that run's `wrfout_*` files (e.g.
        `satoshi_testruns/20250917lassobnfwrfhrrr3/rund1`).
    hrrr_base_dir : directory holding HRRR native-level GRIB2 files, laid
        out as `<hrrr_base_dir>/<YYYYMMDD>/hrrr.tHHz.wrfnatf00.grib2`.
        Defaults to `satoshi_forcing_data/hrrr/hrrrnat_data` relative to
        the repo root.
    mrms_dir : directory holding downloaded MRMS composite reflectivity
        files (as saved by `mrms.download_mrms_file`). Defaults to
        `mrms_data` relative to wherever this is run from.
    domain : WRF domain to plot, e.g. "d01" (default) or "d02".
    output_base_dir : if given (and `out_file` isn't passed explicitly via
        `**kwargs`), the figure is auto-saved to
        `<output_base_dir>/wrf_hrrr_mrms_refl_comparison_<run_dir's last
        two path components>_<time>.png`, creating the directory if
        needed -- using the directory names rather than `run_name`, since
        `run_name` is meant as a free-form label and isn't filename-safe.
    auto_download_mrms : if the matching MRMS file isn't already in
        `mrms_dir`, fetch it via `mrms.download_mrms_file` instead of
        raising `FileNotFoundError`. False by default, so plotting never
        triggers a network download unless you explicitly ask for it.
    **kwargs : forwarded to `plot_refl_comparison` (e.g. `refl_levels`,
        `out_file`).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    wrf_file = _find_wrf_file(run_dir, domain, time)
    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)

    mrms_file = _resolve_mrms_file(time, mrms_dir, auto_download_mrms)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(
            output_base_dir, "wrf_hrrr_mrms_refl_comparison", run_dir, time
        )

    suptitle = kwargs.pop("suptitle", f"{run_name} -- {time:%Y-%m-%d %H:%M} UTC")
    return plot_refl_comparison(wrf_file, hrrr_file, mrms_file, suptitle=suptitle, **kwargs)


def plot_tb_comparison(
    wrf_file: str | Path,
    hrrr_file: str | Path,
    goes_file: str | Path,
    ctt_levels: np.ndarray = DEFAULT_CTT_LEVELS,
    cmap: str | Colormap = DEFAULT_TB_CMAP,
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (16, 5.5),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
    wrf_tb_stride: int = 1,
    wrf_tb_kwargs: dict | None = None,
):
    """Compare brightness temperature across WRF, HRRR, and GOES
    (observations), at the single time each of the given files holds.

    OLR and satellite brightness temperature aren't directly comparable
    (OLR is a broadband, all-wavelength flux; Tb is a narrowband radiance
    at one wavelength converted via the Planck function), so this compares
    three things that *are* the same physical quantity throughout, all via
    the same forward-model approach (CRTM):

    - WRF: `crtm.read_simulated_brightness_temperature` -- this project's
      own CRTM forward-model calculation run directly on WRF's state, the
      same *approach* UPP uses for HRRR's `SBT114` (see below), applied
      straight to our own model output rather than relying on
      `wrf.read_cloud_top_temperature`'s simpler `ctt` diagnostic (WRF's
      own purpose-built approximation, kept in `wrf.py` for reference/
      comparison but no longer used here).
    - HRRR: `hrrr.read_simulated_brightness_temperature` -- HRRR/UPP's own
      CRTM-derived simulated brightness temperature (the same forward-model
      approach NCEP uses operationally to compare models against
      satellite).
    - GOES: `goes.read_brightness_temperature` -- the actual observed
      channel 13 (10.3 micron, the clean IR window channel) brightness
      temperature.

    Together this is about as close to an apples-to-apples three-way
    comparison as is practical: all three are either observed Tb or a
    CRTM-derived simulated Tb for the same channel.

    `cmap` (default `DEFAULT_TB_CMAP`, see `_tb_enhanced_colormap`) is a
    grayscale-above/rainbow-below IR enhancement: black (warm) fading to
    white at 240 K following the conventional satellite IR display
    convention (cold renders bright, warm renders dark), then a rainbow
    for the coldest (highest, most vigorous convective) cloud tops below
    that.

    Parameters otherwise mirror `plot_refl_comparison` -- see there for
    `domain_pad_deg`/`out_file` behavior.

    Parameters
    ----------
    wrf_file : path to a wrfout_d0X_* file (single time per file).
    hrrr_file : path to a HRRR native-level analysis GRIB2 file
        (hrrr.tHHz.wrfnatf00.grib2).
    goes_file : path to a GOES ABI CMIP NetCDF file (see `goes.py`) for a
        scan time close to `wrf_file`'s.
    suptitle : if given, drawn as a big figure-level title above all three
        panels; see `plot_run_tb_comparison` for a wrapper that fills this
        in automatically.
    wrf_tb_stride : forwarded to `crtm.read_simulated_brightness_temperature`
        as `stride` (default 1: full WRF resolution, no spatial
        subsampling). CRTM runs one profile per WRF column with no
        batching, so a naive whole-domain call can be too large for a
        login-node's memory budget -- `read_simulated_brightness_temperature`
        handles that itself by chunking the domain into memory-bounded
        row-bands and stitching the results back together (numerically
        exact -- see its docstring), and caches the result to disk (see
        `wrf_tb_kwargs`'s `cache_dir`/`use_cache`/`recompute`) since a
        full-resolution run still takes a couple of minutes the first time.
        Raise `wrf_tb_stride` for a quicker, coarser WRF panel (e.g. for a
        fast preview) if that's ever preferable to full resolution.
    wrf_tb_kwargs : extra keyword arguments forwarded to
        `crtm.read_simulated_brightness_temperature` (e.g. `sensor_id`,
        `sat_lon`/`sat_height` from a real GOES file's own projection
        metadata, `coefficient_path`, `max_profiles_per_chunk`, `cache_dir`,
        `use_cache`, `recompute`). `stride` goes through `wrf_tb_stride`
        above, not here.

    Returns
    -------
    matplotlib.figure.Figure
    """
    domain_lon, domain_lat, _, _ = wrf_reader.read_field(wrf_file, "HGT")
    proj = wrf_reader.get_lambert_projection(wrf_file)

    wrf_lon, wrf_lat, wrf_tb, wrf_time = crtm_reader.read_simulated_brightness_temperature(
        wrf_file, stride=wrf_tb_stride, **(wrf_tb_kwargs or {})
    )
    hrrr_lon, hrrr_lat, hrrr_tb, hrrr_time = hrrr_reader.read_simulated_brightness_temperature(
        hrrr_file
    )
    goes_lon, goes_lat, goes_tb, goes_time = goes_reader.read_brightness_temperature(goes_file)

    extent = _domain_extent(domain_lon, domain_lat, domain_pad_deg)
    hrrr_lon_c, hrrr_lat_c, hrrr_tb_c = _crop_to_extent(hrrr_lon, hrrr_lat, hrrr_tb, extent)
    goes_lon_c, goes_lat_c, goes_tb_c = _crop_to_extent(goes_lon, goes_lat, goes_tb, extent)
    domain_outline = _domain_outline(domain_lon, domain_lat)

    fig, axes = plt.subplots(
        1, 3, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)

    _plot_row(
        fig, axes,
        [
            (wrf_lon, wrf_lat, wrf_tb, "WRF simulated brightness temp. (CRTM) (K)"),
            (hrrr_lon_c, hrrr_lat_c, hrrr_tb_c, "HRRR simulated brightness temp. (K)"),
            (
                goes_lon_c, goes_lat_c, goes_tb_c,
                f"GOES brightness temperature (K)\n{goes_time:%Y-%m-%d %H:%M} UTC scan",
            ),
        ],
        ctt_levels, cmap, extent, domain_outline, extend="both", colorbar_label="K",
    )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight="bold")

    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")

    return fig


def plot_run_tb_comparison(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    goes_dir: str | Path = "goes_data",
    channel: int = 13,
    satellite: str = goes_reader.DEFAULT_SATELLITE,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_goes: bool = False,
    **kwargs,
):
    """Wrapper around `plot_tb_comparison` that takes a time plus a WRF
    run's name and directory, locates the matching wrfout/HRRR/GOES files
    itself, and labels the figure with a big suptitle showing `run_name`
    and `time`.

    Parameters
    ----------
    time : the output/analysis time to plot; see `plot_run_refl_comparison`
        for accepted forms.
    run_name : label for the run, used in the figure's suptitle.
    run_dir : directory holding that run's `wrfout_*` files.
    hrrr_base_dir : directory holding HRRR native-level GRIB2 files.
        Defaults to `satoshi_forcing_data/hrrr/hrrrnat_data` relative to
        the repo root.
    goes_dir : directory holding downloaded GOES CMIP files (as saved by
        `goes.download_goes_file`, e.g. named
        `OR_ABI-L2-CMIPC-M6C13_G19_s...nc`). Defaults to `goes_data`
        relative to wherever this is run from.
    channel : GOES ABI channel to compare against (13, the clean IR window
        channel, by default; see `goes.find_goes_file`).
    satellite : which GOES satellite/bucket to look in (see
        `goes.DEFAULT_SATELLITE` for the GOES-16 -> GOES-19 operational
        switch this defaults around).
    domain : WRF domain to plot, e.g. "d01" (default) or "d02".
    output_base_dir : if given (and `out_file` isn't passed explicitly via
        `**kwargs`), the figure is auto-saved to
        `<output_base_dir>/wrf_hrrr_goes_tb_comparison_<run_dir's last two
        path components>_<time>.png`, creating the directory if needed --
        same convention as `plot_run_refl_comparison`'s `output_base_dir`.
    auto_download_goes : if the matching GOES file isn't already in
        `goes_dir`, fetch it via `goes.download_goes_file` instead of
        raising `FileNotFoundError`. False by default, so plotting never
        triggers a network download unless you explicitly ask for it.
    **kwargs : forwarded to `plot_tb_comparison` (e.g. `ctt_levels`,
        `out_file`, `wrf_tb_stride`, `wrf_tb_kwargs`).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    wrf_file = _find_wrf_file(run_dir, domain, time)
    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)

    goes_file = _resolve_goes_file(time, goes_dir, channel, satellite, auto_download_goes)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(
            output_base_dir, "wrf_hrrr_goes_tb_comparison", run_dir, time
        )

    suptitle = kwargs.pop("suptitle", f"{run_name} -- {time:%Y-%m-%d %H:%M} UTC")
    return plot_tb_comparison(wrf_file, hrrr_file, goes_file, suptitle=suptitle, **kwargs)


def plot_tb_comparison_4panel(
    wrf_file: str | Path,
    hrrr_file: str | Path,
    goes_file: str | Path,
    ctt_levels: np.ndarray = DEFAULT_CTT_LEVELS,
    cmap: str | Colormap = DEFAULT_TB_CMAP,
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (10.5, 7.5),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
    wrf_tb_stride: int = 1,
    wrf_tb_kwargs: dict | None = None,
    wrf_panel: str = "ctt",
):
    """Like `plot_tb_comparison`, but with a second WRF-side panel next to
    the CRTM-derived one, laid out as a 2x2 grid: <extra WRF panel> | WRF
    simulated Tb (CRTM) on top, HRRR simulated Tb | GOES observed Tb below.

    `wrf_panel` selects what the extra top-left panel shows:

    - ``"ctt"`` (default): WRF's `wrf.read_cloud_top_temperature` (`ctt`)
      diagnostic -- WRF's own purpose-built optical-depth cloud-top
      temperature.
    - ``"olr"``: WRF's quick OLR-derived brightness temperature
      (`wrf.read_brightness_temperature_from_olr`) -- a one-line algebraic
      conversion of the broadband TOA `OLR` field, the cheap
      no-forward-model alternative.

    Either way the panel is a direct side-by-side look at how that simpler
    WRF-side method compares with the full CRTM forward model.

    See `plot_tb_comparison` for the rationale behind the CRTM-based panels
    and for what all the other parameters do; `wrf_tb_stride`/
    `wrf_tb_kwargs` apply only to the CRTM panel (the extra panel always
    runs at WRF's native resolution -- both `ctt` and the OLR fit are
    lightweight, not per-column CRTM forward-model calls, so neither needs
    the same memory-driven subsampling).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if wrf_panel == "ctt":
        wrf_x_lon, wrf_x_lat, wrf_x_val, _ = wrf_reader.read_cloud_top_temperature(wrf_file)
        wrf_x_title = "WRF cloud-top temperature (ctt) (K)"
    elif wrf_panel == "olr":
        wrf_x_lon, wrf_x_lat, wrf_x_val, _ = wrf_reader.read_brightness_temperature_from_olr(
            wrf_file
        )
        wrf_x_title = "WRF simple brightness temp. (OLR fit) (K)"
    else:
        raise ValueError(f"wrf_panel must be 'ctt' or 'olr', got {wrf_panel!r}")
    proj = wrf_reader.get_lambert_projection(wrf_file)

    wrf_lon, wrf_lat, wrf_tb, wrf_time = crtm_reader.read_simulated_brightness_temperature(
        wrf_file, stride=wrf_tb_stride, **(wrf_tb_kwargs or {})
    )
    hrrr_lon, hrrr_lat, hrrr_tb, hrrr_time = hrrr_reader.read_simulated_brightness_temperature(
        hrrr_file
    )
    goes_lon, goes_lat, goes_tb, goes_time = goes_reader.read_brightness_temperature(goes_file)

    extent = _domain_extent(wrf_x_lon, wrf_x_lat, domain_pad_deg)
    hrrr_lon_c, hrrr_lat_c, hrrr_tb_c = _crop_to_extent(hrrr_lon, hrrr_lat, hrrr_tb, extent)
    goes_lon_c, goes_lat_c, goes_tb_c = _crop_to_extent(goes_lon, goes_lat, goes_tb, extent)
    domain_outline = _domain_outline(wrf_x_lon, wrf_x_lat)

    fig, axes = plt.subplots(
        2, 2, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)

    _plot_row(
        fig, axes.flatten(),
        [
            (wrf_x_lon, wrf_x_lat, wrf_x_val, wrf_x_title),
            (wrf_lon, wrf_lat, wrf_tb, "WRF simulated brightness temp. (CRTM) (K)"),
            (hrrr_lon_c, hrrr_lat_c, hrrr_tb_c, "HRRR simulated brightness temp. (K)"),
            (
                goes_lon_c, goes_lat_c, goes_tb_c,
                f"GOES brightness temperature (K)\n{goes_time:%Y-%m-%d %H:%M} UTC scan",
            ),
        ],
        ctt_levels, cmap, extent, domain_outline, extend="both", colorbar_label="K",
    )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight="bold")

    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")

    return fig


def plot_run_tb_comparison_4panel(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    goes_dir: str | Path = "goes_data",
    channel: int = 13,
    satellite: str = goes_reader.DEFAULT_SATELLITE,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_goes: bool = False,
    wrf_panel: str = "ctt",
    **kwargs,
):
    """Wrapper around `plot_tb_comparison_4panel`, otherwise identical to
    `plot_run_tb_comparison` -- see there for what every parameter does, and
    `plot_tb_comparison_4panel` for what `wrf_panel` (`"ctt"` | `"olr"`)
    selects for the extra top-left panel.

    Saved output filenames use the `wrf_ctt_crtm_hrrr_goes_tb_comparison`
    prefix for `wrf_panel="ctt"` and `wrf_simple_tb_crtm_hrrr_goes_tb_comparison`
    for `wrf_panel="olr"` (vs. `plot_run_tb_comparison`'s
    `wrf_hrrr_goes_tb_comparison`), so the three don't overwrite each other
    when run for the same run/time.

    **kwargs : forwarded to `plot_tb_comparison_4panel` (e.g. `ctt_levels`,
        `out_file`, `wrf_tb_stride`, `wrf_tb_kwargs`).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    wrf_file = _find_wrf_file(run_dir, domain, time)
    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)

    goes_file = _resolve_goes_file(time, goes_dir, channel, satellite, auto_download_goes)

    _prefix = {
        "ctt": "wrf_ctt_crtm_hrrr_goes_tb_comparison",
        "olr": "wrf_simple_tb_crtm_hrrr_goes_tb_comparison",
    }[wrf_panel]
    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(output_base_dir, _prefix, run_dir, time)

    suptitle = kwargs.pop("suptitle", f"{run_name} -- {time:%Y-%m-%d %H:%M} UTC")
    return plot_tb_comparison_4panel(
        wrf_file, hrrr_file, goes_file, suptitle=suptitle, wrf_panel=wrf_panel, **kwargs
    )


# ---------------------------------------------------------------------------
# d2 (500-m ndown) 4-panel comparisons
#
# One frame per output time, laid out as a 2x2 grid:
#
#     d2 (500 m, WRF)   |   d1 (2.5 km, WRF)
#     HRRR              |   independent observation
#
# All four panels are drawn on the d2 domain's Lambert projection and the
# same map extent -- the d2 footprint padded by `domain_pad_deg` (so the d2
# panel shows a small margin around its own grid, like the WRF panel in the
# d1 comparisons) -- with d1 / HRRR / obs cropped to that extent so every
# panel covers the same area as the d2 panel. The d2 domain's footprint is
# outlined in red and the BNF site marked with a red star on every panel
# (both handled by `_plot_row` -> `_plot_panel`).
#
# Brightness temperature uses WRF's / HRRR's cheap OLR-derived Tb
# (`*.read_brightness_temperature_from_olr`, the Yang & Slingo 2001
# `OLR = sigma*Tf**4` fit) for all three model panels, so no CRTM /
# `crtm_cache/` is involved; reflectivity uses column-max `REFL_10CM` (WRF)
# / `refc` (HRRR) / `MergedReflectivityQCComposite` (MRMS), the same three
# quantities `plot_refl_comparison` compares.
# ---------------------------------------------------------------------------

# d2 ndown run subdir -> its d1 (2.5-km parent) counterpart in the same case
# directory. The d1 panel reads this run's `wrfout_d01_*` (the 2.5-km outer
# domain, which fully contains d2 and shares its exact Lambert projection).
_D2_TO_D1_RUN = {
    "rund2": "rund1",
    "rund2-dynlit": "rund1-dynlit",
}


def _nested_four_panel(
    fig,
    axes,
    inner_panel,
    other_panels,
    levels,
    cmap,
    domain_pad_deg,
    extend,
    colorbar_label,
    lon_grid_step=None,
):
    """Shared body of the d2 and d3 `plot_d*_comparison_4panel` plotters.

    `inner_panel` is the (lon, lat, values, title) of the finest domain
    shown (d2 for the d2 plots, d3 for the d3 plots): its footprint, padded
    by `domain_pad_deg`, sets the map extent and is outlined in red.
    `other_panels` are the remaining three (lon, lat, values, title)
    panels, cropped to that extent. Draws the 2x2 grid (inner panel top-
    left, then the others in order) with one shared colorbar and the BNF
    marker on every panel."""
    inner_lon, inner_lat, _, _ = inner_panel
    extent = _domain_extent(inner_lon, inner_lat, domain_pad_deg)
    domain_outline = _domain_outline(inner_lon, inner_lat)

    cropped = []
    for lon, lat, values, title in other_panels:
        lon, lat, values = _crop_to_extent(lon, lat, values, extent)
        cropped.append((lon, lat, values, title))

    _plot_row(
        fig, axes.flatten(), [inner_panel, *cropped],
        levels, cmap, extent, domain_outline, extend=extend, colorbar_label=colorbar_label,
        lon_grid_step=lon_grid_step,
    )


def plot_d2_tb_comparison_4panel(
    d2_file: str | Path,
    d1_file: str | Path,
    hrrr_file: str | Path,
    goes_file: str | Path,
    ctt_levels: np.ndarray = DEFAULT_CTT_LEVELS,
    cmap: str | Colormap = DEFAULT_TB_CMAP,
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (12.6, 7.9),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """2x2 OLR-derived brightness-temperature comparison for a 500-m ndown
    run: d2 (500 m WRF) | d1 (2.5 km WRF) on top, HRRR | GOES ABI ch.13
    (observed) below.

    The three model panels all use the cheap OLR-fit Tb
    (`wrf.read_brightness_temperature_from_olr` /
    `hrrr.read_brightness_temperature_from_olr`) -- no CRTM. GOES is the
    real observed channel-13 (10.3 micron) brightness temperature.

    See the module comment above "d2 (500-m ndown) 4-panel comparisons" for
    the projection / extent / outline / BNF-marker conventions; parameters
    otherwise mirror `plot_tb_comparison_4panel`.
    """
    d2_lon, d2_lat, d2_tb, _ = wrf_reader.read_brightness_temperature_from_olr(d2_file)
    d1_lon, d1_lat, d1_tb, _ = wrf_reader.read_brightness_temperature_from_olr(d1_file)
    hrrr_lon, hrrr_lat, hrrr_tb, _ = hrrr_reader.read_brightness_temperature_from_olr(hrrr_file)
    goes_lon, goes_lat, goes_tb, goes_time = goes_reader.read_brightness_temperature(goes_file)
    proj = wrf_reader.get_lambert_projection(d2_file)

    fig, axes = plt.subplots(
        2, 2, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)

    _nested_four_panel(
        fig, axes,
        (d2_lon, d2_lat, d2_tb, "d2 (500 m) WRF OLR brightness temp. (K)"),
        [
            (d1_lon, d1_lat, d1_tb, "d1 (2.5 km) WRF OLR brightness temp. (K)"),
            (hrrr_lon, hrrr_lat, hrrr_tb, "HRRR OLR brightness temp. (K)"),
            (
                goes_lon, goes_lat, goes_tb,
                f"GOES ABI ch.13 brightness temp. (K)\n{goes_time:%Y-%m-%d %H:%M} UTC scan",
            ),
        ],
        ctt_levels, cmap, domain_pad_deg, "both", "K",
    )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight="bold")
    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
    return fig


def plot_d2_refl_comparison_4panel(
    d2_file: str | Path,
    d1_file: str | Path,
    hrrr_file: str | Path,
    mrms_file: str | Path,
    refl_levels: np.ndarray = DEFAULT_REFL_LEVELS,
    refl_cmap: str = "turbo",
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (12.6, 7.9),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """2x2 column-max radar-reflectivity comparison for a 500-m ndown run:
    d2 (500 m WRF) | d1 (2.5 km WRF) on top, HRRR (`refc`) | MRMS
    (`MergedReflectivityQCComposite`, observed) below.

    All four are the same physical quantity (dBZ), a direct apples-to-apples
    comparison -- see `plot_refl_comparison`. Projection / extent / outline /
    BNF-marker conventions are described in the module comment above.
    """
    d2_lon, d2_lat, d2_refl, _ = wrf_reader.read_column_max_reflectivity(d2_file)
    d1_lon, d1_lat, d1_refl, _ = wrf_reader.read_column_max_reflectivity(d1_file)
    hrrr_lon, hrrr_lat, hrrr_refl, _ = hrrr_reader.read_composite_reflectivity(hrrr_file)
    mrms_lon, mrms_lat, mrms_refl, mrms_time = mrms_reader.read_composite_reflectivity(mrms_file)
    proj = wrf_reader.get_lambert_projection(d2_file)

    fig, axes = plt.subplots(
        2, 2, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)

    _nested_four_panel(
        fig, axes,
        (d2_lon, d2_lat, d2_refl, "d2 (500 m) WRF max reflectivity (dBZ)"),
        [
            (d1_lon, d1_lat, d1_refl, "d1 (2.5 km) WRF max reflectivity (dBZ)"),
            (hrrr_lon, hrrr_lat, hrrr_refl, "HRRR composite reflectivity (dBZ)"),
            (
                mrms_lon, mrms_lat, mrms_refl,
                f"MRMS composite reflectivity (dBZ)\n{mrms_time:%Y-%m-%d %H:%M} UTC",
            ),
        ],
        refl_levels, refl_cmap, domain_pad_deg, "max", "dBZ",
    )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight="bold")
    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
    return fig


def _resolve_d1_run_dir(d2_run_dir: Path, d1_run_dir: str | Path | None) -> Path:
    if d1_run_dir is not None:
        return Path(d1_run_dir)
    try:
        return d2_run_dir.parent / _D2_TO_D1_RUN[d2_run_dir.name]
    except KeyError:
        raise ValueError(
            f"no known d1 counterpart for d2 run {d2_run_dir.name!r}; "
            f"pass d1_run_dir= explicitly (known: {sorted(_D2_TO_D1_RUN)})"
        ) from None


def plot_run_d2_tb_comparison_4panel(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    d1_run_dir: str | Path | None = None,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    goes_dir: str | Path = "goes_data",
    channel: int = 13,
    satellite: str = goes_reader.DEFAULT_SATELLITE,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_goes: bool = False,
    **kwargs,
):
    """Wrapper around `plot_d2_tb_comparison_4panel` that takes a time plus a
    500-m ndown run's name/directory and locates the matching d2 wrfout, the
    d1 (2.5-km parent) wrfout, and the HRRR / GOES files itself.

    `run_dir` is the d2 run (e.g.
    `satoshi_testruns/20250917lassobnfwrfhrrr3/rund2`); the d1 panel comes
    from `d1_run_dir`, defaulting to the same case directory's `rund1` (for
    `rund2`) or `rund1-dynlit` (for `rund2-dynlit`). Both the d2 and d1
    runs write the domain as `wrfout_d01_*`, so `domain` applies to both.

    Saved filenames use the `d2_olr_tb_4panel` prefix. Other parameters
    mirror `plot_run_tb_comparison_4panel`.
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    d1_dir = _resolve_d1_run_dir(run_dir, d1_run_dir)
    d2_file = _find_wrf_file(run_dir, domain, time)
    d1_file = _find_wrf_file(d1_dir, domain, time)
    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)
    goes_file = _resolve_goes_file(time, goes_dir, channel, satellite, auto_download_goes)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(output_base_dir, "d2_olr_tb_4panel", run_dir, time)

    suptitle = kwargs.pop("suptitle", f"{run_name} -- {time:%Y-%m-%d %H:%M} UTC")
    return plot_d2_tb_comparison_4panel(
        d2_file, d1_file, hrrr_file, goes_file, suptitle=suptitle, **kwargs
    )


def plot_run_d2_refl_comparison_4panel(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    d1_run_dir: str | Path | None = None,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    mrms_dir: str | Path = "mrms_data",
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_mrms: bool = False,
    **kwargs,
):
    """Wrapper around `plot_d2_refl_comparison_4panel`; the reflectivity
    counterpart of `plot_run_d2_tb_comparison_4panel` -- see there for the
    d1-counterpart resolution and the file-layout conventions. Saved
    filenames use the `d2_refl_4panel` prefix.
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    d1_dir = _resolve_d1_run_dir(run_dir, d1_run_dir)
    d2_file = _find_wrf_file(run_dir, domain, time)
    d1_file = _find_wrf_file(d1_dir, domain, time)
    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)
    mrms_file = _resolve_mrms_file(time, mrms_dir, auto_download_mrms)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(output_base_dir, "d2_refl_4panel", run_dir, time)

    suptitle = kwargs.pop("suptitle", f"{run_name} -- {time:%Y-%m-%d %H:%M} UTC")
    return plot_d2_refl_comparison_4panel(
        d2_file, d1_file, hrrr_file, mrms_file, suptitle=suptitle, **kwargs
    )


# ---------------------------------------------------------------------------
# d2 (500-m) free-run- vs. HRRR-driven 4-panel comparisons
#
# Two 500-m runs on the same grid that differ only in what drives their
# lateral boundaries: the `*hrrr3/rund2` run is ndown-forced (15-min) from
# the free-running 2.5-km d1 run, the `*hrrr4/rund2` run directly by hourly
# HRRRv4 analyses. One frame per output time, laid out as a 2x2 grid:
#
#     free-run-driven d2 (500 m)   |   HRRR-driven d2 (500 m)
#     HRRR                         |   independent observation
#
# Same projection / extent / outline / BNF-marker conventions as the d2
# plots above (the two d2 grids agree to float32 round-off), and like them
# Tb is the OLR-derived Tb for all three model panels -- no CRTM.
# ---------------------------------------------------------------------------

# HRRR-driven case-dir suffix -> the free-run-driven one holding the
# matching `rund2` (same case day, same 500-m grid).
_HRRR_DRIVEN_TO_FREE_CASE = {"hrrr4": "hrrr3"}


def plot_d2_forcing_tb_comparison_4panel(
    free_d2_file: str | Path,
    hrrr_d2_file: str | Path,
    hrrr_file: str | Path,
    goes_file: str | Path,
    ctt_levels: np.ndarray = DEFAULT_CTT_LEVELS,
    cmap: str | Colormap = DEFAULT_TB_CMAP,
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (12.6, 7.9),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """2x2 OLR-derived brightness-temperature comparison of two 500-m runs:
    free-run-driven d2 | HRRR-driven d2 on top, HRRR | GOES ABI ch.13
    (observed) below.

    The three model panels all use the OLR-fit Tb
    (`wrf.read_brightness_temperature_from_olr` /
    `hrrr.read_brightness_temperature_from_olr`) -- no CRTM. See the module
    comment above "d2 (500-m) free-run- vs. HRRR-driven 4-panel
    comparisons"; parameters otherwise mirror `plot_d2_tb_comparison_4panel`.
    """
    free_lon, free_lat, free_tb, _ = wrf_reader.read_brightness_temperature_from_olr(free_d2_file)
    hd2_lon, hd2_lat, hd2_tb, _ = wrf_reader.read_brightness_temperature_from_olr(hrrr_d2_file)
    hrrr_lon, hrrr_lat, hrrr_tb, _ = hrrr_reader.read_brightness_temperature_from_olr(hrrr_file)
    goes_lon, goes_lat, goes_tb, goes_time = goes_reader.read_brightness_temperature(goes_file)
    proj = wrf_reader.get_lambert_projection(free_d2_file)

    fig, axes = plt.subplots(
        2, 2, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)

    _nested_four_panel(
        fig, axes,
        (free_lon, free_lat, free_tb, "free-run-driven d2 (500 m) WRF OLR brightness temp. (K)"),
        [
            (hd2_lon, hd2_lat, hd2_tb, "HRRR-driven d2 (500 m) WRF OLR brightness temp. (K)"),
            (hrrr_lon, hrrr_lat, hrrr_tb, "HRRR OLR brightness temp. (K)"),
            (
                goes_lon, goes_lat, goes_tb,
                f"GOES ABI ch.13 brightness temp. (K)\n{goes_time:%Y-%m-%d %H:%M} UTC scan",
            ),
        ],
        ctt_levels, cmap, domain_pad_deg, "both", "K",
    )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight="bold")
    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
    return fig


def plot_d2_forcing_refl_comparison_4panel(
    free_d2_file: str | Path,
    hrrr_d2_file: str | Path,
    hrrr_file: str | Path,
    mrms_file: str | Path,
    refl_levels: np.ndarray = DEFAULT_REFL_LEVELS,
    refl_cmap: str = "turbo",
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (12.6, 7.9),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """2x2 column-max radar-reflectivity comparison of two 500-m runs:
    free-run-driven d2 | HRRR-driven d2 on top, HRRR (`refc`) | MRMS
    (`MergedReflectivityQCComposite`, observed) below. See
    `plot_d2_forcing_tb_comparison_4panel` and `plot_d2_refl_comparison_4panel`.
    """
    free_lon, free_lat, free_refl, _ = wrf_reader.read_column_max_reflectivity(free_d2_file)
    hd2_lon, hd2_lat, hd2_refl, _ = wrf_reader.read_column_max_reflectivity(hrrr_d2_file)
    hrrr_lon, hrrr_lat, hrrr_refl, _ = hrrr_reader.read_composite_reflectivity(hrrr_file)
    mrms_lon, mrms_lat, mrms_refl, mrms_time = mrms_reader.read_composite_reflectivity(mrms_file)
    proj = wrf_reader.get_lambert_projection(free_d2_file)

    fig, axes = plt.subplots(
        2, 2, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)

    _nested_four_panel(
        fig, axes,
        (free_lon, free_lat, free_refl, "free-run-driven d2 (500 m) WRF max reflectivity (dBZ)"),
        [
            (hd2_lon, hd2_lat, hd2_refl, "HRRR-driven d2 (500 m) WRF max reflectivity (dBZ)"),
            (hrrr_lon, hrrr_lat, hrrr_refl, "HRRR composite reflectivity (dBZ)"),
            (
                mrms_lon, mrms_lat, mrms_refl,
                f"MRMS composite reflectivity (dBZ)\n{mrms_time:%Y-%m-%d %H:%M} UTC",
            ),
        ],
        refl_levels, refl_cmap, domain_pad_deg, "max", "dBZ",
    )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight="bold")
    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
    return fig


def _resolve_free_d2_run_dir(hrrr_d2_run_dir: Path, free_run_dir: str | Path | None) -> Path:
    if free_run_dir is not None:
        return Path(free_run_dir)
    case = hrrr_d2_run_dir.parent.name
    for driven, free in _HRRR_DRIVEN_TO_FREE_CASE.items():
        if case.endswith(driven):
            return hrrr_d2_run_dir.parent.parent / (case[: -len(driven)] + free) / hrrr_d2_run_dir.name
    raise ValueError(
        f"no known free-run-driven counterpart for case {case!r}; pass free_run_dir= "
        f"explicitly (known suffixes: {sorted(_HRRR_DRIVEN_TO_FREE_CASE)})"
    )


def plot_run_d2_forcing_tb_comparison_4panel(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    free_run_dir: str | Path | None = None,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    goes_dir: str | Path = "goes_data",
    channel: int = 13,
    satellite: str = goes_reader.DEFAULT_SATELLITE,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_goes: bool = False,
    **kwargs,
):
    """Wrapper around `plot_d2_forcing_tb_comparison_4panel` that takes a time
    plus the HRRR-driven 500-m run's name/directory and locates both d2
    wrfouts and the HRRR / GOES files itself.

    `run_dir` is the HRRR-driven run (e.g.
    `satoshi_testruns/20250917lassobnfwrfhrrr4/rund2`); the free-run-driven
    panel comes from `free_run_dir`, defaulting to the same subdir of the
    matching `*hrrr3` case directory (`_HRRR_DRIVEN_TO_FREE_CASE`). Both
    runs write the domain as `wrfout_d01_*`, so `domain` applies to both.

    Saved filenames use the `d2_free_vs_hrrr_olr_tb_4panel` prefix. Other
    parameters mirror `plot_run_d2_tb_comparison_4panel`.
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    free_dir = _resolve_free_d2_run_dir(run_dir, free_run_dir)
    hrrr_d2_file = _find_wrf_file(run_dir, domain, time)
    free_d2_file = _find_wrf_file(free_dir, domain, time)
    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)
    goes_file = _resolve_goes_file(time, goes_dir, channel, satellite, auto_download_goes)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(
            output_base_dir, "d2_free_vs_hrrr_olr_tb_4panel", run_dir, time
        )

    suptitle = kwargs.pop("suptitle", f"{run_name} -- {time:%Y-%m-%d %H:%M} UTC")
    return plot_d2_forcing_tb_comparison_4panel(
        free_d2_file, hrrr_d2_file, hrrr_file, goes_file, suptitle=suptitle, **kwargs
    )


def plot_run_d2_forcing_refl_comparison_4panel(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    free_run_dir: str | Path | None = None,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    mrms_dir: str | Path = "mrms_data",
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_mrms: bool = False,
    **kwargs,
):
    """Wrapper around `plot_d2_forcing_refl_comparison_4panel`; the
    reflectivity counterpart of `plot_run_d2_forcing_tb_comparison_4panel` --
    see there for the free-run-counterpart resolution. Saved filenames use
    the `d2_free_vs_hrrr_refl_4panel` prefix.
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    free_dir = _resolve_free_d2_run_dir(run_dir, free_run_dir)
    hrrr_d2_file = _find_wrf_file(run_dir, domain, time)
    free_d2_file = _find_wrf_file(free_dir, domain, time)
    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)
    mrms_file = _resolve_mrms_file(time, mrms_dir, auto_download_mrms)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(
            output_base_dir, "d2_free_vs_hrrr_refl_4panel", run_dir, time
        )

    suptitle = kwargs.pop("suptitle", f"{run_name} -- {time:%Y-%m-%d %H:%M} UTC")
    return plot_d2_forcing_refl_comparison_4panel(
        free_d2_file, hrrr_d2_file, hrrr_file, mrms_file, suptitle=suptitle, **kwargs
    )


# ---------------------------------------------------------------------------
# d3 (100-m) 4-panel comparisons
#
# For the `rund2d3` runs: the 500-m domain (ndown-forced from d1) run
# together with a 100-m domain nested in it online (one-way, `feedback=0`).
# WRF numbers the domains within the run, so the 500-m domain is written as
# `wrfout_d01_*` (the project's "d2") and the 100-m domain as `wrfout_d02_*`
# (the project's "d3"). The d3 domain starts 2 h after d2.
#
# One frame per output time, laid out as a 2x2 grid:
#
#     d3 (100 m, WRF)   |   d2 (500 m, WRF, the same run's parent domain)
#     HRRR              |   independent observation
#
# Same conventions as the d2 plots above, with d3 in d2's place: every panel
# on the (shared) Lambert projection, the map extent the d3 footprint padded
# by `domain_pad_deg` (0.25 here vs. 0.5 for d1/d2, since d3 is much smaller),
# d2 / HRRR / obs cropped to it, the d3 domain outlined
# in red, the BNF site starred. The d2 panel comes from the same `rund2d3`
# run (d3's actual parent), not the separate `rund2` run -- the two runs'
# 500-m fields diverge once d3 starts.
#
# Like the d2 plots (and unlike d1's CRTM-based `tb4` / `wrf_tb`), Tb is the
# cheap OLR-derived Tb for all three model panels; GOES is the only observed
# Tb. Reflectivity is column-max `REFL_10CM` / HRRR `refc` / MRMS composite.
# ---------------------------------------------------------------------------


def plot_d3_tb_comparison_4panel(
    d3_file: str | Path,
    d2_file: str | Path,
    hrrr_file: str | Path,
    goes_file: str | Path,
    ctt_levels: np.ndarray = DEFAULT_CTT_LEVELS,
    cmap: str | Colormap = DEFAULT_TB_CMAP,
    domain_pad_deg: float = 0.25,
    figsize: tuple[float, float] = (9.3, 7.9),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """2x2 OLR-derived brightness-temperature comparison for a 100-m nested
    run: d3 (100 m WRF) | d2 (500 m WRF) on top, HRRR | GOES ABI ch.13
    (observed) below.

    The three model panels all use the OLR-fit Tb
    (`wrf.read_brightness_temperature_from_olr` /
    `hrrr.read_brightness_temperature_from_olr`) -- no CRTM. See the module
    comment above "d3 (100-m) 4-panel comparisons"; parameters otherwise
    mirror `plot_d2_tb_comparison_4panel`.
    """
    d3_lon, d3_lat, d3_tb, _ = wrf_reader.read_brightness_temperature_from_olr(d3_file)
    d2_lon, d2_lat, d2_tb, _ = wrf_reader.read_brightness_temperature_from_olr(d2_file)
    hrrr_lon, hrrr_lat, hrrr_tb, _ = hrrr_reader.read_brightness_temperature_from_olr(hrrr_file)
    goes_lon, goes_lat, goes_tb, goes_time = goes_reader.read_brightness_temperature(goes_file)
    proj = wrf_reader.get_lambert_projection(d3_file)

    fig, axes = plt.subplots(
        2, 2, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)

    _nested_four_panel(
        fig, axes,
        (d3_lon, d3_lat, d3_tb, "d3 (100 m) WRF OLR brightness temp. (K)"),
        [
            (d2_lon, d2_lat, d2_tb, "d2 (500 m) WRF OLR brightness temp. (K)"),
            (hrrr_lon, hrrr_lat, hrrr_tb, "HRRR OLR brightness temp. (K)"),
            (
                goes_lon, goes_lat, goes_tb,
                f"GOES ABI ch.13 brightness temp. (K)\n{goes_time:%Y-%m-%d %H:%M} UTC scan",
            ),
        ],
        ctt_levels, cmap, domain_pad_deg, "both", "K", lon_grid_step=1.0,
    )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight="bold")
    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
    return fig


def plot_d3_refl_comparison_4panel(
    d3_file: str | Path,
    d2_file: str | Path,
    hrrr_file: str | Path,
    mrms_file: str | Path,
    refl_levels: np.ndarray = DEFAULT_REFL_LEVELS,
    refl_cmap: str = "turbo",
    domain_pad_deg: float = 0.25,
    figsize: tuple[float, float] = (9.3, 7.9),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """2x2 column-max radar-reflectivity comparison for a 100-m nested run:
    d3 (100 m WRF) | d2 (500 m WRF) on top, HRRR (`refc`) | MRMS
    (`MergedReflectivityQCComposite`, observed) below. See
    `plot_d3_tb_comparison_4panel` and `plot_d2_refl_comparison_4panel`.
    """
    d3_lon, d3_lat, d3_refl, _ = wrf_reader.read_column_max_reflectivity(d3_file)
    d2_lon, d2_lat, d2_refl, _ = wrf_reader.read_column_max_reflectivity(d2_file)
    hrrr_lon, hrrr_lat, hrrr_refl, _ = hrrr_reader.read_composite_reflectivity(hrrr_file)
    mrms_lon, mrms_lat, mrms_refl, mrms_time = mrms_reader.read_composite_reflectivity(mrms_file)
    proj = wrf_reader.get_lambert_projection(d3_file)

    fig, axes = plt.subplots(
        2, 2, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)

    _nested_four_panel(
        fig, axes,
        (d3_lon, d3_lat, d3_refl, "d3 (100 m) WRF max reflectivity (dBZ)"),
        [
            (d2_lon, d2_lat, d2_refl, "d2 (500 m) WRF max reflectivity (dBZ)"),
            (hrrr_lon, hrrr_lat, hrrr_refl, "HRRR composite reflectivity (dBZ)"),
            (
                mrms_lon, mrms_lat, mrms_refl,
                f"MRMS composite reflectivity (dBZ)\n{mrms_time:%Y-%m-%d %H:%M} UTC",
            ),
        ],
        refl_levels, refl_cmap, domain_pad_deg, "max", "dBZ", lon_grid_step=1.0,
    )

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight="bold")
    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
    return fig


def plot_run_d3_tb_comparison_4panel(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    goes_dir: str | Path = "goes_data",
    channel: int = 13,
    satellite: str = goes_reader.DEFAULT_SATELLITE,
    d3_domain: str = "d02",
    d2_domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_goes: bool = False,
    **kwargs,
):
    """Wrapper around `plot_d3_tb_comparison_4panel` that takes a time plus a
    d2+d3 run's name/directory (e.g.
    `satoshi_testruns/20250917lassobnfwrfhrrr3/rund2d3`) and locates the d3
    (`wrfout_<d3_domain>_*`) and d2 (`wrfout_<d2_domain>_*`) files in it,
    plus the HRRR / GOES files. Saved filenames use the `d3_olr_tb_4panel`
    prefix.
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    d3_file = _find_wrf_file(run_dir, d3_domain, time)
    d2_file = _find_wrf_file(run_dir, d2_domain, time)
    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)
    goes_file = _resolve_goes_file(time, goes_dir, channel, satellite, auto_download_goes)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(output_base_dir, "d3_olr_tb_4panel", run_dir, time)

    suptitle = kwargs.pop("suptitle", f"{run_name} -- {time:%Y-%m-%d %H:%M} UTC")
    return plot_d3_tb_comparison_4panel(
        d3_file, d2_file, hrrr_file, goes_file, suptitle=suptitle, **kwargs
    )


def plot_run_d3_refl_comparison_4panel(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    mrms_dir: str | Path = "mrms_data",
    d3_domain: str = "d02",
    d2_domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_mrms: bool = False,
    **kwargs,
):
    """Wrapper around `plot_d3_refl_comparison_4panel`; the reflectivity
    counterpart of `plot_run_d3_tb_comparison_4panel`. Saved filenames use
    the `d3_refl_4panel` prefix.
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    d3_file = _find_wrf_file(run_dir, d3_domain, time)
    d2_file = _find_wrf_file(run_dir, d2_domain, time)
    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)
    mrms_file = _resolve_mrms_file(time, mrms_dir, auto_download_mrms)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(output_base_dir, "d3_refl_4panel", run_dir, time)

    suptitle = kwargs.pop("suptitle", f"{run_name} -- {time:%Y-%m-%d %H:%M} UTC")
    return plot_d3_refl_comparison_4panel(
        d3_file, d2_file, hrrr_file, mrms_file, suptitle=suptitle, **kwargs
    )


# ---------------------------------------------------------------------------
# Single-panel plots
#
# One field, one map, on the same Lambert projection / map extent / WRF-domain
# outline + BNF site marker as the comparison plots above, and -- crucially for
# stitching frames into a movie -- the same fixed color scale
# (`DEFAULT_CTT_LEVELS` / `DEFAULT_REFL_LEVELS`) for every frame. Used by
# `bnf_wrf_eval.batch_plot` to render whole time series, one PNG per frame,
# which `bnf_wrf_eval.make_movies` then encodes.
#
# The WRF panels (`plot_run_wrf_tb_single`, `plot_run_wrf_refl_single`) are
# per-run, keyed on a run directory like the comparison wrappers. The HRRR /
# GOES / MRMS panels don't depend on the WRF microphysics, so they're keyed on
# a *case* (a `.../<case_dir>/` holding the per-scheme run subdirs) plus a
# `wrf_ref_dir` -- any one run subdir, used only to pull the shared domain
# projection / outline / extent (all time-independent).
# ---------------------------------------------------------------------------


def _first_wrf_file(run_dir: str | Path, domain: str) -> Path:
    """First (time-sorted) `wrfout_<domain>_*` file in `run_dir` -- for
    pulling the domain's projection/outline/extent, which don't change with
    time, without needing a file at one specific time."""
    files = sorted(Path(run_dir).glob(f"wrfout_{domain}_*"))
    if not files:
        raise FileNotFoundError(f"No wrfout_{domain}_* files in {run_dir}")
    return files[0]


def _case_output_path(
    output_base_dir: str | Path, prefix: str, case_label: str, time: dt.datetime
) -> Path:
    output_base_dir = Path(output_base_dir)
    output_base_dir.mkdir(parents=True, exist_ok=True)
    return output_base_dir / f"{prefix}_{case_label}_{time:%Y%m%d_%H%MZ}.png"


def plot_single_field(
    lon: np.ndarray,
    lat: np.ndarray,
    values: np.ndarray,
    *,
    proj,
    domain_lon: np.ndarray,
    domain_lat: np.ndarray,
    levels: np.ndarray,
    cmap: str | Colormap,
    extend: str,
    colorbar_label: str,
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (7.0, 6.0),
    crop: bool = False,
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """Draw one field on one map panel, using the same projection / extent /
    domain-outline / BNF-site conventions as the comparison plots.

    There is no per-panel axes title -- the only text label is `suptitle`
    (the model/obs name plus the field's actual valid/scan time, built by
    the caller). The colorbar is drawn to exactly the height of the
    rendered map: cartopy shrinks a `set_extent` axes to the data's aspect
    ratio on draw, so the colorbar axes is placed from `ax.get_position()`
    *after* a first draw rather than via `constrained_layout` (which can't
    match that shrunk box).

    `proj` / `domain_lon` / `domain_lat` come from the WRF domain (via
    `wrf.get_lambert_projection` and `wrf.read_field(..., "HGT")` or an
    equivalent). `crop=True` slices `lon/lat/values` to the map extent first
    (do this for the wide HRRR / GOES / MRMS grids; leave it False for WRF's
    own grid, which already matches).
    """
    extent = _domain_extent(domain_lon, domain_lat, domain_pad_deg)
    if crop:
        lon, lat, values = _crop_to_extent(lon, lat, values, extent)
    domain_outline = _domain_outline(domain_lon, domain_lat)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.04, 0.04, 0.84, 0.88], projection=proj)
    mesh = _plot_panel(
        ax, lon, lat, values, levels, cmap, extent, "",
        extend=extend, domain_outline=domain_outline,
    )

    # First draw so cartopy applies the map aspect and `ax` settles to the
    # box it actually occupies; then match the colorbar to that box.
    fig.canvas.draw()
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.015, pos.y0, 0.022, pos.height])
    fig.colorbar(mesh, cax=cax, label=colorbar_label)

    if suptitle is not None:
        fig.text(
            0.5, min(pos.y1 + 0.035, 0.99), suptitle,
            ha="center", va="bottom", fontsize=13, fontweight="bold",
        )
    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")

    return fig


def plot_run_wrf_tb_single(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    wrf_tb_kwargs: dict | None = None,
    **kwargs,
):
    """Single-panel WRF CRTM-simulated ABI ch.13 brightness temperature for
    one run/time. Reads from `crtm_cache/` -- pass
    `wrf_tb_kwargs={"cache_dir": ..., "require_cache": True}` so a cache miss
    is a fast, loud error rather than a minutes-long inline CRTM run."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    run_dir = Path(run_dir)
    wrf_file = _find_wrf_file(run_dir, domain, time)

    domain_lon, domain_lat, _, _ = wrf_reader.read_field(wrf_file, "HGT")
    proj = wrf_reader.get_lambert_projection(wrf_file)
    lon, lat, tb, valid_time = crtm_reader.read_simulated_brightness_temperature(
        wrf_file, **(wrf_tb_kwargs or {})
    )

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _output_path(output_base_dir, "wrf_crtm_tb", run_dir, time)
    suptitle = kwargs.pop("suptitle", f"{run_name} -- {valid_time:%Y-%m-%d %H:%M} UTC")
    return plot_single_field(
        lon, lat, tb, proj=proj, domain_lon=domain_lon, domain_lat=domain_lat,
        levels=DEFAULT_CTT_LEVELS, cmap=DEFAULT_TB_CMAP, extend="both",
        colorbar_label="Brightness Temp. (K)",
        suptitle=suptitle, out_file=out_file, **kwargs,
    )


def plot_run_wrf_refl_single(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    **kwargs,
):
    """Single-panel WRF column-max `REFL_10CM` reflectivity for one run/time."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    run_dir = Path(run_dir)
    wrf_file = _find_wrf_file(run_dir, domain, time)

    proj = wrf_reader.get_lambert_projection(wrf_file)
    lon, lat, refl, valid_time = wrf_reader.read_column_max_reflectivity(wrf_file)

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _output_path(output_base_dir, "wrf_refl", run_dir, time)
    suptitle = kwargs.pop("suptitle", f"{run_name} -- {valid_time:%Y-%m-%d %H:%M} UTC")
    levels = kwargs.pop("refl_levels", DEFAULT_REFL_LEVELS)
    return plot_single_field(
        lon, lat, refl, proj=proj, domain_lon=lon, domain_lat=lat,
        levels=levels, cmap=kwargs.pop("refl_cmap", "turbo"), extend="max",
        colorbar_label="Comp. Refl (dBZ)",
        suptitle=suptitle, out_file=out_file, **kwargs,
    )


def plot_run_wrf_olr_tb_single(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    **kwargs,
):
    """Single-panel WRF OLR-derived IR-window brightness temperature for one
    run/time.

    Uses `wrf.read_brightness_temperature_from_olr` -- the broadband
    Yang & Slingo (2001) ``OLR = sigma*Tf**4`` inversion of WRF's TOA `OLR`
    field, the same cheap conversion the 4-panel ``wrf_panel="olr"``
    comparison and the `tb4_simple` frames use. No forward model and no
    `crtm_cache/`; shares the `DEFAULT_TB_CMAP` / `DEFAULT_CTT_LEVELS`
    scale with the CRTM and `ctt` Tb panels so the series is comparable."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    run_dir = Path(run_dir)
    wrf_file = _find_wrf_file(run_dir, domain, time)

    proj = wrf_reader.get_lambert_projection(wrf_file)
    lon, lat, tb, valid_time = wrf_reader.read_brightness_temperature_from_olr(wrf_file)

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _output_path(output_base_dir, "wrf_olr_tb", run_dir, time)
    suptitle = kwargs.pop("suptitle", f"{run_name} -- {valid_time:%Y-%m-%d %H:%M} UTC")
    return plot_single_field(
        lon, lat, tb, proj=proj, domain_lon=lon, domain_lat=lat,
        levels=DEFAULT_CTT_LEVELS, cmap=DEFAULT_TB_CMAP, extend="both",
        colorbar_label="Brightness Temp. (K)",
        suptitle=suptitle, out_file=out_file, **kwargs,
    )


def plot_case_hrrr_tb_single(
    time: dt.datetime | str,
    case_label: str,
    wrf_ref_dir: str | Path,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    **kwargs,
):
    """Single-panel HRRR CRTM-simulated brightness temperature (`SBT114`) for
    one case/time, cropped to the WRF domain's map extent."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    wrf_ref_file = _first_wrf_file(wrf_ref_dir, domain)
    proj = wrf_reader.get_lambert_projection(wrf_ref_file)
    domain_lon, domain_lat, _, _ = wrf_reader.read_field(wrf_ref_file, "HGT")

    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)
    lon, lat, tb, valid_time = hrrr_reader.read_simulated_brightness_temperature(hrrr_file)

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _case_output_path(output_base_dir, "hrrr_tb", case_label, time)
    suptitle = kwargs.pop("suptitle", f"HRRR -- {valid_time:%Y-%m-%d %H:%M} UTC")
    return plot_single_field(
        lon, lat, tb, proj=proj, domain_lon=domain_lon, domain_lat=domain_lat,
        levels=DEFAULT_CTT_LEVELS, cmap=DEFAULT_TB_CMAP, extend="both",
        colorbar_label="Brightness Temp. (K)",
        crop=True, suptitle=suptitle, out_file=out_file, **kwargs,
    )


def plot_case_hrrr_refl_single(
    time: dt.datetime | str,
    case_label: str,
    wrf_ref_dir: str | Path,
    hrrr_base_dir: str | Path = "satoshi_forcing_data/hrrr/hrrrnat_data",
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    **kwargs,
):
    """Single-panel HRRR composite reflectivity (`refc`) for one case/time."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    wrf_ref_file = _first_wrf_file(wrf_ref_dir, domain)
    proj = wrf_reader.get_lambert_projection(wrf_ref_file)
    domain_lon, domain_lat, _, _ = wrf_reader.read_field(wrf_ref_file, "HGT")

    hrrr_file = _find_hrrr_file(hrrr_base_dir, time)
    lon, lat, refl, valid_time = hrrr_reader.read_composite_reflectivity(hrrr_file)

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _case_output_path(output_base_dir, "hrrr_refl", case_label, time)
    suptitle = kwargs.pop("suptitle", f"HRRR -- {valid_time:%Y-%m-%d %H:%M} UTC")
    levels = kwargs.pop("refl_levels", DEFAULT_REFL_LEVELS)
    return plot_single_field(
        lon, lat, refl, proj=proj, domain_lon=domain_lon, domain_lat=domain_lat,
        levels=levels, cmap=kwargs.pop("refl_cmap", "turbo"), extend="max",
        colorbar_label="Comp. Refl (dBZ)",
        crop=True, suptitle=suptitle, out_file=out_file, **kwargs,
    )


def plot_case_goes_tb_single(
    time: dt.datetime | str,
    case_label: str,
    wrf_ref_dir: str | Path,
    goes_dir: str | Path = "goes_data",
    channel: int = 13,
    satellite: str = goes_reader.DEFAULT_SATELLITE,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_goes: bool = False,
    **kwargs,
):
    """Single-panel GOES ABI observed brightness temperature for the scan
    closest to `time`. The saved filename carries the *requested* `time`
    (a regular grid), not the scan time, so a movie's frames stay evenly
    spaced even where two grid steps snap to the same scan."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    wrf_ref_file = _first_wrf_file(wrf_ref_dir, domain)
    proj = wrf_reader.get_lambert_projection(wrf_ref_file)
    domain_lon, domain_lat, _, _ = wrf_reader.read_field(wrf_ref_file, "HGT")

    goes_file = _resolve_goes_file(time, goes_dir, channel, satellite, auto_download_goes)
    lon, lat, tb, scan_time = goes_reader.read_brightness_temperature(goes_file)

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _case_output_path(output_base_dir, "goes_tb", case_label, time)
    sat_num = satellite.lstrip("Gg")
    suptitle = kwargs.pop(
        "suptitle", f"GOES-{sat_num} (#{channel}) -- {scan_time:%Y-%m-%d %H:%M} UTC"
    )
    return plot_single_field(
        lon, lat, tb, proj=proj, domain_lon=domain_lon, domain_lat=domain_lat,
        levels=DEFAULT_CTT_LEVELS, cmap=DEFAULT_TB_CMAP, extend="both",
        colorbar_label="Brightness Temp. (K)",
        crop=True, suptitle=suptitle, out_file=out_file, **kwargs,
    )


def plot_case_mrms_refl_single(
    time: dt.datetime | str,
    case_label: str,
    wrf_ref_dir: str | Path,
    mrms_dir: str | Path = "mrms_data",
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_mrms: bool = False,
    **kwargs,
):
    """Single-panel MRMS composite reflectivity for the scan closest to
    `time`; filename carries the requested `time` (see
    `plot_case_goes_tb_single`)."""
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)
    wrf_ref_file = _first_wrf_file(wrf_ref_dir, domain)
    proj = wrf_reader.get_lambert_projection(wrf_ref_file)
    domain_lon, domain_lat, _, _ = wrf_reader.read_field(wrf_ref_file, "HGT")

    mrms_file = _resolve_mrms_file(time, mrms_dir, auto_download_mrms)
    lon, lat, refl, scan_time = mrms_reader.read_composite_reflectivity(mrms_file)

    out_file = kwargs.pop("out_file", None)
    if out_file is None and output_base_dir is not None:
        out_file = _case_output_path(output_base_dir, "mrms_refl", case_label, time)
    suptitle = kwargs.pop("suptitle", f"MRMS -- {scan_time:%Y-%m-%d %H:%M} UTC")
    levels = kwargs.pop("refl_levels", DEFAULT_REFL_LEVELS)
    return plot_single_field(
        lon, lat, refl, proj=proj, domain_lon=domain_lon, domain_lat=domain_lat,
        levels=levels, cmap=kwargs.pop("refl_cmap", "turbo"), extend="max",
        colorbar_label="Comp. Refl (dBZ)",
        crop=True, suptitle=suptitle, out_file=out_file, **kwargs,
    )


def plot_lightning_cg_comparison(
    wrf_file_t1: str | Path,
    wrf_file_t2: str | Path,
    mrms_dir: str | Path,
    auto_download_mrms: bool = True,
    count_levels: np.ndarray = DEFAULT_LIGHTNING_COUNT_LEVELS,
    cmap: str = "YlGn_r",
    dark_background: bool = True,
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (11, 5.5),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """Compare cloud-to-ground flash counts between WRF's Dynamic Lightning
    Scheme and MRMS/NLDN observations, over the interval between two
    consecutive wrfout times.

    WRF panel: `cg_pos + cg_neg` from `wrf.read_dyn_lightning_flash_counts`
    (the difference of `LPOS`/`LNEG` between `wrf_file_t1` and
    `wrf_file_t2`), on the run's native grid.

    MRMS panel: `lightning.mrms_cg_counts_on_wrf_grid` -- every ~1km NLDN
    cell's value in `(t1, t2]`, binned and summed into the *same* WRF grid
    cells, so both panels share one grid/projection/extent and are directly
    differenceable, not just visually side-by-side.

    Flash counts are sparse (mostly 0) and non-negative, unlike the smooth
    dBZ/K fields the other comparison plots use, so cells with a count of
    exactly 0 are left unfilled (masked to NaN before contouring) rather
    than colored as the bottom `count_levels` bin, and `count_levels`
    defaults to a small set of discrete count thresholds rather than an
    evenly-spaced range.

    This is a first-cut, not-yet-validated comparison -- see the caveat in
    `mrms.read_nldn_cg_density` about the AvgDensity-to-count conversion,
    and `prompts/lightning_evaluation.md` for the broader plan.

    Parameters
    ----------
    wrf_file_t1, wrf_file_t2 : two consecutive wrfout_d0X_* files from a
        `dyn_lightning_option=1` run (e.g. `rund1-dynlit`).
    mrms_dir : local directory to look for/download NLDN files into.
    suptitle : if given, drawn as a big figure-level title.
    out_file : if given, the figure is saved there.

    Returns
    -------
    matplotlib.figure.Figure
    """
    wrf_lon, wrf_lat, cg_pos, cg_neg, _ic, (t1, t2) = wrf_reader.read_dyn_lightning_flash_counts(
        wrf_file_t1, wrf_file_t2
    )
    wrf_cg = cg_pos + cg_neg
    proj = wrf_reader.get_lambert_projection(wrf_file_t1)

    mrms_lon, mrms_lat, mrms_cg, mrms_scan_times = lightning_reader.mrms_cg_counts_on_wrf_grid(
        wrf_file_t1, t1, t2, mrms_dir, auto_download=auto_download_mrms
    )

    extent = _domain_extent(wrf_lon, wrf_lat, domain_pad_deg)
    domain_outline = _domain_outline(wrf_lon, wrf_lat)

    wrf_cg_masked = np.where(wrf_cg > 0, wrf_cg, np.nan)
    mrms_cg_masked = np.where(mrms_cg > 0, mrms_cg, np.nan)

    fig, axes = plt.subplots(
        1, 2, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)
    if dark_background:
        fig.patch.set_facecolor("black")

    n_scans = len(mrms_scan_times)
    mrms_title = f"MRMS/NLDN\n({n_scans} 1-min scans summed)"
    _plot_row(
        fig, axes,
        [
            (wrf_lon, wrf_lat, wrf_cg_masked, "WRF"),
            (mrms_lon, mrms_lat, mrms_cg_masked, mrms_title),
        ],
        count_levels, cmap, extent, domain_outline, extend="max",
        colorbar_label="CG flashes per grid cell", plot_type="pcolormesh",
        dark_background=dark_background,
    )

    if suptitle is not None:
        fig.suptitle(
            suptitle, fontsize=16, fontweight="bold", y=0.93,
            color="white" if dark_background else "black",
        )

    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")

    return fig


def plot_lightning_total_comparison(
    wrf_file_t1: str | Path,
    wrf_file_t2: str | Path,
    goes_dir: str | Path,
    satellite: str = goes_reader.DEFAULT_SATELLITE,
    auto_download_glm: bool = True,
    count_levels: np.ndarray = DEFAULT_LIGHTNING_TOTAL_COUNT_LEVELS,
    cmap: str = "YlGn_r",
    dark_background: bool = True,
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (11, 5.5),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """Compare total-lightning (CG + intracloud) flash counts between WRF's
    Dynamic Lightning Scheme and GOES GLM observations, over the interval
    between two consecutive wrfout times.

    WRF panel: `cg_pos + cg_neg + ic` from
    `wrf.read_dyn_lightning_flash_counts` (the difference of
    `LPOS`/`LNEG`/`LNEU` between `wrf_file_t1` and `wrf_file_t2`), on the
    run's native grid.

    GLM panel: `lightning.glm_total_counts_on_wrf_grid` -- every flash
    centroid from GLM-L2-LCFA files in `(t1, t2]`, binned and summed into
    the *same* WRF grid cells. GLM (like WRF's scheme here) reports total
    lightning without distinguishing CG from intracloud, so this is a
    genuinely apples-to-apples comparison, unlike the CG-only comparison in
    `plot_lightning_cg_comparison`.

    Same rendering conventions as `plot_lightning_cg_comparison` --
    pcolormesh (not contourf, which silently drops isolated sparse cells),
    zero-count cells left unfilled, dark background by default. See that
    function's docstring for why.

    This is a first-cut, not-yet-validated comparison -- see
    `prompts/lightning_evaluation.md` for the broader plan.

    Parameters
    ----------
    wrf_file_t1, wrf_file_t2 : two consecutive wrfout_d0X_* files from a
        `dyn_lightning_option=1` run (e.g. `rund1-dynlit`).
    goes_dir : local directory to look for/download GLM-L2-LCFA files into.
    suptitle : if given, drawn as a big figure-level title.
    out_file : if given, the figure is saved there.

    Returns
    -------
    matplotlib.figure.Figure
    """
    wrf_lon, wrf_lat, cg_pos, cg_neg, ic, (t1, t2) = wrf_reader.read_dyn_lightning_flash_counts(
        wrf_file_t1, wrf_file_t2
    )
    wrf_total = cg_pos + cg_neg + ic
    proj = wrf_reader.get_lambert_projection(wrf_file_t1)

    glm_lon, glm_lat, glm_total, glm_scan_times = lightning_reader.glm_total_counts_on_wrf_grid(
        wrf_file_t1, t1, t2, goes_dir, satellite=satellite, auto_download=auto_download_glm
    )

    extent = _domain_extent(wrf_lon, wrf_lat, domain_pad_deg)
    domain_outline = _domain_outline(wrf_lon, wrf_lat)

    wrf_total_masked = np.where(wrf_total > 0, wrf_total, np.nan)
    glm_total_masked = np.where(glm_total > 0, glm_total, np.nan)

    fig, axes = plt.subplots(
        1, 2, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0.02, hspace=0.02)
    if dark_background:
        fig.patch.set_facecolor("black")

    n_scans = len(glm_scan_times)
    glm_title = f"GLM ({satellite})\n({n_scans} 20-s scans summed)"
    _plot_row(
        fig, axes,
        [
            (wrf_lon, wrf_lat, wrf_total_masked, "WRF"),
            (glm_lon, glm_lat, glm_total_masked, glm_title),
        ],
        count_levels, cmap, extent, domain_outline, extend="max",
        colorbar_label="Total flashes per grid cell", plot_type="pcolormesh",
        dark_background=dark_background,
    )

    if suptitle is not None:
        fig.suptitle(
            suptitle, fontsize=16, fontweight="bold", y=0.93,
            color="white" if dark_background else "black",
        )

    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")

    return fig


def plot_run_lightning_cg_comparison(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    mrms_dir: str | Path = "mrms_data",
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_mrms: bool = False,
    **kwargs,
):
    """Wrapper around `plot_lightning_cg_comparison` that takes an hourly
    window's start time plus a WRF run's name/directory, locates the
    matching pair of consecutive-hour wrfout files itself, and labels the
    figure with a suptitle showing `run_name` and the window.

    Parameters
    ----------
    time : the window's *start* time; the run's `wrfout_<domain>_*` files
        at both `time` and `time + 1h` must exist (so this must be an
        on-the-hour wrfout time, and not the run's very last one).
    run_name : label for the run, used in the figure's suptitle.
    run_dir : directory holding that run's `wrfout_*` files -- must be a
        `dyn_lightning_option=1` run (e.g. `.../rund1-dynlit`).
    mrms_dir : directory holding downloaded NLDN files (see
        `mrms.download_mrms_file`). Defaults to `mrms_data`.
    domain : WRF domain to plot, e.g. "d01" (default).
    output_base_dir : if given (and `out_file` isn't passed explicitly via
        `**kwargs`), auto-saves to
        `<output_base_dir>/lightning_cg_wrf_mrms_<run_dir's last two path
        components>_<time>.png`.
    auto_download_mrms : fetch missing NLDN files from S3 instead of
        skipping them. False by default, so plotting never triggers a
        network call unless explicitly asked -- important under a SLURM
        array on a compute node with no outbound internet.
    **kwargs : forwarded to `plot_lightning_cg_comparison`.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    time_end = time + dt.timedelta(hours=1)
    wrf_file_t1 = _find_wrf_file(run_dir, domain, time)
    wrf_file_t2 = _find_wrf_file(run_dir, domain, time_end)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(output_base_dir, "lightning_cg_wrf_mrms", run_dir, time)

    suptitle = kwargs.pop(
        "suptitle",
        f"CG lightning: {run_name}\n{time:%Y-%m-%d %H:%M}-{time_end:%H:%M} UTC",
    )
    return plot_lightning_cg_comparison(
        wrf_file_t1, wrf_file_t2, mrms_dir,
        auto_download_mrms=auto_download_mrms, suptitle=suptitle, **kwargs,
    )


def plot_run_lightning_total_comparison(
    time: dt.datetime | str,
    run_name: str,
    run_dir: str | Path,
    goes_dir: str | Path = "goes_data",
    satellite: str = goes_reader.DEFAULT_SATELLITE,
    domain: str = "d01",
    output_base_dir: str | Path | None = None,
    auto_download_glm: bool = False,
    **kwargs,
):
    """Wrapper around `plot_lightning_total_comparison` -- GLM's
    counterpart of `plot_run_lightning_cg_comparison`; see that function's
    docstring for the shared conventions (hourly window, output-path
    naming, why `auto_download_glm` defaults to False).

    **kwargs : forwarded to `plot_lightning_total_comparison`.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if isinstance(time, str):
        time = dt.datetime.fromisoformat(time)

    run_dir = Path(run_dir)
    time_end = time + dt.timedelta(hours=1)
    wrf_file_t1 = _find_wrf_file(run_dir, domain, time)
    wrf_file_t2 = _find_wrf_file(run_dir, domain, time_end)

    if output_base_dir is not None and "out_file" not in kwargs:
        kwargs["out_file"] = _output_path(output_base_dir, "lightning_total_wrf_glm", run_dir, time)

    suptitle = kwargs.pop(
        "suptitle",
        f"Total lightning: {run_name}\n{time:%Y-%m-%d %H:%M}-{time_end:%H:%M} UTC",
    )
    return plot_lightning_total_comparison(
        wrf_file_t1, wrf_file_t2, goes_dir, satellite=satellite,
        auto_download_glm=auto_download_glm, suptitle=suptitle, **kwargs,
    )


def plot_lightning_nalma_comparison(
    wrf_file_t1: str | Path,
    wrf_file_t2: str | Path,
    nalma_dir: str | Path,
    token_file: str | Path = nalma_reader.DEFAULT_TOKEN_FILE,
    auto_download_nalma: bool = True,
    wrf_count_levels: np.ndarray = DEFAULT_LIGHTNING_TOTAL_COUNT_LEVELS,
    nalma_source_levels: np.ndarray = DEFAULT_LIGHTNING_SOURCE_COUNT_LEVELS,
    cmap: str = "YlGn_r",
    dark_background: bool = True,
    domain_pad_deg: float = 0.5,
    figsize: tuple[float, float] = (11.5, 5.5),
    suptitle: str | None = None,
    out_file: str | Path | None = None,
):
    """Compare WRF's total-lightning flash count against NALMA's raw VHF
    source count, over the interval between two consecutive wrfout times.

    WRF panel: `cg_pos + cg_neg + ic` from
    `wrf.read_dyn_lightning_flash_counts`, same as
    `plot_lightning_total_comparison`'s WRF panel.

    NALMA panel: `lightning.nalma_source_counts_on_wrf_grid` -- every VHF
    source from NALMA granules in `(t1, t2]`, binned and summed into the
    *same* WRF grid cells.

    **Not a like-for-like count comparison** -- unlike GLM, NALMA's raw
    data are individual VHF sources, not flash-clustered totals (a flash
    produces many tens to hundreds of sources; see `nalma.py`'s module
    docstring). Checked against a real active-storm hour, NALMA sources
    outnumbered WRF flashes ~51:1 domain-wide -- sharing one color scale
    between the two panels would wash out the WRF side entirely, so
    (unlike the CG/GLM comparisons) **each panel gets its own levels and
    colorbar**: `wrf_count_levels` (same flash-count scale as
    `plot_lightning_total_comparison`) and `nalma_source_levels` (a much
    wider range for raw source counts). Read this as a spatial/qualitative
    check -- does NALMA's activity line up with WRF's and GLM's -- rather
    than a magnitude comparison. No chi^2/power quality filtering is
    applied to the raw sources yet (see `nalma.read_nalma_sources`).

    Requires a NASA Earthdata Login token to download NALMA granules
    (unlike GLM/MRMS's anonymous public buckets) -- see `nalma.py`'s module
    docstring.

    Same rendering conventions as `plot_lightning_cg_comparison`/
    `plot_lightning_total_comparison` -- pcolormesh, zero-count cells left
    unfilled, dark background by default.

    Parameters
    ----------
    wrf_file_t1, wrf_file_t2 : two consecutive wrfout_d0X_* files from a
        `dyn_lightning_option=1` run (e.g. `rund1-dynlit`).
    nalma_dir : local directory to look for/download NALMA granules into.
    suptitle : if given, drawn as a big figure-level title.
    out_file : if given, the figure is saved there.

    Returns
    -------
    matplotlib.figure.Figure
    """
    wrf_lon, wrf_lat, cg_pos, cg_neg, ic, (t1, t2) = wrf_reader.read_dyn_lightning_flash_counts(
        wrf_file_t1, wrf_file_t2
    )
    wrf_total = cg_pos + cg_neg + ic
    proj = wrf_reader.get_lambert_projection(wrf_file_t1)

    nalma_lon, nalma_lat, nalma_counts, granule_names = (
        lightning_reader.nalma_source_counts_on_wrf_grid(
            wrf_file_t1, t1, t2, nalma_dir, token_file=token_file,
            auto_download=auto_download_nalma,
        )
    )

    extent = _domain_extent(wrf_lon, wrf_lat, domain_pad_deg)
    domain_outline = _domain_outline(wrf_lon, wrf_lat)

    wrf_total_masked = np.where(wrf_total > 0, wrf_total, np.nan)
    nalma_counts_masked = np.where(nalma_counts > 0, nalma_counts, np.nan)

    fig, axes = plt.subplots(
        1, 2, figsize=figsize, subplot_kw={"projection": proj}, constrained_layout=True
    )
    fig.get_layout_engine().set(w_pad=0.05, h_pad=0.02, wspace=0.05, hspace=0.02)
    if dark_background:
        fig.patch.set_facecolor("black")

    n_granules = len(granule_names)
    mesh_wrf = _plot_panel(
        axes[0], wrf_lon, wrf_lat, wrf_total_masked, wrf_count_levels, cmap, extent,
        "WRF (flashes)", extend="max", domain_outline=domain_outline,
        plot_type="pcolormesh", dark_background=dark_background,
    )
    mesh_nalma = _plot_panel(
        axes[1], nalma_lon, nalma_lat, nalma_counts_masked, nalma_source_levels, cmap, extent,
        f"NALMA\n({n_granules} granules summed)", extend="max", domain_outline=domain_outline,
        plot_type="pcolormesh", dark_background=dark_background,
    )
    cbar_wrf = fig.colorbar(mesh_wrf, ax=axes[0], label="Flashes per grid cell", shrink=0.5, pad=0.02)
    cbar_nalma = fig.colorbar(
        mesh_nalma, ax=axes[1], label="VHF sources per grid cell", shrink=0.5, pad=0.02
    )
    if dark_background:
        for cbar, label in (
            (cbar_wrf, "Flashes per grid cell"), (cbar_nalma, "VHF sources per grid cell")
        ):
            cbar.set_label(label, color="white")
            cbar.ax.yaxis.set_tick_params(color="white")
            plt.setp(plt.getp(cbar.ax, "yticklabels"), color="white")

    if suptitle is not None:
        fig.suptitle(
            suptitle, fontsize=16, fontweight="bold", y=0.93,
            color="white" if dark_background else "black",
        )

    if out_file is not None:
        fig.savefig(out_file, dpi=150, bbox_inches="tight")

    return fig

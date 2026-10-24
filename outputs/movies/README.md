# BNF WRF evaluation movies

Animated comparisons of the WRF LASSO BNF test runs (Satoshi's `hrrr3` runs — the
d1 (2.5-km) runs, the d2 (500-m) runs nested into them via `ndown`, and the d3 (100-m)
runs nested online inside d2 — plus the `hrrr4` 500-m runs driven directly by hourly
HRRRv4 analyses) against the
HRRR analyses that drove them and against independent satellite / radar
observations, for the DOE ARM Bankhead National Forest (BNF) site in northern
Alabama.

All 102 movies (72 d1 + 12 d2 + 6 d2 free-run- vs. HRRR-driven + 6 d3 + 6 lightning)
are H.264 (`yuv420p`), no audio.
Frame rate depends on the time step between frames:

- **Hourly movies, 2 fps** (half a second per hour): the d1 comparison movies (73
  frames, ~36 s), the HRRR single-panel movies (73 frames, ~36 s), the d2 (15–24
  frames — those runs cover only 14–23 h; see below), d2 free-run- vs. HRRR-driven
  (same 15–24 frames) and d3 (13–21 frames) comparison movies (~6.5–12 s), and the
  lightning movies (72 frames, ~36 s — one shorter than a
  normal hourly movie, since each frame differences *two* consecutive hours and the
  run's last hour has no successor to pair with; see below).
- **15-min and 5-min movies, 12 fps**: the WRF and MRMS single-panel movies (289
  frames, ~24 s) and the GOES movies (~865 frames, ~72 s).

## The runs

Three case periods, each a 72-hour WRF forecast on a single 2.5-km domain
(550 x 390 grid, Lambert Conformal), initialized from HRRR at 00 UTC and writing
output every 15 minutes:

| Case label | Simulated period (UTC) |
|------------|------------------------|
| `20250502` | 2025-05-01 00Z – 2025-05-04 00Z |
| `20250520` | 2025-05-19 00Z – 2025-05-22 00Z |
| `20250917` | 2025-09-16 00Z – 2025-09-19 00Z |

Each case was run with three microphysics schemes (everything else identical):

| Run tag       | `mp_physics` | Scheme |
|---------------|--------------|--------|
| `rund1`       | 28 | Thompson aerosol-aware |
| `rund1-mp52`  | 52 | P3 (1 ice category) |
| `rund1-mp53`  | 53 | P3 (2 ice categories) |

A fourth run per case, `rund1-dynlit`, is identical to `rund1` (Thompson,
`mp_physics=28`) but with WRF's dynamics-based lightning-threat diagnostic
enabled (`dyn_lightning_option=1`, `lightning_option=0`). That diagnostic is
passive — it does not feed back into microphysics or radiation — so its
brightness-temperature and reflectivity fields are essentially identical to
`rund1`; the point of evaluating it is to confirm that. Only the two cheap
`wrfout`-derived single-panel series (`wrf_olr_tb`, `wrf_refl`) were produced
for `rund1-dynlit`; it has no CRTM Tb and no comparison movies.

### d2 (500-m) runs

Each case also has one or more **500-m runs** nested into the d1 domain with WRF's
`ndown` one-way nesting (the ndown output domain is written as `wrfout_d01_*`, so
these too are single-domain from the plotting side). The movies here cover two per
case:

| Run tag        | `mp_physics` | Notes |
|----------------|--------------|-------|
| `rund2`        | 28 | Thompson aerosol-aware, 500 m |
| `rund2-dynlit` | 28 | as `rund2` + passive dynamics-based lightning-threat diagnostic (`dyn_lightning_option=1`); Tb / reflectivity essentially identical to `rund2` |

Unlike the 72-hour d1 runs, these d2 runs cover only part of each case period
(roughly its convective window):

| Case | d2 period (UTC) | hourly frames |
|------|-----------------|---------------|
| `20250502` | 05-02 07Z – 05-03 06Z | 24 |
| `20250520` | 05-20 12Z – 05-21 05Z | 18 |
| `20250917` | 09-17 10Z – 09-18 00Z | 15 |

### HRRR-driven d2 (500-m) runs

The `hrrr4` case directories each have a `rund2` on **the same 500-m grid** as the
`hrrr3` `rund2` (same projection, Thompson `mp_physics=28`, same output times and
period as the table above), but with its lateral boundaries **driven directly by the
hourly HRRRv4 analyses** instead of by the free-running 2.5-km d1 run via `ndown`.
In these movies "free-run-driven d2" is `hrrr3/rund2` and "HRRR-driven d2" is
`hrrr4/rund2`.

### d3 (100-m) runs

`rund2d3` (one per case) runs the 500-m d2 domain (ndown-forced, as `rund2`) together
with a **100-m domain nested in it online** (one-way, `feedback = 0`). WRF numbers
domains within the run, so the 500-m domain is `wrfout_d01_*` and the 100-m domain
(the project's "d3") is `wrfout_d02_*`. d3 starts 2 h after d2. The d2 panel in these
movies is this run's own 500-m domain (d3's actual parent), **not** `rund2`: the two
runs' 500-m fields are identical until d3 starts and then diverge (cause not yet
understood).

| Case | d3 period (UTC) | hourly frames |
|------|-----------------|---------------|
| `20250502` | 05-02 09Z – 05-03 06Z | 21 (no 11Z: d2 output missing at a restart) |
| `20250520` | 05-20 14Z – 05-21 05Z | 15 (no 22Z: same reason) |
| `20250917` | 09-17 12Z – 09-18 00Z | 13 |

## Reference data

- **GOES**: GOES-19 ABI channel 13 (10.3 µm "clean window" IR) brightness
  temperature, CONUS scans (~5 min).
- **HRRR**: hourly HRRR analyses — `refc` composite reflectivity and `SBT114`
  (UPP CRTM–simulated ABI ch.13 brightness temperature).
- **MRMS**: `MergedReflectivityQCComposite` composite reflectivity (sampled every
  15 min here).

## What the WRF quantities are

- **WRF CRTM Tb**: ABI ch.13 brightness temperature computed from WRF's own state
  with the CRTM forward model (the same approach UPP uses for HRRR's `SBT114`), so
  the WRF/HRRR/GOES Tb comparison is genuinely like-for-like (same forward model,
  same channel).
- **WRF `ctt`**: WRF's simpler built-in cloud-top-temperature diagnostic
  (`wrf-python`), shown alongside the CRTM panel for reference.
- **WRF simple Tb (OLR fit)**: brightness temperature obtained by inverting the
  broadband Yang & Slingo (2001) `OLR = σ·Tf⁴`, `Tf = Tb·(a + b·Tb)` relation on
  WRF's top-of-atmosphere `OLR` field (the same conversion PyFLEXTRKR uses). No
  forward model — a cheap alternative to the CRTM panel. Being a column-integrated
  flux it runs slightly warm relative to a true window channel and smooths out the
  coldest convective cloud tops. Also produced as a stand-alone single-panel
  series (`wrf_olr_tb`) for every run, on the same 180–315 K Tb colour scale as
  the CRTM and `ctt` panels. **The d2 comparison movies use this OLR-fit Tb for
  every model panel — d2, d1 *and* HRRR** (`hrrr.read_brightness_temperature_from_olr`,
  the same conversion applied to HRRR's `OLR`) — so no CRTM is involved there. This
  differs from the d1 `*_tb_comparison` movies, which use CRTM-derived Tb for WRF
  against HRRR's own CRTM `SBT114`.
- **WRF reflectivity**: column maximum of the simulated `REFL_10CM` field (dBZ),
  directly comparable to the HRRR and MRMS composite reflectivity.

## Common plotting conventions

Every panel is reprojected onto the **same** Lambert Conformal projection and map
extent (taken from the WRF domain), with the WRF domain footprint outlined in
**red** and the BNF site marked with a **red star**. For the d1 movies the
projection / extent / red outline are the 2.5-km domain's; for the **d2 movies**
they are the **500-m domain's** (its footprint padded slightly), so the d1 panel
there shows the 2.5-km run subset to the d2 area. Colour scales are fixed across
all frames so the animations do not flicker:

- Brightness temperature: 180–315 K, 5-K steps, conventional IR enhancement
  (grayscale for 240–315 K, then a cyan-to-magenta rainbow for the coldest
  convective cloud tops down to 180 K).
- Composite / column-max reflectivity: 5–75 dBZ, 5-dBZ steps.

## Movie inventory

### Comparison movies (hourly, 73 frames each, 2 fps) — 27 files

| File pattern | Layout |
|--------------|--------|
| `wrf_ctt_crtm_hrrr_goes_tb_comparison_<case>_<run>.mp4` | 2×2 brightness temperature: WRF `ctt` \| WRF CRTM \| HRRR CRTM \| GOES observed |
| `wrf_simple_tb_crtm_hrrr_goes_tb_comparison_<case>_<run>.mp4` | 2×2 brightness temperature, as above but with the top-left panel replaced by **WRF simple Tb (OLR fit)** |
| `wrf_hrrr_mrms_refl_comparison_<case>_<run>.mp4` | 3-panel reflectivity: WRF \| HRRR \| MRMS |

9 of each of the three (3 cases × 3 runs).

### Single-panel movies — 45 files

One field per frame, at its native cadence.

| File pattern | Field | Cadence / frames | Count |
|--------------|-------|------------------|-------|
| `wrf_crtm_tb_<case>_<run>.mp4` | WRF CRTM ABI ch.13 Tb | 15 min / 289 | 9 |
| `wrf_refl_<case>_<run>.mp4` | WRF column-max reflectivity | 15 min / 289 | 12 |
| `wrf_olr_tb_<case>_<run>.mp4` | WRF OLR-fit brightness temperature | 15 min / 289 | 12 |
| `hrrr_tb_<case>.mp4` | HRRR CRTM Tb (`SBT114`) | hourly / 73 | 3 |
| `hrrr_refl_<case>.mp4` | HRRR composite reflectivity (`refc`) | hourly / 73 | 3 |
| `goes_tb_<case>.mp4` | GOES-19 ABI ch.13 Tb | ~5 min / ~865 | 3 |
| `mrms_refl_<case>.mp4` | MRMS composite reflectivity | 15 min / 289 | 3 |

`wrf_refl` and `wrf_olr_tb` have 12 files each: `rund1` / `rund1-mp52` /
`rund1-mp53` / `rund1-dynlit` × 3 cases. `wrf_crtm_tb` has 9 (no `rund1-dynlit`).

The observation-only single-panel movies (`hrrr_*`, `goes_tb`, `mrms_refl`) depend
only on the case, not on the microphysics run, so their filenames have no run tag.
Their suptitle carries the field's *actual* valid/scan time (e.g. a GOES frame
labelled 20:02 is the scan nearest the requested 20:00).

### d2 (500-m) comparison movies (hourly, 15–24 frames each, 2 fps) — 12 files

2×2 layout: **d2 (500 m WRF) | d1 (2.5 km WRF)** on top, **HRRR | independent
observation** below. All on the d2 domain's projection / extent, d2 footprint
outlined in red, BNF star on every panel.

| File pattern | Layout |
|--------------|--------|
| `d2_olr_tb_4panel_<case>_<run>.mp4` | Brightness temperature, **OLR-fit for all three model panels** (d2 \| d1 \| HRRR), GOES observed bottom-right — no CRTM |
| `d2_refl_4panel_<case>_<run>.mp4` | Reflectivity: WRF d2 column-max \| WRF d1 column-max \| HRRR `refc` \| MRMS composite |

6 of each (3 cases × `rund2` / `rund2-dynlit`). Frame counts per case: 24
(`20250502`), 18 (`20250520`), 15 (`20250917`).

### d2 free-run- vs. HRRR-driven comparison movies (hourly, 15–24 frames each, 2 fps) — 6 files

2×2 layout: **free-run-driven d2 (500 m, `hrrr3/rund2`) | HRRR-driven d2 (500 m,
`hrrr4/rund2`)** on top, **HRRR | independent observation** below. Same projection /
extent / red d2 outline / BNF star as the d2 movies above.

| File pattern | Layout |
|--------------|--------|
| `d2_free_vs_hrrr_olr_tb_4panel_<case>lassobnfwrfhrrr4_rund2.mp4` | Brightness temperature, **OLR-fit for all three model panels** (free-run-driven d2 \| HRRR-driven d2 \| HRRR), GOES observed bottom-right — no CRTM |
| `d2_free_vs_hrrr_refl_4panel_<case>lassobnfwrfhrrr4_rund2.mp4` | Reflectivity: free-run-driven d2 column-max \| HRRR-driven d2 column-max \| HRRR `refc` \| MRMS composite |

3 of each (one per case). Frame counts per case: 24 (`20250502`), 18 (`20250520`),
15 (`20250917`).

### d3 (100-m) comparison movies (hourly, 13–21 frames each, 2 fps) — 6 files

2×2 layout: **d3 (100 m WRF) | d2 (500 m WRF, same `rund2d3` run)** on top, **HRRR |
independent observation** below (no d1 panel). All on the d3 domain's projection,
extent = d3 footprint padded by 0.25° (half the d1/d2 padding), d3 footprint outlined
in red, BNF star on every panel.

| File pattern | Layout |
|--------------|--------|
| `d3_olr_tb_4panel_<case>_rund2d3.mp4` | Brightness temperature, **OLR-fit for all three model panels** (d3 \| d2 \| HRRR), GOES observed bottom-right — no CRTM |
| `d3_refl_4panel_<case>_rund2d3.mp4` | Reflectivity: WRF d3 column-max \| WRF d2 column-max \| HRRR `refc` \| MRMS composite |

### Lightning movies (hourly, 72 frames each, 2 fps) — 6 files

WRF's dynamics-based lightning-threat diagnostic (`rund1-dynlit`,
`dyn_lightning_option=1`; Lynn et al. 2012) against two independent observed
lightning sources, for each case's `rund1-dynlit` run only (the other three
d1 schemes don't have the diagnostic enabled). Unlike the Tb/reflectivity
comparisons above, lightning is inherently an *accumulated-count* quantity,
not an instantaneous field: each frame differences WRF's cumulative flash
counts (`LPOS`/`LNEG`/`LNEU`, confirmed to be running totals since simulation
start, not per-output-interval values) between two consecutive on-the-hour
`wrfout` times, and sums the observation source's raw data over that same
hour, so the two sides are counting flashes over an identical window.

| File pattern | Layout | WRF quantity | Obs. quantity |
|--------------|--------|--------------|----------------|
| `lightning_cg_wrf_mrms_<case>_rund1-dynlit.mp4` | 2-panel: WRF \| MRMS/NLDN | Cloud-to-ground flash count (`LPOS + LNEG` difference) | MRMS `NLDN_CG_001min_AvgDensity` (ground-based NLDN CG network), every 1-min sample in the hour binned onto the WRF grid and summed |
| `lightning_total_wrf_glm_<case>_rund1-dynlit.mp4` | 2-panel: WRF \| GLM | Total (CG + intracloud) flash count (`LPOS + LNEG + LNEU` difference) | GOES-19 GLM-L2-LCFA flash centroids (satellite optical total-lightning detection), every ~20-s scan in the hour binned onto the WRF grid and summed |

3 of each (one per case). Both use a **black background** (`YlGn_r`: dark
green for common low counts, brighter yellow-green for higher counts, forced
to white for any cell exceeding the top color level, so the most intense
storm-core cells can never blend into the background) rather than this
project's usual white-background map style — flash counts are sparse and
mostly-zero, and pcolormesh (not contourf, which silently drops isolated
single-cell values surrounded by no-data neighbors) renders each grid cell on
its own. `lightning_total`'s color levels run higher (up to 200) than
`lightning_cg`'s (up to 100), since total lightning (CG + intracloud) is a
larger count than CG alone.

A raw-VHF-source (not flash-clustered) comparison against NALMA also exists
(`plotting.plot_lightning_nalma_comparison`) but isn't part of this movie set
yet — NALMA source counts run roughly 50x flash counts (a flash produces many
VHF sources as its channel propagates), and converting sources to flash/
cluster counts (matching the standard LMA flash-clustering algorithm) hasn't
been done; see `prompts/lightning_evaluation.md`.

## Filename key

    <type>_<case>lassobnfwrfhrrr3[_<run>].mp4
    d2_free_vs_hrrr_<field>_4panel_<case>lassobnfwrfhrrr4_rund2.mp4

- `<case>` — `20250502`, `20250520`, or `20250917` (see tables above)
- `<run>` — for the d1 movies: `rund1`, `rund1-mp52`, `rund1-mp53`, or
  `rund1-dynlit` (WRF single-panel movies; d1 comparison movies are `rund1` /
  `rund1-mp52` / `rund1-mp53` only, and `rund1-dynlit` has only `wrf_olr_tb` /
  `wrf_refl`, plus both `lightning_*` movies, which are `rund1-dynlit` only).
  For the `d2_*_4panel` movies: `rund2` or `rund2-dynlit`; for the `d3_*_4panel`
  movies: `rund2d3`. The `d2_free_vs_hrrr_*` movies are named after the
  HRRR-driven run (`hrrr4/rund2`); their free-run-driven panel is the same case's
  `hrrr3/rund2`.

# Hand-off: "fulltest: refl movies"

Written 2026-11-15 at the end of session "fulltest: comparison plan" (in `bnf_wrf_da`). The plan behind it: `../bnf_wrf_da/notes/fulltest_comparison_plan.md` (Sec. 1 decisions, Sec. 3 how the fields were made, Sec. 4 what they show at 05-20 19:00).

## Task

Three movies comparing the reflectivity of the DA run's ensemble-mean analyses with HRRR, the free run and MRMS, in the style of `outputs/movies/wrf_hrrr_mrms_refl_comparison_20250520lassobnfwrfhrrr3_rund1.mp4` (`batch_plot` type `refl3`, `plotting.plot_run_refl_comparison`): each frame has four panels, in this order:

1. the DA analysis (one reflectivity definition per movie, below);
2. HRRR composite reflectivity (`refc`, f00 analysis);
3. the free run's column-max `REFL_10CM`;
4. MRMS composite reflectivity.

Everything else as `refl3`: WRF Lambert projection and extent from the domain, `DEFAULT_REFL_LEVELS` (5–75 dBZ by 5), `turbo`, `contourf`, `extend="max"`, the domain outline in red, the BNF star, one shared colorbar, a suptitle with the valid time, MRMS's scan time in its panel title.

| movie | DA panel | variable in the DA files |
|---|---|---|
| (i) | DART's mean of `REFL_10CM` | `REFL_COLMAX_ANL_DART` |
| (ii) | mean of the members' recomputed reflectivity, in dBZ¹ | `REFL_COLMAX_ANL_MEMDBZ` |
| (iii) | recomputed reflectivity of the mean state | `REFL_COLMAX_ANL_STATE` |

¹ The mean in dBZ, not in linear Z (user, 2026-11-15).

**Frames**: hourly, 05-20 01:00 to 05-21 06:00 UTC, 30 frames per movie (HRRR has analyses only on the hour; 05-20 00:00 is the DA cold start, with no analysis). 2 fps, as the hourly movies (`make_movies`, `HOURLY_TYPES`).

## The DA run, in short

50-member WRF-DART EAKF run `_conv_refl_f2` (`../bnf_wrf_da/cases/20250520/da_theta_m0_n50/`), 2.5 km, the same domain and grid as the free run, 15-min cycles from a 05-20 00Z cold start, 120 analyses 05-20 00:15 to 05-21 06:00. Assimilated: MADIS conventional obs, mesonet and MRMS 3-D reflectivity. The evaluation looks at the analyses, i.e. the states WRF starts each 15-min segment from. All three fields are ensemble means; the column maximum is always taken of the 3-D mean, never a mean of the members' column maxima (user, 2026-11-15).

## What the three DA fields are

- **(i) `ANL_DART`**: column maximum of DART's `postassim_mean.nc` `REFL_10CM`, the mean over the members of DART's regression-updated reflectivity, averaged in dBZ. It is not computed from the hydrometeors and WRF never sees it. It can lie below −35 and above 75 dBZ: over the 120 times 74 columns exceed 75 dBZ, at most 87.6 dBZ (05-21 00:15).
- **(ii) `ANL_MEMDBZ`**: each member's reflectivity recomputed from its analysed hydrometeors (WRF's own Thompson `calc_refl10cm`, run offline), averaged over the members in dBZ, with members without echo at WRF's floor of −35 dBZ. It shows where nearly all members agree: with 45-dBZ echo in some members and none in the others, the mean exceeds 30 dBZ only if at least 81 % of them have the echo.
- **(iii) `ANL_STATE`**: the reflectivity of the mean of the members' analysed fields (the state an LES forced by the mean analysis would start from). It spreads every member's cells, at a fraction of their strength, over all their positions: a broad smear wherever the members disagree on where the cells are.
- Context: HRRR's f00 `refc` is nearly the MRMS picture (HRRR assimilates radar, so it is not independent of MRMS); the free run has far too much convection.

## Data

- **DA**: `/gpfs/wolf2/arm/cli120/proj-shared/hengxiao80/bnf_wrf_da/cases/20250520/da_theta_m0_n50/diag_conv_refl_f2/thompson_refl_mean/thompson_refl_mean_d01_<YYYY-MM-DD_HH_MM_SS>.nc`, one file per analysis time (120 files, every 15 min; use the hourly ones). netCDF-4, 2-D float32 fields (`Time` = 1, `south_north` = 389, `west_east` = 549) in dBZ, with `Times`, `XLAT`, `XLONG` and a member wrfout's global attributes, so `wrf.get_lambert_projection` works on these files (checked: `STAND_LON`, `MOAD_CEN_LAT`, `TRUELAT1`, `TRUELAT2` present). The column maxima are already computed. Eight fields, `REFL_COLMAX_<ANL|FCST>_<STATE|MEMDBZ|MEMZ|DART>`; these movies need only the three above (`FCST_*` is the 15-min forecast valid at the same time, `MEMZ` the members' mean in linear Z). A reader in this repo's convention, returning `(lon, lat, values, valid_time)`, is needed.
- **Free run**: `/gpfs/wolf2/arm/cli120/proj-shared/hengxiao80/bnf_wrf_da/cases/20250520/theta_m0/wrfout_d01_<YYYY-MM-DD_HH_MM_SS>` (15-min, 05-19 00Z to 05-22 00Z; run `theta_m0` in `configs/da_20250520.json`); `wrf.read_column_max_reflectivity`. It is the DA's control: the members' setup (`hybrid_opt = 0`, `use_theta_m = 0`), started 05-19 00Z, so a 25- to 54-h forecast in these frames. It is not sey's `rund1` of the reference movie (`hybrid_opt = 2`, `use_theta_m = 1`). Its `REFL_10CM` is 0 instead of −35 dBZ on the outermost row and column (WRF doesn't compute it there); invisible with levels from 5 dBZ.
- **HRRR**: `satoshi_forcing_data/hrrr/hrrrnat_data/<YYYYMMDD>/hrrr.t<HH>z.wrfnatf00.grib2`, `hrrr.read_composite_reflectivity` (as `refl3`).
- **MRMS**: `mrms_data/`, `plotting._resolve_mrms_file` (local files; 05-20 19:00 resolves to the 19:00:41 scan).

## Checks for the test images

Coverage above 30 dBZ at 05-20 19:00, measured in `bnf_wrf_da` (`run_compare`'s method: MRMS and HRRR pixels inside the WRF domain box via `_inside`; the DA fields over the interior columns, the free run over all columns): MRMS 2.13 %, HRRR 2.17 %, free run 5.66 %, (i) 3.24 %, (ii) 4.11 %, (iii) 10.66 %. Test figures with these fields at 19:00, for comparison: `../bnf_wrf_da/cases/20250520/da_theta_m0_n50/diag_conv_refl_f2/figures/test_mean_of_refl_vs_refl_of_mean_202505201900.png` and `test_refl_stage_202505201900.png`.

## To settle with the user

- Layout: one row of four (like `refl3` widened) or 2 × 2.
- Panel titles (e.g. "DA analysis: DART mean", "DA analysis: mean of members' reflectivity", "DA analysis: reflectivity of the mean state"); whether to print coverage numbers in the titles.
- Movie and frame names, and where they go (e.g. `outputs/compare/...` like `run_compare`, or `outputs/movies/`).
- Hourly only, or also 15-min frames (DA, free run and MRMS every 15 min, HRRR only on the hour); hourly is what was asked for.
- Whether to add the members' mean in linear Z (`REFL_COLMAX_ANL_MEMZ`) or the forecast stage as further movies; they are in the files.
- Not part of this hand-off: the Tb comparison (decided: the OLR fit on the members' forecast OLR, labelled as the forecast; the ensemble mean of OLR is not settled), surface and upper-air fields (not compared for now).

## Working rules

- `../bnf_wrf_da` is read only from here (never write into its `output_conv_refl_f2/`, `run_conv_refl_f2/` or `diag_conv_refl_f2/`).
- Load the `dataviz` skill before chart code. Show a few rendered test images (e.g. 19:00 of each movie) before any batch; production frames and movies, any `sbatch`, and commits each need the user's go-ahead.
- No quantitative guesses (derive from a measurement or say unknown); don't propose checks whose logic can't separate what they claim to. Status lines before tool calls are full sentences. No LaTeX in chat; notes per `~/.claude/rules/typora-markdown.md`.

## Starting the session

Name it `fulltest: refl movies` and name this file in the first prompt (`notes/handoffs/2026-11-15_fulltest_refl_movies.md`). The `fulltest` group was added to this repo's `.claude-sessions.toml` on 2026-11-15 (user), so the name files the session there.

# Hand-off: "fulltest: 09-17 comparisons"

Written 2026-11-23 at the end of session "fulltest_0917: products for comparisons and ndown" (in `bnf_wrf_da`). That session made the ensemble-mean files of the second DA run, case 2025-09-17, as for 05-20. Its record is `../bnf_wrf_da/notes/fulltest_0917_steps.md`: Part D (the run and its radar misfit) and Part E (the files below).

## Task

The DA reflectivity comparisons made for 05-20 (commit `5f7bf35`; `CLAUDE.md`, section "DA (WRF-DART) ensemble-mean reflectivity 4-panel comparisons"), now for 09-17, with the same layouts and definitions:
- three 4-panel movies, DA analysis mean | HRRR `refc` (f00) | free run column-max `REFL_10CM` | MRMS, one per DA definition (`da_refl_dart`, `da_refl_memdbz`, `da_refl_state`);
- the stages movie, DART's prior | posterior mean and forecast | analysis mean state (`da_refl_stages`).

Background, unchanged from 05-20: `../bnf_wrf_da/notes/fulltest_comparison_plan.md` (Sec. 1 decisions, Sec. 3 how the fields are made).

## The 09-17 run and its files

All under `../bnf_wrf_da/cases/20250917/` (read-only from here).
- **Run**: `da_theta_m0_n50/`, `_conv_refl_f2`, the 05-20 setup on the same 2.5-km grid (550 × 390, the d1 grid). Cold start 09-17 04:00, 92 analyses every 15 min, 04:15 to 09-18 03:00.
- **Reflectivity means**: `da_theta_m0_n50/diag_conv_refl_f2/thompson_refl_mean/thompson_refl_mean_d01_<YYYY-MM-DD_HH_MM_SS>.nc`, 92 files, made by the same code as 05-20's, so with the same eight fields `REFL_COLMAX_<ANL|FCST>_<STATE|MEMDBZ|MEMZ|DART>`, names and layout. Checked: all finite, −35.0 to 62.9 dBZ; only `ANL_DART` goes above 75 dBZ (3 columns over all times, at most 83.8).
- **Free run** (the DA's control, as 05-20's `theta_m0`): `theta_m0/`, `hybrid_opt = 0`, `use_theta_m = 0`, from 09-16 00Z, 289 `wrfout_d01` files every 15 min to 09-19 00Z.
- **Ensemble-mean analysis `wrfout`s** (for ndown; not needed for these movies): `da_theta_m0_n50/diag_conv_refl_f2/ens_mean_wrfout/wrfout_d01_<t>`, 93 files (the 04:00 cold start and the analyses), 3.95 GB each.
- **Observations here** (checked 2026-11-23): MRMS composites in `mrms_data/` every 15 min on 09-17 and 09-18 (96 files per day); HRRR native-level files in `satoshi_forcing_data/hrrr/hrrrnat_data/20250917/` and `20250918/`.

## What needs changing here

- `batch_plot.DA_RUNS` holds only the 05-20 run (`_DA_CASE_20250520`). The 09-17 entry: the mean directory above, with `theta_m0/` as the free run.
- **The frame and movie names don't carry the case.** `plotting._output_path` labels a DA frame by the last two components of its diagnostics directory, `da_theta_m0_n50_diag_conv_refl_f2`, the same in both cases, and one batch call goes through every `DA_RUNS` entry into one `OUTPUT_ROOT`. So a run as for 05-20 would put the 09-17 frames into 05-20's frames root, and the movies would get 05-20's names. How to keep the cases apart (e.g. a per-case `OUTPUT_ROOT` such as `outputs/compare/20250917lassobnfwrfhrrr3` with `--start`/`--end`, or a case choice in the batch) is to be settled before the batch run.
- Frame times: the hourly analyses the free run also has, 09-17 05:00 to 09-18 03:00, i.e. 23 frames per type (04:00 is the cold start, with no reflectivity file).
- `CLAUDE.md`'s DA section describes 05-20 only.

## What to expect in the plots

- **Far less echo than on 05-20.** Share of interior columns above 30 dBZ in the analysis, from the files' logs:

| time | mean state | members, dBZ |
|---|---|---|
| 09-17 10:00 | 0.04 % | 0.01 % |
| 15:00 | 0.11 % | 0.03 % |
| 21:00 | 6.21 % | 2.44 % |
| 09-18 00:00 | 2.52 % | 0.62 % |
| 03:00 | 0.37 % | 0.04 % |

- **The radar obs fit poorly where echoes are few and isolated**, i.e. most of the LES window (500-m window 09-17 10Z to 09-18 00Z). From 10Z to ~17Z, 37–73 % of the rain obs per cycle are rejected as outliers, and nearly all of the rejected ones are stronger than every member: small, mostly elevated showers the ensemble doesn't have. In the afternoon convection (19Z–00Z) the rejected share falls to 17 %. So MRMS should show small cells in the morning and the LES window that the DA panels lack.

## Other open items for this repo (none started; the user decides)

From the list at the end of Part D of `../bnf_wrf_da/notes/fulltest_0917_steps.md`:
- **Item 1**: an independent check that the missed cells are rain (GOES cloud tops, lightning), at 09:30–10:00 and through the LES window. Where they are: `../bnf_wrf_da/cases/20250917/da_theta_m0_n50/diag_conv_refl_f2/figures/d2_refl_outlier_map.png` (08:00, 09:30, 10:00), and in Part D, "They look like rain" (Alabama and Georgia, 82.5–86.5° W, 32.8–35.9° N; two cells near 93.5° W 32.5° N and 92° W 31° N). Only if they aren't rain would the observation sets in `bnf_wrf_da` change.
- **Item 5**: check (a) of `bnf_wrf_da`'s 09-17 `hybrid2` against sey's `20250917lassobnfwrfhrrr3/rund1-dynlit`: `run_compare.py` `timeseries` with a 09-17 config (as `configs/da_20250520.json`), at all 289 times. And the correction of `CLAUDE.md` line 579: "bit-identical" holds for the model state and for the two fields compared there, not for the dynamic-lightning fields: the energies `LIGHTNING_PE`, `LIGHTNING_NE`, `LIGHTNING_NEU` differ from +15 min, and by +72 h the flash densities too (`LNEU` on 09-17; `LNEG` and `LNEU` on 05-20). Source: `fulltest_0917_steps.md` P4, check (b). The lightning evaluation here reads those flash counts.

## Working rules

- The user's go-ahead before every `sbatch` and every production run of plots; when plots change, a few test frames first. Commits only after the user has reviewed; push only when asked.
- `bnf_wrf_da` is read-only from here.
- Status lines are full sentences (what, on what, why). No LaTeX in chat. Notes per `~/.claude/rules/typora-markdown.md`.
- No guessed times or sizes: use measured values or say "unknown".

## Starting the session

Name it `fulltest: 09-17 comparisons` (group `fulltest` in `.claude-sessions.toml`), and name this file in the first prompt (`notes/handoffs/2026-11-23_fulltest_0917_comparisons.md`).

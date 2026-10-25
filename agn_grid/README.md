# Recipe: CLOUDY AGN grids

**Status: done (24 Oct 2026).** All 19 grids are built, checked and installed (23 Oct), and
two grids rebuilt from the kept CLOUDY runs are identical to the installed ones (§5).

## 1. Product
synthesizer AGN grids from CLOUDY c23.01: the standard-BLR and gas-enshrouded cores, plus 17 one-axis
variation grids, each per erg/s of incident L_bol (`WeightVariable = bolometric_luminosities`).
They are installed in `SYNTHESIZER_DIR/grids/agn/photoionised/agn-cloudy/`, with the
incident grids they were built from in `grids/agn/incident/agn-cloudy/` (synthesizer needs only
the former; the latter are there to rebuild or extend the grids). The design (axes, counts, run
settings, timing pilot) is in `cloudy_grid_plan.txt`, and the record of what was run and why
is in `cloudy_grid_report.md`. Published in smock's Google Drive folder (`../README.md`).

**Names.** The grids ship as `..._cloudy-standard_blr_*` and `..._cloudy-gas_enshrouded_*`
(author, 22 Oct 2026; both are BLR-type grids in synthesizer's terms). Everything upstream
keeps the internal run-set names `blr` and `dense`: the YAMLs in `params/`, `make_cloudy_inputs.py`'s
run sets, and the work and build directories in `cloudy_runs/`. `build_agn_grids.py` renames
each finished file (step 7, `SHIPPED_RUN_SET_NAMES`) and records the internal name in its
`smock_internal_name` attribute.

**Raw data.** The author keeps the CLOUDY runs behind the shipped grids on Vera (about 53 GB;
not public; ask if you need them): the whole work directory, with every run a shipped grid
uses, the successful runs at cut grid points, and one example per failure class and rerun type
(the rest were pruned on 22 Oct 2026, `archive/2026-10-22_cloudy_runs_prune/`).

**Finding the run behind a grid point.** Each build directory, `cloudy_runs/build/<grid>/`, has
a `MANIFEST.tsv` with one row per model of the shipped grid: `incident_index`, `index` (its
place in the shipped, cut grid), `orig_incident_index`, `orig_index` (its place in the full run
set), `source` and `flags` (1 rerun, 2 nudged, 4 thick_es, 8 hole, 16 masers_off, 32
limit_cycle, 64 filled; the same values as the grid's `smock_flags`). The run is at
`grids/<grid>/<orig_incident_index>/<orig_index>.*` when `source` is `production`, at
`reruns/<source>/<grid>/<orig_incident_index>/<orig_index>.*` when it is a rerun stage (e.g.
`fix_test_4/nudge`), and nowhere when it is `filled` (every run of that model failed, and the
grid holds a fill along T). Or follow the link `build/<grid>/final/<grid>/<incident_index>/<index>.out`.
`<grid>` is the internal name (see "Names").

## 2. Requirements
- smock's env (`sharable_mock_obs` on Vera), with smock installed (`pip install --no-deps -e /path/to/smock`).
- syncretize (https://github.com/synthesizer-project/syncretize, imported as
  `synthesizer_grids`), installed at commit `e646ac9` (30 Jul 2026) with
  `pip install --no-deps .`. Its dependency list includes `fsps`, which needs Fortran and
  `SPS_HOME`; the recipe doesn't use it.
- CLOUDY c23.01; set `CLOUDY_DIR` to its folder (the one holding `source/cloudy.exe` and
  `data/`) before submitting `run_cloudy.slurm`.
- Compute: production took ~7,000 core-hours (21,000 runs); see the report for the reruns.

## 3. Steps
**Work directory.** Every script defaults to `cloudy_runs/` relative to the current
directory (`--work-dir` overrides it where a script has one), so run the commands from the
top of this repository; new CLOUDY work starts in a fresh `cloudy_runs/` (not committed). To
read a kept copy of old runs without writing to it, point `--work-dir` at it; the builder then
needs `--build-dir` so that it writes elsewhere. The kept run lists (`*_runs.txt`,
`reruns/*/runs*.txt`) name each model directory by its absolute path when it ran, a path that
no longer exists. The builder only uses their part below `reruns/`. To re-parse an old run
list, swap that prefix for the current work directory first (e.g. with `sed` into a new file).

1. `make_incident_grid.py`: the incident grids (`--dry-run`, then `--apply`).
2. `make_cloudy_inputs.py <set>`: stage the CLOUDY decks for a run set with syncretize.
3. `run_cloudy.slurm`: run every deck in a run list (**submits SLURM jobs**; ask first;
   throttle arrays with `%10`).
4. `parse_cloudy_runs.py`: classify every run (`ok`, failure class, `thick_es`).
5. `make_rerun_inputs.py`: stage edited copies of failed runs in `cloudy_runs/reruns/<name>/`,
   then step 3 on their run list. The decks already carry every fix that is a no-op for healthy
   runs (`make_cloudy_inputs.py`: no `turbulence 0`, `failures 100`, `set pressure ionize
   300000`, `iterate to convergence 100`). Of production's failures inside the kept T ranges,
   those defaults would leave two kinds that need a rerun (report §1a, §4): maser aborts
   (`--add-command "database H-like masers off"`) and deterministic floating-point crashes
   (`--replace-line` nudging `hden` by +0.001 dex, the BLR one also with masers off).
   Limit cycles and runs that fail every fix are handled by the builder; the rest lie outside
   the kept T ranges and are cut. These reruns are manual stages for now; a script that selects
   and stages them is a possible future item (author, 24 Oct 2026).
6. `build_agn_grids.py`: production plus the latest ok reruns → synthesizer grids, with
   per-model flags, the remaining holes filled along T (dry run by default; the dry run
   lists each hole's stencil). The core grids need ~6 GB (all 19 grids took 26 min in one
   job, 1391044): run them through SLURM (`build_agn_grids.slurm`, **submits a job**;
   `GRIDS=` builds a subset).
7. `install_agn_grids.py`: copy the built grids and their incident grids into
   smock's grids folder (dry run by default; never replaces a file unless `--overwrite`).

The full commands for each stage are in each script's usage lines (and, for the shipped grids, in the ROADMAP status log).

## 4. Problems solved
See `cloudy_grid_report.md`: failure classes (§4), the fixes and their no-op checks (§5), and
the two syncretize problems handled by the builder (§6: the per-L☉ normalisation and
the NZONE-blind failure test).

Solved in the builder since:
- **Holes left after every fix** (4 dense-core runs; synthesizer would interpolate their
  zero spectra): filled by a cubic in log T, flagged `filled`, with a leave-one-out error in
  `smock_fill_error` (author, 18 Oct 2026). Three details matter, each found by a leave-one-out
  test (`archive/2026-10-21_hole_fill_components/`): zeros are propagated, not floored (a
  floored zero with a negative cubic weight blew one fill up by 113%); the transmitted is
  clipped to its two bracketing neighbours (it overshoots where heavily absorbed); and the
  nebular is filled as nebular continuum + linecont, so the continuum can't go negative.
  Everything is filled per ionising photon, then rescaled, so each filled model gets its
  own syncretize `normalisation`.
- **syncretize's model order is not C order** (21 Oct 2026). With several photoionisation
  axes, unravelling the model number in C order puts per-model datasets at the wrong grid
  points. The builder takes each model's grid index from syncretize's own index list
  (`model_grid`). Grids with one photoionisation axis were never affected.

## 5. Verification
**Checks on the built grids** (report §6; `archive/2026-10-21_core_grid_checks/`,
`archive/2026-10-23_gas_enshrouded_core_v5_checks/`), all passed on all 19 grids: the
normalisation is within 1.4×10⁻⁴ of the incident grid; 1,100 sampled models (every rerun and
limit cycle included) equal their own CLOUDY runs to 2×10⁻¹⁶; `smock_flags` matches
`MANIFEST.tsv`; there are no failures or empty models; the fill's leave-one-out errors are
1.5–3.5%. The installed files are byte-identical to the builds.

**Reproduction.** CLOUDY is not rerun for this (author, 24 Oct 2026): the grid may change with
the author's work on the CLOUDY AGN incident shape, and the whole grid could be rerun in about
a week of jobs if needed. Instead, a variation grid and a core grid are rebuilt from the kept
copy of the runs and compared with the installed files (author, 24 Oct 2026). **Result (24 Oct
2026, job 1392405; `archive/2026-10-24_agn_grid_reproduction/`):** the gas-enshrouded core
(117 reruns, 4 filled holes) and the standard-BLR turbulence variation (15 reruns) rebuilt from
the kept runs are identical to the installed grids in every dataset and attribute, apart from the
two build stamps (`smock_build`, syncretize's `date_created`).

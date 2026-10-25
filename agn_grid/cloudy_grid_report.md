# CLOUDY AGN grid report

**Status: updated 24 Oct 2026.** Every run and decision is in; every grid is built, checked and installed. The recipe reproduces the installed grids from the kept copy of the runs (§7).

Times are US Eastern (EDT), as SLURM prints them. Paths starting `cloudy_runs/` are in the work directory, which the author keeps on Vera (not public; see the README, "Raw data").

This is the record of how the AGN photoionisation grids were made: what was run, what went wrong, what was done about it, and the evidence for each decision. The recipe that rebuilds the grids is in this folder, `smock-recipes/agn_grid/`. The session-by-session detail is in the ROADMAP status log of the author's build record (14–23 Oct).

## 1. Summary
- **Runs:** 21,000 CLOUDY c23.01 models in four run sets (BLR core, dense core, and the variations of each), about 7,000 core-hours. 20,437 (97.3%) succeeded first time.
- **Failures:** the 563 flagged runs fall into seven classes (§4). Most sit at log T ≲ 4.2, Z☉ and high column density.
- **Cuts and fixes:**
  - The grids drop the lowest T values: BLR keeps log T 4.25–6.5, dense keeps 4.2–6.5.
  - Inside those ranges, failures were rerun with targeted deck changes (§5).
- **Unfixable runs:**
  - 68 runs never converge, because their line optical depths cycle. Each keeps its production spectrum, flagged, with the cycle's amplitude stored as its error: BLR median 8.7%, dense 2.8% (§4).
  - 4 dense-core runs failed every fix, and are filled from their T neighbours, flagged (§6).
- **Final grids:** production outputs, with failed runs replaced by successful reruns. Built by `build_agn_grids.py`, which runs syncretize's builder on a symlink view, fixes a normalisation bug in it and fills the last holes (§6).

## 1a. Issues at a glance
Each issue, what was tried, and where it ended. The details and evidence are in §3–§6. `archive/2026-10-22_cloudy_runs_prune/history.csv` lists every affected grid point.
- **v_turb = 0 refused** (20 runs). CLOUDY rejects `turbulence 0`. *Tried:* omitting the command. *Solved:* 20/20 ok; it is now the recipe default.
- **Not converged in 10 iterations** (289). *Tried:* 30 iterations (5/10 converged), then 100 (69 more). *Ended:* the remaining 68 are limit cycles, where a line optical depth oscillates for ever. The cycle test confirmed it. They keep their production spectrum, flagged, with the cycle's half-amplitude stored as the error.
- **Maser aborts** (101). *Tried:* `database H-like masers off`. *Solved:* 3/3 inside the kept range; it is rerun-only, because it changes healthy runs at the 1% level. The rest lie at log T ≤ 4.2 and are cut.
- **NZONE limit** (77). *Ended:* all outside the kept T range, so cut. (The pilot's NZONE failures led to `nend` 6000.)
- **Pressure–ionisation aborts** (64). *Tried:* `set pressure ionize` 30000, which fixed 37. *Then:* 300000 with `failures 100`, which fixed 21 of the 25 left. Both commands are now defaults (no-op test). *Ended:* the last 4 are holes, filled along T.
- **Floating-point crashes** (10). *Tried:* a plain rerun, which is deterministic and crashed again. *Solved:* nudging n_H by 0.001 dex (BLR also with masers off): 2/2 inside the kept range. The rest are cut.
- **ConvFail, the 20-failure limit** (2). *Solved:* `failures 100`, 2/2; now a default.
- **Holes** (4 dense-core points). *Ended:* a cubic in log T through 4 neighbours, with every part clipped to its two bracketing neighbours.
- **syncretize's double normalisation.** Its output came out 3.83×10³³ too bright. *Solved:* the builder works around it.
- **syncretize's failure test is too weak.** It accepts NZONE-truncated runs. *Solved:* only `ok` runs enter the build.
- **Builder index bug.** syncretize's model order is not C order. *Solved:* the builder takes syncretize's own indices.
- **Operational problems.** A SLURM TaskProlog failure killed a few tasks: rerun them from their own run list. `set presioniz` didn't parse: parse-check new commands with a short run. The QOS caps jobs at 2 days.

## 2. The grids
- **Incident SED:** CLOUDY's `agn` continuum, ported to `make_incident_grid.py` (`parse_agn.cpp`, `ffun1`). Validated against the old `.con` incident column to <0.05%. One file per (α_ox, α_uv).
  - The big-bump temperature axis has 20 values: log T = 4.0–4.25 in 0.05 steps, 4.3–5.0 in 0.1 steps, 5.0–6.5 in 0.25 steps.
  - Normalised per erg/s of L_bol.
- **Photoionisation settings:**
  - U fixed per grid point;
  - GalacticConcordance abundances (Z☉ = 0.015), no grains, cosmic rays on, CMB off;
  - `iterate to convergence` (CLOUDY's default limit is 10);
  - a 100 K floor and `stop temperature off`, with `nend` 6000 (3000 in the pilot);
  - `set continuum resolution 1.0`, plane parallel.

  Design and axis choices: `cloudy_grid_plan.txt`.

| Run set | Varied axes (besides T) | Models |
|---|---|---|
| BLR core (v_turb 10, Z☉) | log U −2…−0.5 (4), log n_H 8…12 (9), log N_H 22…24 (5) | 3,600 |
| Dense core | log U −2…−0.5 (4), log n_H 9…12 (7), log N_H 23.5…25.5 (5), v_turb 100/400, Z 0.01/0.1/1 Z☉ | 16,800 |
| BLR variations (about U −1.5, n_H 10^10, N_H 10^23) | log U −3/−2.5/0; v_turb 0/30/100; Z 0.1/3/10 Z☉; α_uv −1/−0.1/0.33 | 240 |
| Dense variations (about U −1.5, n_H 10^9.5, N_H 10^24.5, v 400, 0.1 Z☉) | log U −3/−2.5/0; log N_H 22/22.5/23/26; v_turb 10/250/500; α_ox −1/−1.8/−2.2/−3/−5; α_uv −1/−0.1/0.33 | 360 |

## 3. What was run

| Stage | Job | Runs | Changes from production | Outcome |
|---|---|---|---|---|
| Smoke test | 1388268 | 1 | — | converged in 4 iterations; ~70 s fixed start-up per run |
| Timing pilot | 1388277 | 70 | `nend` 3000 | 66 ok; median 14 min, max 2.2 h; 22 core-h |
| Production | 1388297 (40 × 28 cores, ≤ 10 nodes) | 21,000 | — | 15 Oct 04:08 → 16 Oct 10:26; 6,987 core-h (median 16.6 min, max 8.5 h); 20,437 ok |
| Production retry | 1388560 | 525 | — | task 39 died in a SLURM TaskProlog failure before starting; rerun from its own run list |
| fix_test | 1388600 | 64 | omit `turbulence 0`; `set presioniz` (invalid syntax, results void) | v_turb 0: 13/13 ok; the BLR log T 4.25 crash recurred; task 2 (TaskProlog) never started |
| fix_test_2 | 1388604 | 16 | masers off; 30 iterations; n_H nudge | masers: 3/3 ok; 30 iterations: 5/10 converged; dense nudge ok; BLR nudge masered |
| fix_test_3 | 1388616 | 43 | `set pressure ionize 30000` | 26 ok; 17 failed again (cap or 20-failure limit) |
| fix_test_4 | 1388671, 1388672 | 165 | 100 iterations; cap at high T; v_turb 0; nudges | v_turb 0: 7/7; high-T cap: 11/19; 100 iterations: 69 converged, 66 limit cycles, and 2 slow runs (BLR 6/167, dense 5/251) that later joined the limit cycles (§4); nudges ok |
| fix_test_5 | 1388748 | 26 | `failures 100` + `set pressure ionize 300000` + 100 iterations | 22 ok. 4 failed: 3 aborted on 100 convergence failures within iteration 2, and 1 hit the 47 h per-run cap. All 4 are filled (§6). Ended 19 Oct 06:43, after 47 h of its 48 h limit |
| fix_test_6 | 1388756 | 1 | `failures 100` | ok: 24 convergence failures (over the old limit of 20), full column, converged in 4 iterations, 12 min. Smooth between its T neighbours (continuum ≤ 0.003 dex, Hβ/Hα ≤ 0.02 dex; the Lyα pixel lies 0.05 dex beyond the 4.9 value) |
| No-op test | 1388753 | 48 | each fix on 8 healthy runs | done (§5): `failures`, the pressure cap and 100 iterations change nothing; masers off changes every run and made one abort |
| Cycle test | 1388749 | 68 | restarted from production with 30 iterations, continuum saved each iteration | none converged; cycles stationary; amplitudes in `cycle_amplitudes.csv` (§4, limit cycles) |
| Grid build | 1391044 | all 19 grids | — | 26 min; every check passed (§6) |
| Dense core rebuild, final fill | 1391692 | dense core | every filled part clipped | 27 min (22 Oct 19:12–19:39); every check passed (§6) |

Operational lessons:
- A task that dies within seconds and leaves no timing file is a TaskProlog failure. Rerun it from its own run list, never with `--array=N`: `SLURM_ARRAY_TASK_COUNT` becomes 1 and the slice is empty.
- The account's RM QOS caps jobs at 2 days, although the partition allows 7.
- CLOUDY matches keywords only at the start of a word, so `set presioniz` fails to parse. Parse-check new deck commands with a short run first.

## 4. Failure classes (production)

| Class | Runs | Where | Outcome |
|---|---|---|---|
| Not converged in 10 iterations ("soft": full column reached) | 289 | across the grid; worst at dense N_H 10^23.5 | 142 inside the kept T ranges were rerun. 5 converged with 30 iterations, and 69 with 100. The other **68 are limit cycles** (below): kept, flagged, with their amplitude |
| Maser aborts (`RT_line_escape` large negative τ; 1 `escmase`) | 101 | log T ≤ 4.2, 88% Z☉, N_H ≥ 10^24 | 3 inside the kept range (dense, log T 4.2) were fixed with `database H-like masers off` (3/3). The rest are cut |
| NZONE (6000-zone limit) | 77 | low T, high n_H, high U, the largest N_H, Z☉ | all outside the kept T ranges, so cut |
| Pressure–ionisation aborts (`nPres2Ioniz exceeds limPres2Ioniz`) | 64 | dense core only: Z☉, N_H 10^25–25.5, mostly v_turb 400, low U; rising with T | 2 cut (log T < 4.2). Cap 30000: 37 of 62 fixed. Of the 25 leftovers, `failures 100` with cap 300000 (fix_test_5) fixed 21. The last 4 (Z☉, v_turb 400, N_H 10^25.5) are filled along T |
| v_turb = 0 (`turbulence 0` refused by CLOUDY) | 20 | the BLR turbulence variation | the command is omitted, so only thermal broadening applies: 20/20 ok |
| Floating-point crashes | 10 | same regime as the masers | 2 inside the kept ranges (BLR and dense, log T 4.25, log U −0.5, N_H 10^24). Both are deterministic, so their n_H was nudged by +0.001 dex (BLR also with masers off): both ok |
| ConvFail (20-convergence-failure limit) | 2 | dense log T 6.5 (incident 19, index 828); BLR Z = 3 Z☉ variation, log T 4.8 | Both fixed by `failures 100` (BLR: fix_test_6; dense: fix_test_5) |

Also, 181 runs are electron-scattering thick (τ_es ≈ 5–8; N_H 10^25–25.5). They complete at the requested column, so they are accepted and flagged (author, 15 Oct). 39 of the ok ones lie inside the kept dense range.

**Limit cycles.** CLOUDY calls an iteration unconverged when a checked line optical depth (H- or He-like Lyα, He-like subordinate lines) changes by more than 20% from the previous iteration (`conv_itercheck.cpp`).
- The 66 runs still unconverged after 100 iterations fail a check in every iteration, either on one check each time or in a repeating period-3 pattern. That is a limit cycle, not slow convergence. More iterations won't help. Two more that hadn't converged in 100 iterations (BLR 6/167, dense 5/251) joined them: 68 in all, BLR 41 and dense 27.
- Their NIRCam bands moved ≤ 3.9% between 10 and 100 iterations, but those snapshots could all share a phase of the cycle.
- **The cycle test (1388749)** restarted all 68 from production with a continuum save every iteration. None converged, and every one fails the H- or He-like Lyα check. The cycles are stationary: the window mean drifts ≤ 2.1% over iterations 11–30. Iteration 10 reproduces production exactly.
- **Half-amplitude** in the 7 NIRCam wide bands at z = 5–8, outside the Lyα band (the biggest, > 10%, sit at n_H ≥ 10^11 and log T ≤ 4.9):

  | | Models | Median | 90th percentile | Max |
  |---|---|---|---|---|
  | BLR core | 41 | 8.7% | 22% | 25% |
  | Dense core | 27 | 2.8% | 12% | 15% |

- **Decision (author, 18 Oct):** keep the production spectrum, flag it `limit_cycle`, and store the half-amplitude (the largest over the bands, Lyα included) in `smock_cycle_half_amplitude`. Production is within the cycle (BLR median 5.9% from the cycle mean, dense 1.9%), and the per-iteration saves hold only the continuum, so a cycle-mean spectrum would not match its lines.

## 5. Fixes and how they were checked
Each fix was applied only to runs that had failed, so production's successful runs are untouched. **No-op test (1388753, 17 Oct; `cloudy_runs/reruns/noop_test/compare_noop.txt`):** an unchanged rerun reproduces production exactly (8/8), so the comparison is exact. The runs: dense 11/828 (pressure–ionisation corner, log T 4.8), dense 18/741 (log T 6.25), dense 4/320 (log T 4.2, maser regime), BLR 5/155 (log T 4.25, maser regime), dense 14/0 (typical), BLR 13/25 and 14/91 (9 iterations), and one BLR variation run. Whether each fix could become a recipe default depends on whether it changes a healthy run. The no-op test (1388753) checks this directly: 8 healthy production runs from the risky regimes, each rerun with each fix alone and with all of them together, compared exactly with an unchanged rerun.

| Fix | Acts on | Evidence it helps | No-op for healthy runs? |
|---|---|---|---|
| Omit `turbulence 0` | v_turb = 0 decks only | 20/20 ok | yes, by construction (it only touches those decks) |
| `set pressure ionize 30000` → 300000 | runs that reach the cap | 37/62 at 30000 | **yes**: identical, 8/8 |
| `failures 100` | runs that reach 20 convergence failures | 2/2: BLR Z 3 Z☉ variation (fix_test_6) and dense 19/828 (fix_test_5) | **yes**: identical, 8/8 |
| `iterate to convergence 100` | runs not converged in 10 iterations | 69 of 137 converged; 67 of them moved up to 28% outside Lyα, so their 10-iteration spectra are replaced | **yes**: identical, 8/8, including two runs that needed 9 iterations |
| `database H-like masers off` | H-like masing levels | 3/3 dense log T 4.2 maser aborts ok; the BLR nudge ok | **no**: it changes all 8 (continuum up to 0.8% of its peak, lines up to 0.1%, zone counts shift), and it made the healthy dense 11/828 abort on the pressure–ionisation cap. Rerun-only |
| n_H nudge (±0.001 dex) | deterministic crashes only | 10^11.501 and 10^11.499 agree to ≤ 0.001 dex in every band; leave-one-out fits along n_H show nudged points as smooth as their neighbours | a manual stage, not a default |

## 6. Building the synthesizer grids
`build_agn_grids.py` (dry run by default) builds each grid in its own directory, `cloudy_runs/build/<grid>/`:

1. **Choose each model's source.** That's the production run, unless a later rerun of it succeeded. Success is `parse_cloudy_runs.py`'s `ok` flag. Each run's location is read from its stage's run list. The dense 10^11.499 nudge is ignored (it was a cross-check).
2. **Cut the axes before syncretize sees them:** a copy of the incident grid with only the kept T values, plus trimmed YAMLs. syncretize only ever sees the intended grid (author).
3. **Link** the chosen outputs into a final view of relative symlinks, re-indexed by parameter value. `MANIFEST.tsv` records every model's source and flags.
4. **Run syncretize's `create_synthesizer_grid.py`** unchanged.
5. **Fill the holes** along T (below).
6. **Finish** each grid:
   - restore `WeightVariable`;
   - add `smock_flags` (1 rerun, 2 nudged, 4 thick_es, 8 hole, 16 masers_off, 32 limit_cycle, 64 filled), `smock_cycle_half_amplitude`, `smock_fill_error`, `smock_build_note` and `smock_build`;
   - write `MANIFEST.tsv`;
   - check the normalisation against the incident grid.

**Filling the last holes** (author, 18 Oct). synthesizer 1.3.1 reads `failures` only as metadata, so its interpolation would use a hole's zero spectrum. Each of the 4 is filled by a cubic in log T through T−2, T−1, T+1, T+2. 13/824 and 14/824 are adjacent along T, so each uses the other's gap: log T 4.8, 4.9, 5.5 and 5.75.
- **Choosing the scheme:** a leave-one-out test over the 444 good models of the hole slab (18 Oct; archive `2026-10-18_cloudy_hole_fill_loo/`). Each model was hidden and predicted from its neighbours. Near the holes, in the bright NIRCam bands of the net transmitted spectrum:

  | Scheme | Median | Max |
  |---|---|---|
  | **Cubic in log T** | **1.2%** | **4.7%** |
  | Linear in log T | 4.4% | 12% |
  | Linear in n_H | 12% | 41% |
  | Linear in U | 46% | 54% |

  Mixing axes makes it worse. Lines: cubic median 2.6%, 90th percentile 10%.
- **Each stored spectrum** (21 Oct; archive `2026-10-21_hole_fill_components/`). This test settled three details:
  - zeros are propagated instead of floored in log (a floored zero with a negative cubic weight put one gap fill off by 113%);
  - every part is clipped to its two bracketing neighbours (the final fill, author 22 Oct). At first only the transmitted was clipped (the first fill), because it overshoots where heavily absorbed. The built dense core then showed the line continuum overshooting too, at 1658 Å in the adjacent pair: 5–9× above every T neighbour, because the real runs rise steeply there. The final fill costs a little in the leave-one-out (plain cubic, rebuilt total: median 0.9% → 1.5%, max unchanged) and removes the catastrophic outliers in absorbed bands. It matters most when the covering fraction is ~1, where nothing escapes to swamp the error. The line luminosities are clipped the same way, to stay consistent with the line continuum;
  - the nebular is filled as nebular continuum + linecont, so the continuum stays positive.

  Result, the rebuilt total (transmitted + nebular): plain cubic max 4.7%, gap cubics max 7.3% / 9.0%.
- **Per ionising photon:** CLOUDY's incident Q is constant along T at fixed U and n_H (1.8×10⁻⁴ dex), so each quantity is filled divided by the grid's Q and then rescaled. Each filled model gets its own syncretize `normalisation`. The incident is copied from a model at the same T, and the line continua are recomputed from the filled spectra as syncretize does.
- **Flags and error:** `filled` replaces `hole`, `failures` is set to 0, and `smock_fill_error` holds a leave-one-out error measured on the filled model's converged neighbours.
- **Not reliable:** the heavily absorbed rest-UV bands (transmitted < 1% of the incident). There, every scheme is off by tens of percent or more. They matter only when the covering fraction is close to 1.

**Two syncretize problems, both handled:**
- **Normalisation.** syncretize loads the incident grid through synthesizer, which converts AGN grids from per erg/s to per L☉ on every load. An unmodified build therefore comes out per L☉ and is converted again when used: 3.83×10³³ too bright.
  - The builder gives syncretize an incident copy with `WeightVariable "None"`, then sets `bolometric_luminosities` on the output.
  - Every build checks that the loaded grid's incident spectrum equals the incident grid's own. First builds: within 1.3×10⁻⁴.
  - **If AGN normalisation ever looks wrong, check this first.**
- **Failure test.** syncretize treats a model as failed only if its line list is missing or tiny, so it accepts runs cut short at NZONE (it flagged 1 of the pilot's 4 failures). The final view therefore contains only `ok` runs. A model with no ok run gets only its `.out`, so syncretize records it as a hole (zero spectra, `failures = 1`), and the builder then fills it.

**And one in the builder itself (fixed 21 Oct):** syncretize's model order is not C order. With several photoionisation axes, unravelling the model number in C order puts the per-model datasets at the wrong grid points. The builder now takes each model's grid index from syncretize's own index list. Grids with one photoionisation axis were never affected, and the cores had not been built yet.

**Final grids** (job 1391044, 21 Oct):

| Grid | Shape (T, then sorted photoionisation axes) | Models | From reruns | Limit cycles | Filled |
|---|---|---|---|---|---|
| BLR core | (15, 9, 4, 5) | 2,700 | 21 | 41 | 0 |
| Dense core | (16, 7, 4, 3, 5, 2) | 13,440 | 117 | 27 | 4 |
| 17 variation grids | (15 or 16[, 3 or 4]) | 15–64 each | 16 (BLR v_turb 15, BLR Z 1) | 0 | 0 |

Build checks (`archive/2026-10-21_core_grid_checks/`), all passed on all 19 grids:
- normalisation within 1.4×10⁻⁴ of the incident grid;
- 1,100 sampled models (every rerun and limit cycle included) equal their own CLOUDY runs to 2×10⁻¹⁶;
- `smock_flags` equals MANIFEST everywhere, and the flag counts are as planned;
- no failures and no empty models;
- fill leave-one-out errors 1.5–3.4% (first fill).

**The dense core rebuilt with the final fill** (job 1391692; `archive/2026-10-23_gas_enshrouded_core_v5_checks/`, 23 Oct): the checks above pass again (normalisation within 1.3×10⁻⁴; 340 sampled models equal their runs to 2×10⁻¹⁶; flags as MANIFEST; no failures, no empty models; leave-one-out 1.5–3.5%). Against the first fill:
- the 13,436 models that aren't filled are identical in every dataset, so only the fill changed;
- every filled part (transmitted, nebular continuum, linecont, line luminosities) now lies between its two T neighbours at every wavelength; the first fill's linecont reached 1.45× the brighter neighbour;
- in the bands at z = 6 only F115W moves noticeably, in the adjacent filled pair (13/824 and 14/824): 0.31× and 0.62× of the first fill. Every other band changes by < 1%.
- The parts are clipped separately, so their sum can still sit above both neighbours: F115W at z = 6 in that pair is 0.06 and 0.19 dex above the brighter one (first fill: up to 0.51 dex). The band is line-dominated and heavily absorbed there.

The shipped file (`..._cloudy-gas_enshrouded_core.hdf5`) is installed in `synthesizer_data/grids/agn/photoionised/agn-cloudy/` with the other 18 grids.

## 7. Open items
- ~~Parse fix_test_5, fix_test_6 and `fix_test_4/runs_slow.txt`; decide the remaining holes and the limit-cycle runs~~ done (§4, §6).
- ~~The no-op test~~ done. `failures 100`, `set pressure ionize 300000` and `iterate to convergence 100` are now defaults in `make_cloudy_inputs.py` (author, 17 Oct). Their text matches the tested decks exactly (8/8). Masers off stays rerun-only, and its 4 grid points carry the `masers_off` flag.
- ~~Build every grid and check it~~ done (21 Oct). ~~Rebuild the dense core with the final fill and re-check it, then install the grids~~ done (22–23 Oct).
- ~~Check that a sample of grid points reproduces from scratch with the recipe defaults, including one from each failure class~~ dropped (author, 24 Oct 2026): no CLOUDY reruns for reproduction. The grid may change with the author's work on the CLOUDY AGN incident shape, and the whole grid could be rerun in about a week of jobs if needed.
- Reproduction instead (author, 24 Oct 2026): after the raw runs move to their long-term place, rebuild a variation grid and a core grid from the kept copy and compare them with the installed files. **Done (job 1392405):** the gas-enshrouded core and the standard-BLR turbulence variation are identical to the installed grids apart from the build stamps (`archive/2026-10-24_agn_grid_reproduction/`). The evidence so far: the build checks (§6), in which sampled models equal their own CLOUDY runs to 2×10⁻¹⁶, and the install, identical to the builds.
- ~~The pilot band comparison (§8)~~ done: the dense v_turb 500 variation is kept (author, 21 Oct).

## 8. Pilot band comparisons
The question was whether the dense v_turb 500 variation is worth keeping, and what N_H 10^26.5 adds. The pilot runs at the dense fiducial (log T 5, log U −1.5, n_H 10^9.5, 0.1 Z☉) were compared in the NIRCam wide bands at z = 5–8 (net transmitted, photon-weighted; archive `2026-10-21_pilot_band_comparison/`):
- **v_turb 500 vs 400:** median 3.5%, max 11% (F277W at z = 7), mostly fainter at 500. That is a small but real change. The variation is already run and built, so it is kept (author, 21 Oct 2026).
- **N_H 10^26.5 vs 10^26:** the bands drop by 39–100% (the rest-UV bands by ~100%). 10^26 vs 10^24.5 is similar. At these columns the slab is very Compton-thick, so each half-dex in N_H is a large change, not a refinement. The grids stop at 10^25.5 (cores) and 10^26 (dense variation).

## 9. Where the evidence is
- Parse CSVs:
  - `cloudy_runs/pilot_runs.csv`;
  - `cloudy_runs/production_runs.csv`;
  - `cloudy_runs/reruns/<stage>/runs*.csv`.
- Failure-map scripts (16 Oct): kitchen archive `2026-10-16_cloudy_failure_map/`.
- Convergence:
  - `fix_test_2/compare_iter30_final.txt`;
  - `fix_test_4/compare_iter100.txt`;
  - `fix_test_4/reason_sequences.tsv`;
  - `cycle_test/analyse_cycles.py` → `cycle_amplitudes.csv`.
- Hole fill: archive `2026-10-18_cloudy_hole_fill_loo/` and `2026-10-21_hole_fill_components/`.
- Nudges: `fix_test_4/nudge/nudge_checks.txt`.
- Builds: `cloudy_runs/build/<grid>/MANIFEST.tsv` and `syncretize.log`.

# walk devel vs modified (speed-up-slim): benchmark report

Generated 2026-11-10 11:32 by `bench_walk_devel_vs_modified.R`. Run time 103.6 min.

## Setup

- **devel:** https://github.com/dsjohnson/walk.git @ `devel` (b39f00b).
- **modified:** https://github.com/bmcclintock/walk.git @ speed-up-slim (b48d30d).
- **Versions loaded:** walk 0.0.1.9005 (devel) and 0.0.1.9005 (modified; OpenMP TRUE). Each ran in its own R process.
- **Simulated data:** langevinSSM 0.0.2, `simLangevin` (underdamped, beta = (-4, 6, 5, -0.1), sigma = 5, gamma = 0.5) on `exampleCovs`.
  - 4 animals, 1200 fixes, subsampled to irregular times (median dt 0.070).
  - Argos-type error ellipses (median smaj 2.38).
- **Grids:** `exampleCovs` aggregated by 4, 2, 1, i.e. 10,000, 40,000, 160,000 cells.
- **Model:** movement `~cov1 + cov2 + cov3`, log links, `norm = FALSE`. −2ll evaluated at 3 parameter vectors.
- **Timing:** median of 3 run(s) after one warm-up. Machine: Darwin, 12 cores, R 4.5.2.

## Checks

- **L, devel vs modified `buffer` (abs tol 1e-10): PASS.**
  - grid4: 15097 vs 15097 non-zeros, 0/0 cells only in one version, max abs diff 4.44e-16.
  - grid2: 42657 vs 42657 non-zeros, 0/0 cells only in one version, max abs diff 9.16e-16.
  - grid1: 138651 vs 138651 non-zeros, 0/0 cells only in one version, max abs diff 3.15e-15.
- **−2ll, devel vs modified (`buffer` L, and the modified likelihood on devel's own L), rel tol 1e-08: PASS.**
  - Max relative difference 8.64e-10 over 36 comparisons.
- **Joint call vs sum of per-animal values (modified):** max abs diff 0.
- **Full fit of animal 1** (10,000 cells, nlminb, same start):
  - −2ll: devel 897.3738143, modified 897.3738143 (difference 2.96e-12).
  - Max abs parameter difference 9.15e-07; convergence codes 0/0.
- **`disc` vs `buffer` (not expected to be identical):**
  - Extra cells in `disc` that the polygonised buffer misses: 10,000 cells: 15; 40,000 cells: 68; 160,000 cells: 277.
  - Its −2ll differs from devel by 1.7e-07 to 1.1e-06 (relative).

## Speed

Speed-up = devel time / modified time. All times are for all animals together.

### `proc_Langevin`

| cells | devel | buffer, per animal, 1 thr | buffer, one call, 2 thr | disc, per animal, 1 thr | disc, one call, 2 thr |
|---|---|---|---|---|---|
| 10,000 | 2.40 s | 0.130 s (18×) | 0.110 s (22×) | 0.019 s (126×) | 0.008 s (300×) |
| 40,000 | 7.57 s | 0.225 s (34×) | 0.187 s (40×) | 0.028 s (270×) | 0.013 s (582×) |
| 160,000 | 26.34 s | 0.614 s (43×) | 0.551 s (48×) | 0.067 s (393×) | 0.032 s (823×) |

### −2 log-likelihood (one evaluation; range over parameter vectors)

| cells | devel | per animal, devel's L | per animal, buffer | joint, 1 thr | joint, 2 thr |
|---|---|---|---|---|---|
| 10,000 | 0.71–0.99 s | 0.031–0.046 s (22–23×) | 0.031–0.046 s (22–23×) | 0.011–0.025 s (40–65×) | 0.009–0.017 s (58–79×) |
| 40,000 | 2.33–4.22 s | 0.127–0.218 s (18–19×) | 0.131–0.212 s (18–20×) | 0.040–0.131 s (32–58×) | 0.037–0.086 s (49–67×) |
| 160,000 | 9.74–22.05 s | 0.561–2.634 s (8–20×) | 0.545–2.651 s (8–19×) | 0.167–2.290 s (10–58×) | 0.147–1.375 s (16–66×) |

**Full fit, animal 1:** devel 15.1 s → modified 0.16 s (94×).

## Notes

- **Devel has no joint likelihood,** so the devel figure is `ctmc_n2ll` per animal, summed.
- **"per animal, devel's L"** runs the modified likelihood on exactly the inputs devel used. It isolates the likelihood code from the `proc_Langevin` step.
- **`disc`** uses every cell that touches the exact error circle, while `buffer` uses `terra::buffer` + `terra::cells(touches = TRUE)`. Only `buffer` is expected to match devel.
- **Files:** `walkbench_summary.txt` (all tables), `walkbench_*.csv`, `walkbench_results.rds`.

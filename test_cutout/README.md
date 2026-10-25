# Recipe: ASTRID test cutout

**Status: done** (24 Oct 2026). The product was made and validated on 15 Oct 2026, and the
recipe reproduces it byte for byte (§5).

## 1. Product
A small BigFile with a PIG's block layout, holding 5 FoF groups (264, 557, 784, 9732, 736267;
194 subhalos) from ASTRID PIG_047 (z = 6): gas, star and BH fields, no dark matter, 17 MB. It
goes to `smock/test_data/ASTRID/PIG_047_subfind/` (`smock.config.TEST_CUTOUT_PIG_DIR`),
with provenance in its `README.txt`. ASTRID is public, so the cutout can be shared
(data-sharing terms checked by the author). Published (24 Oct 2026) in the Google Drive folder `smock_data` as one
file, `test_data.zip`, made from smock's folder with `cd smock && zip -r -X test_data.zip test_data`
(one file instead of 228, so command-line downloads stay under Google's limit; users unzip it inside
their smock folder). No raw data kept.

## 2. Requirements
- The `sharable_mock_obs` env with smock installed.
- Read access to the ASTRID PIG (`smock.config.EXTERNAL["astrid_pig_dir"]`; on Vera
  `/hildafs/datasets/Asterix/PIG2/`, read-only, owned by others).
- Compute: seconds to minutes on the login node.

## 3. Steps
From the top of this repository:

```bash
python test_cutout/make_test_cutout.py              # dry run
python test_cutout/make_test_cutout.py --apply      # write, then validate
python test_cutout/make_test_cutout.py --validate   # re-check an existing cutout
```

## 4. Problems solved
- **Which groups** (author, 15 Oct 2026): one subhalo per category, from a survey of PIG_047
  (`archive/2026-10-15_test_group_survey/`). 151,889 subhalos pass the run_5 cut. At z = 6,
  ASTRID BHs are light and overmassive bright BHs are absent, so "LRD-like" is ranked by
  L_bol/M* among compact (r½* < 1.05 pkpc) systems.
  - typical star-forming galaxy with an accreting BH: 303604;
  - LRD-like: 54437;
  - several subhalos: 24900 (two comparable galaxies 1.8″ apart);
  - small galaxy: 4283857 (102 stars);
  - massive group, for chunking: 42595 (54k group stars).
- **Finding a subhalo's group:** offset containment.
  `SubhaloGroupNr` is not always the FoF row.
- **Offsets** are renumbered so the cutout reads like a full PIG; the original rows are kept
  in `*/OriginalIndex`. Provenance goes in `README.txt`, because string header attributes
  don't read back under bigfile 0.1.52 + numpy 2.

## 5. Verification
`--validate` reads every subhalo and group from the full PIG and the cutout and asserts the
arrays are identical, then checks units.
- **Shipped cutout, 24 Oct 2026:** 1,644 particle arrays (5 groups, 194 subhalos) and every
  catalogue block not renumbered are identical to the PIG; Σ particle mass = `SubhaloMassType`
  to ≤ 5.7×10⁻⁸ (gas, stars, BHs); `SubhaloHalfmassRadType` / particle r½* = 0.9999 (0.9977–1.0046,
  16 subhalos with ≥ 100 stars). The first validation (15 Oct) counted 1,539 arrays: that was
  before `TotalMassReturned` was added to the star fields the same day (one array more for each
  of the 105 objects with stars) and the cutout rewritten.
- **Reproduction, 24 Oct 2026:** `--apply --out-dir runs/test_cutout_repro` (10 s on the
  login node) rewrote the cutout, which then validated as above.
  All 228 files are byte-identical to the shipped cutout's, except the first line of
  `README.txt` (the date, and the script's path before the repo split).

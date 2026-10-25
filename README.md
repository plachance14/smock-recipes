# smock-recipes

How the data that [smock](https://github.com/plachance14/smock) downloads were made. smock users
don't need anything here: the products below come ready-made in smock's
[Google Drive folder](https://drive.google.com/drive/folders/1aruCXByevDyer_LeufUUbo_knz311sDV?usp=sharing).
This repository is for rebuilding them, changing them, or making new ones the same way.

There is one folder per product. Each README has the same sections: the product, what you need,
the steps, the problems met along the way and how the recipe deals with them, and how the result
was checked. The "problems solved" sections are the reason the recipes exist: please read them
before simplifying a recipe, or an old problem may come back.

| Product | Recipe | Where it goes | Checked |
|---|---|---|---|
| CLOUDY AGN grids (standard-BLR and gas-enshrouded cores, 17 variations) | [`agn_grid/`](agn_grid/) | `SYNTHESIZER_DIR/grids/agn/` | two grids rebuilt from the kept CLOUDY runs are identical to the shipped ones |
| NIRCam PSFs at 0.03″, 7 filters | [`psfs/`](psfs/) | `SYNTHESIZER_DIR/instrument_cache/psfs/` | reproduced bit for bit |
| ASTRID PIG_047 test cutout | [`test_cutout/`](test_cutout/) | `smock/test_data/ASTRID/PIG_047_subfind/` | reproduced byte for byte |

`SYNTHESIZER_DIR` is synthesizer's own data folder, which smock shares (see smock's `SETUP.txt`
and `smock/config.py`). The Drive folder also holds the BPASS 2.2.1 grid the walkthrough uses
(synthesizer's published grid) and the seven NIRCam filter curves (synthesizer's copies from
the SVO Filter Profile Service).

## Running the recipes
- Run every command from the top of this repository. The AGN grid's work directory
  (`cloudy_runs/`) and the PSF recipe's STPSF tools (`external/`) go there, and are not committed.
- The recipes import smock, so install it first (smock's `SETUP.txt`, then
  `pip install --no-deps -e /path/to/smock`).
- The SLURM scripts carry the account and partition used on PSC's Vera; change them for your
  cluster.

## References to the build record
The READMEs, comments and the AGN grid report mention dated folders (e.g.
`archive/2026-10-21_hole_fill_components/`) and the ROADMAP. These are in my private record of
the build, which holds the tests and evidence behind each decision; they are named here so each
decision can be traced. If you need one of them, ask me (Patrick LaChance, p.lachance@soton.ac.uk) or open an
[issue](https://github.com/plachance14/smock-recipes/issues).

## Licence
MIT; see [LICENSE](LICENSE).

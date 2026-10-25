"""Install the built AGN grids into the synthesizer data directory (dry run by default).

Copies each build's shipped grid (cloudy_runs/build/<grid>/grids/*_cloudy-{standard_blr,
gas_enshrouded}_*.hdf5; build_agn_grids.py step 7) to

    <GRID_DIR>/agn/photoionised/agn-cloudy/

and the incident grids they were built from (cloudy_runs/incident_grids/agn-cloudy_T20_*.hdf5,
the full 20-point T axis; not the pilot's) to

    <GRID_DIR>/agn/incident/agn-cloudy/

where GRID_DIR comes from smock.config: synthesizer's own grids folder (SYNTHESIZER_GRID_DIR,
else $SYNTHESIZER_DIR/grids).
synthesizer needs only the photoionised grids; the incident grids are there to rebuild or
extend them. An existing file is never replaced unless --overwrite is given.

Usage (from the top of smock-recipes, which holds cloudy_runs/):
    python agn_grid/install_agn_grids.py                 # dry run
    python agn_grid/install_agn_grids.py --apply [--grids gas_enshrouded_core]
"""
import argparse
import shutil
import sys
from pathlib import Path

from smock import config

PHOTOIONISED = Path("agn/photoionised/agn-cloudy")
INCIDENT = Path("agn/incident/agn-cloudy")
SHIPPED_MARKERS = ("_cloudy-standard_blr_", "_cloudy-gas_enshrouded_")


def files_to_install(work, grid_filter):
    """(source, destination relative to the grid directory) for every file to install."""
    pairs = []
    for src in sorted((work / "build").glob("*/grids/*.hdf5")):
        if any(m in src.name for m in SHIPPED_MARKERS) and (not grid_filter or any(g in src.name for g in grid_filter)):
            pairs.append((src, PHOTOIONISED / src.name))
    if not grid_filter:
        pairs += [(src, INCIDENT / src.name) for src in sorted((work / "incident_grids").glob("agn-cloudy_T20_*.hdf5"))]
    return pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--work-dir", default="cloudy_runs")
    parser.add_argument("--grids", nargs="+", help="substrings of shipped grid names (default: all, plus incident grids)")
    parser.add_argument("--overwrite", action="store_true", help="replace files already installed")
    parser.add_argument("--apply", action="store_true", help="copy (default: dry run)")
    args = parser.parse_args()

    if not config.GRID_DIR.is_dir():
        sys.exit(f"no grids folder at {config.GRID_DIR}: set SYNTHESIZER_DIR (or SYNTHESIZER_GRID_DIR)")
    pairs = files_to_install(Path(args.work_dir), args.grids)
    for src, rel in pairs:
        dst = config.GRID_DIR / rel
        action = "replace" if dst.exists() and args.overwrite else "exists, skip" if dst.exists() else "copy"
        print(f"{action:13s} {src.stat().st_size / 1e6:9.1f} MB  {rel}")
        if args.apply and action != "exists, skip":
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
    print(f"{len(pairs)} files -> {config.GRID_DIR}")
    if not args.apply:
        print("Dry run: nothing copied. Re-run with --apply.")


if __name__ == "__main__":
    main()

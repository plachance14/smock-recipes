"""Generate CLOUDY input files for one run set with syncretize.

A run set is a list of (incident grid, base parameter YAML, optional extra YAML).
syncretize runs every incident point against every combination of the list-valued
parameters, so one-at-a-time variations are a fiducial base plus an extra YAML
that overrides a single axis (``--cloudy-paramfile-extra``); alpha_ox and alpha_uv
variations swap the incident grid instead (see make_incident_grid.py).

syncretize's create_cloudy_input_grid.py reads ``params/<file>`` and copies
``linelist*.dat`` relative to its working directory, so this script stages a
work directory (default ``cloudy_runs/``, git-ignored) with those in place.

It also writes ``<work>/<set>_runs.txt`` (one model directory and index per line),
which the SLURM scripts read.

Usage (sharable_mock_obs env, from the top of smock-recipes, which holds cloudy_runs/; see README.md):
    python agn_grid/make_cloudy_inputs.py pilot            # dry run: list the steps
    python agn_grid/make_cloudy_inputs.py pilot --apply
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import synthesizer_grids

HERE = Path(__file__).resolve().parent
PARAMS = HERE / "params"
SYNCRETIZE_CLOUDY = Path(synthesizer_grids.__file__).parent / "cloudy"

FID = "agn-cloudy_T20_aox-1.40_auv-0.50"
T20_AOX = [f"agn-cloudy_T20_aox{a}_auv-0.50" for a in ("-1.00", "-1.80", "-2.20", "-3.00", "-5.00")]
T20_AUV = [f"agn-cloudy_T20_aox-1.40_auv{a}" for a in ("-1.00", "-0.10", "0.33")]

# (incident grid, base YAML, extra YAML or None)
RUN_SETS = {
    "pilot": [
        ("agn-cloudy_pilotTcorners_aox-1.40_auv-0.50", "pilot_corners.yaml", None),  # 64
        ("agn-cloudy_pilotT5_aox-1.40_auv-0.50", "blr_fiducial.yaml", None),  # 1
        ("agn-cloudy_pilotT5_aox-1.40_auv-0.50", "dense_fiducial.yaml", None),  # 1 (also v_turb 400)
        ("agn-cloudy_pilotT5_aox-5.00_auv-0.50", "dense_fiducial.yaml", None),  # 1
        ("agn-cloudy_pilotT5_aox-1.40_auv-0.50", "dense_fiducial.yaml", "pilot/turbulence_500.yaml"),  # 1
        ("agn-cloudy_pilotT5_aox-1.40_auv-0.50", "dense_fiducial.yaml", "pilot/column_26.yaml"),  # 2
    ],
    "blr_core": [(FID, "blr_core.yaml", None)],
    "dense_core": [(FID, "dense_core.yaml", None)],
    "blr_variations": (
        [(FID, "blr_fiducial.yaml", f"blr_variations/{v}.yaml")
         for v in ("ionisation_parameter", "turbulence", "metallicities")]
        + [(g, "blr_fiducial.yaml", None) for g in T20_AUV]
    ),
    "dense_variations": (
        [(FID, "dense_fiducial.yaml", f"dense_variations/{v}.yaml")
         for v in ("ionisation_parameter", "stop_column_density", "turbulence")]
        + [(g, "dense_fiducial.yaml", None) for g in T20_AOX + T20_AUV]
    ),
}


def stage_work_dir(work):
    """Put params/ and the line lists where syncretize expects them."""
    work.mkdir(parents=True, exist_ok=True)
    link = work / "params"
    if not link.exists():
        link.symlink_to(PARAMS)
    for name in ("linelist.dat", "linelist-standard.dat"):
        shutil.copyfile(SYNCRETIZE_CLOUDY / name, work / name)


def syncretize_grid_name(incident, base, extra):
    """The output directory name create_cloudy_input_grid.py will use."""
    name = f"{incident}_cloudy-{base.removesuffix('.yaml')}"
    if extra:
        name += "-" + Path(extra).name.removesuffix(".yaml")
    return name


ZERO_TURBULENCE = re.compile(r"^turbulence\s+0+(\.0*)?\s*\n", re.MULTILINE)


def omit_zero_turbulence(text):
    """Drop a "turbulence 0.0" command from a deck.

    CLOUDY c23.01 refuses a zero turbulent velocity ("the turbulent velocity needs to
    be > 0", parse_commands.cpp), so the whole v_turb = 0 slice of the first production
    run failed. Without the command CLOUDY uses thermal broadening only, which is what
    v_turb = 0 means (author, 16 Oct 2026). The .yaml keeps turbulence 0 as the axis value.
    """
    return ZERO_TURBULENCE.sub("", text)


ROBUSTNESS_LIMITS = ("failures 100", "set pressure ionize 300000")


def add_robustness_limits(text):
    """Raise three CLOUDY limits that production runs hit (cloudy_grid_report.md §5).

    "iterate to convergence" -> "... 100" (default limit 10), plus "failures 100"
    (convergence-failure abort, default 20) and "set pressure ionize 300000" (pressure-
    ionisation calls per zone, default 3000). Each acts only on a run that reaches its
    limit: in the no-op test (17 Oct 2026, job 1388753) 8 healthy production runs came
    out byte-identical with each one. Production (15 Oct) was generated before these
    became defaults; its failures got them in reruns. "database H-like masers off" is
    deliberately not here: it changes healthy runs and is for maser aborts only.
    """
    text, n = re.subn(r"^iterate to convergence\s*$", "iterate to convergence 100", text,
                      flags=re.MULTILINE)
    if n != 1:
        sys.exit("expected exactly one bare 'iterate to convergence' line")
    return text + "".join(f"{cmd}\n" for cmd in ROBUSTNESS_LIMITS)


def add_performance_save(in_file):
    """Append CLOUDY's per-zone timing output (pilot only) to an .in file."""
    text = in_file.read_text()
    if "save performance" not in text:
        in_file.write_text(text + f'save performance "{in_file.stem}.per"\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("run_set", choices=RUN_SETS)
    parser.add_argument("--incident-dir", default="cloudy_runs/incident_grids")
    parser.add_argument("--work-dir", default="cloudy_runs")
    parser.add_argument("--apply", action="store_true", help="generate (default: dry run)")
    parser.add_argument("--save-performance", action=argparse.BooleanOptionalAction,
                        default=None, help="add 'save performance' (default: on for pilot)")
    args = parser.parse_args()

    work = Path(args.work_dir).resolve()
    incident_dir = Path(args.incident_dir).resolve()
    out_dir = work / "grids"
    save_perf = args.run_set == "pilot" if args.save_performance is None else args.save_performance

    steps = RUN_SETS[args.run_set]
    for incident, base, extra in steps:
        name = syncretize_grid_name(incident, base, extra)
        print(f"{incident} x {base}" + (f" + {extra}" if extra else "") + f"\n    -> {out_dir / name}")
        if (out_dir / name).exists():
            sys.exit(f"{out_dir / name} exists; remove it first (inputs are not overwritten).")
    if not args.apply:
        print("Dry run: nothing generated. Re-run with --apply.")
        return

    stage_work_dir(work)
    runs = []
    for incident, base, extra in steps:
        cmd = [sys.executable, "-m", "synthesizer_grids.cloudy.create_cloudy_input_grid",
               "--grid-dir", str(incident_dir), "--incident-grid", incident,
               "--cloudy-output-dir", str(out_dir),
               # without ".yaml": syncretize builds the grid name before adding it
               "--cloudy-paramfile", base.removesuffix(".yaml")]
        if extra:
            cmd += ["--cloudy-paramfile-extra", extra.removesuffix(".yaml")]
        subprocess.run(cmd, cwd=work, check=True, stdout=subprocess.DEVNULL)

        grid = out_dir / syncretize_grid_name(incident, base, extra)
        for in_file in sorted(grid.glob("*/*.in")):
            in_file.write_text(add_robustness_limits(omit_zero_turbulence(in_file.read_text())))
            if save_perf:
                add_performance_save(in_file)
            runs.append(f"{in_file.parent} {in_file.stem}")

    run_list = work / f"{args.run_set}_runs.txt"
    run_list.write_text("\n".join(runs) + "\n")
    print(f"{len(runs)} CLOUDY inputs; run list: {run_list}")


if __name__ == "__main__":
    main()

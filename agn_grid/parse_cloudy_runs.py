"""Summarise CLOUDY runs (runtime, zones, iterations, stop reason) into a CSV.

Reads a run list written by make_cloudy_inputs.py ("<model_dir> <index>" per line),
the parameters syncretize saved beside each input (<index>.yaml), and each run's
.out. Optionally merges the wall times and exit codes the SLURM job recorded.

The .out lines parsed are CLOUDY c23.01's (prt_comment.cpp, maincl.cpp):
    "Calculation stopped because <reason> Iteration i of n"
    " Cloudy ends: <z> zones, <i> iterations ... ExecTime(s) <t>"
    " C-Iterate to convergence did not converge in <i> iterations."  (default limit 10)
    "Cloudy exited OK" (replaced by "something went wrong" if there were warnings)
    "ABORT DISASTER" (aborted, e.g. a convergence failure)
    " W-The model is optically thick to electron scattering" (tau_es > ~5)

A run is flagged ``ok`` only if it ended without aborting, reached the requested
column, did not hit nend, and its last iteration converged. Warnings alone do not
fail a run: the electron-scattering-thick warning is the physical regime at
N_H ~ 10^26, not a numerical failure (author, 15 Oct 2026), so it gets its own
``thick_es`` column for masking later. syncretize's own check only tests that the
line-list file exists, so this flag is what to use before building a grid.

Usage:
    python agn_grid/parse_cloudy_runs.py cloudy_runs/pilot_runs.txt \
        --timing cloudy_runs/pilot_timing_*.tsv --csv cloudy_runs/pilot_runs.csv
"""

import argparse
import csv
import re
from pathlib import Path

import numpy as np
import yaml

STOP_RE = re.compile(r"Calculation stopped because (.+?) Iteration\s+(\d+) of\s+(\d+)")
END_RE = re.compile(r"Cloudy ends:\s*(\d+) zones?, (\d+) iterations?.*ExecTime\(s\)\s*([\d.]+)")
PARAM_KEYS = ["big_bump_temperature", "ionisation_parameter", "hydrogen_density",
              "stop_column_density", "turbulence", "metallicities"]


def parse_out(path):
    """Pull runtime, zones, iterations, stop reason and warnings from one .out."""
    info = {"exited_ok": False, "ended": False, "aborted": False, "stop_reason": "",
            "converged": True, "zones": None, "iterations": None, "exec_s": None,
            "n_warnings": 0, "thick_es": False}
    if not path.exists():
        info["stop_reason"] = "no .out"
        return info
    for line in path.read_text(errors="replace").splitlines():
        if m := STOP_RE.search(line):  # keep the last iteration's
            info["stop_reason"] = m.group(1).strip()
        elif m := END_RE.search(line):
            info["zones"], info["iterations"] = int(m.group(1)), int(m.group(2))
            info["exec_s"] = float(m.group(3))
            info["ended"] = True
        elif "Iterate to convergence did not converge" in line:
            info["converged"] = False
        elif "Cloudy exited OK" in line:
            info["exited_ok"] = True
        elif "ABORT DISASTER" in line:
            info["aborted"] = True
        elif line.lstrip().startswith("W-"):
            info["n_warnings"] += 1
            if "optically thick to electron scattering" in line:
                info["thick_es"] = True
    return info


def model_params(model_dir, index):
    """The run's parameters from syncretize's <index>.yaml, as log10 where natural."""
    p = yaml.safe_load((model_dir / f"{index}.yaml").read_text())
    aox = re.search(r"aox(-?[\d.]+)", model_dir.parent.name)
    auv = re.search(r"auv(-?[\d.]+)", model_dir.parent.name)
    return {
        "log_T": round(np.log10(p["big_bump_temperature"]), 3),
        "log_U": round(np.log10(p["ionisation_parameter"]), 3),
        "log_nH": round(np.log10(p["hydrogen_density"]), 3),
        "log_NH": float(p["stop_column_density"]),
        "v_turb": float(p["turbulence"]),
        "Z_over_Zsun": round(float(p["metallicities"]) / 0.015, 4),  # GalacticConcordance
        "alpha_ox": float(aox.group(1)) if aox else None,
        "alpha_uv": float(auv.group(1)) if auv else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("run_list")
    parser.add_argument("--timing", nargs="+", default=[],
                        help="TSV(s) written by the SLURM job (one per array task)")
    parser.add_argument("--csv", required=True)
    args = parser.parse_args()

    timing = {}
    for path in args.timing:
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                timing[(row["model_dir"], row["index"])] = row

    rows = []
    for line in Path(args.run_list).read_text().split("\n"):
        if not line.strip():
            continue
        model_dir, index = line.split()
        d = Path(model_dir)
        info = parse_out(d / f"{index}.out")
        t = timing.get((model_dir, index), {})
        ok = (info["ended"] and not info["aborted"] and info["converged"]
              and "column dens reached" in info["stop_reason"])  # not NZONE (nend) or a T stop
        rows.append({"grid": d.parent.name, "incident_index": d.name, "index": index,
                     **model_params(d, index), **info,
                     "wall_s": t.get("wall_s"), "exit_code": t.get("exit_code"), "ok": ok})

    with open(args.csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    n_ok = sum(r["ok"] for r in rows)
    print(f"{len(rows)} runs, {n_ok} ok, {len(rows) - n_ok} flagged -> {args.csv}")


if __name__ == "__main__":
    main()

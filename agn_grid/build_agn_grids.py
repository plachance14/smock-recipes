"""Build the synthesizer AGN grids from the CLOUDY runs (dry run by default).

Each grid is built by syncretize's own create_synthesizer_grid.py, run unchanged, on a
"final view" of the runs that shows it only the grid we intend to ship (author, 17 Oct
2026). For each syncretize output directory in <work>/grids/, in its own build directory
<work>/build/<grid>/:

 1. Choose each model's source: the production run, or the latest successful rerun of it
    (RERUN_SOURCES, in order). "Successful" is parse_cloudy_runs.py's ``ok``, not
    syncretize's check, which passes runs cut short at NZONE.
    A model with no ok run is a hole, except a limit cycle (LIMIT_CYCLES): its
    production run reached the full column but never converged, and it is kept, flagged,
    with its cycle amplitude. Each hole is given its fill stencil
    along T (step 5).
 2. Cut the axes (KEEP): a copy of the incident grid holding only the kept values, and
    parameter YAMLs trimmed to match, so syncretize never sees a cut model.
 3. Link the chosen outputs into final/<grid>/<incident_index>/<index>.* (relative
    symlinks), re-indexed to the cut grid by matching parameter values. A hole gets only
    its .out, so syncretize records it as a failure (zeros, failures = 1).
 4. Run create_synthesizer_grid.py --include-spectra on that view.
 5. Fill each hole along T: a cubic in log T through
    converged neighbours, flagged filled, with a leave-one-out error estimate.
 6. Finish the HDF5: restore WeightVariable, add the smock_flags,
    smock_cycle_half_amplitude and smock_fill_error datasets and provenance attributes,
    write MANIFEST.tsv (each model's source), then check the normalisation against the
    incident grid (known answer).
 7. Rename the HDF5 to its shipped name (SHIPPED_RUN_SET_NAMES): the run sets "blr" and
    "dense" ship as "standard_blr" and "gas_enshrouded" (author, 22 Oct 2026; both are
    BLR-type grids in synthesizer's terms). Only the shipped file changes name; run sets,
    YAMLs, work and build directories keep the internal names (README, "Names").

**WeightVariable (author approved 17 Oct 2026).** Our incident grids
are per erg/s of L_bol, with WeightVariable = bolometric_luminosities, so synthesizer
rescales them to per Lsun whenever it loads them. syncretize loads the incident grid
through synthesizer while building, so an unmodified build comes out per Lsun and is
then rescaled again on every load: 3.83e33 too bright. So the build's incident copy has
WeightVariable "None" (nothing is rescaled during the build; everything stays per erg/s),
and step 6 sets the output's WeightVariable back to bolometric_luminosities, so
synthesizer converts once, on load. Each grid carries this note in its
``smock_build_note`` attribute. If the normalisation ever looks wrong, start there.

The core grids hold ~6 GB of spectra in memory while building: run those through SLURM.

Usage (from the top of smock-recipes, which holds cloudy_runs/; see README.md):
    python agn_grid/build_agn_grids.py                    # dry run, every grid
    python agn_grid/build_agn_grids.py --grids dense_fiducial-turbulence --apply
    python agn_grid/build_agn_grids.py --work-dir /path/to/kept/cloudy_runs \
        --build-dir runs/agn_grid_rebuild --grids dense_core --apply   # from a kept copy of the runs
"""

import argparse
import csv
import datetime
import os
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml
from synthesizer_grids.cloudy.utils import get_cloudy_params, get_grid_props_cloudy

from make_cloudy_inputs import PARAMS, RUN_SETS, SYNCRETIZE_CLOUDY, syncretize_grid_name

GRID_RUN_SETS = ("blr_core", "dense_core", "blr_variations", "dense_variations")

# Kept axis values, per family (the base YAML's prefix). Only the incident (big-bump
# temperature) axis is cut so far (cloudy_grid_report.md §6).
KEEP = {
    "blr": {"big_bump_temperature": lambda T: np.log10(T) >= 4.25 - 1e-6},
    "dense": {"big_bump_temperature": lambda T: np.log10(T) >= 4.2 - 1e-6},
}

# Rerun stages under <work>/reruns/ as (run list, its parse CSV), oldest first: an ok run
# in a later stage supersedes an earlier source. The CSV rows follow the run list's
# order; the run list says where each run lives (a stage's list can include runs in a
# subdirectory, e.g. fix_test_2/nudge/). Deck changes: cloudy_grid_report.md §5.
RERUN_SOURCES = [
    ("fix_test/runs.txt", "fix_test/runs.csv"),              # v_turb 0 without the turbulence line
    ("fix_test_2/runs.txt", "fix_test_2/runs.csv"),          # masers off; 30 iterations; nudge
    ("fix_test_3/runs.txt", "fix_test_3/runs.csv"),          # set pressure ionize 30000
    ("fix_test_4/runs_main.txt", "fix_test_4/runs_main.csv"),  # 100 iterations; cap; nudges
    ("fix_test_4/runs_slow.txt", "fix_test_4/runs_slow.csv"),
    ("fix_test_5/runs.txt", "fix_test_5/runs.csv"),          # failures 100, pressure ionize 300000
    ("fix_test_6/runs.txt", "fix_test_6/runs.csv"),          # failures 100 (BLR Z variation ConvFail)
]
# Reruns never used, keyed (rerun directory, grid, incident index, index).
IGNORE = {
    ("fix_test_4/nudge", "agn-cloudy_T20_aox-1.40_auv-0.50_cloudy-dense_core", 5, 163):
        "cross-check of the fix_test_2 nudge (10^11.499 vs 10^11.501; equal to <= 0.001 dex)",
}

# Runs that never converge because their line optical depths cycle (cycle test, job
# 1388749; cloudy_grid_report.md §4). Their production spectrum is kept, flagged limit_cycle,
# with the cycle's half-amplitude in the bands as its error (author, 18 Oct 2026): more
# iterations don't help, and the per-iteration saves hold only the continuum, so a
# cycle-mean spectrum would not match its lines. Written by
# reruns/cycle_test/analyse_cycles.py.
LIMIT_CYCLES = "cycle_test/cycle_amplitudes.csv"

# Bits of the smock_flags dataset. masers_off marks runs computed with "database H-like
# masers off", which changes healthy runs at the ~1% level (cloudy_grid_report.md §5).
FLAGS = {"rerun": 1, "nudged": 2, "thick_es": 4, "hole": 8, "masers_off": 16, "limit_cycle": 32,
         "filled": 64}

# Holes are filled along the big-bump temperature (author, 18 Oct 2026): synthesizer 1.3.1 would interpolate a hole's zero spectrum. A cubic in
# log T through T-2, T-1, T+1, T+2, or, when T+1 (T-1) is a hole too, through T-2, T-1,
# T+2, T+3 (T-3, T-2, T+1, T+2). Leave-one-out tests next to the holes:
# archive/2026-10-18_cloudy_hole_fill_loo (the method) and
# archive/2026-10-21_hole_fill_components (each stored spectrum, the zero rule and the
# clip of every part to its bracketing neighbours, author 22 Oct 2026).
FILL_AXIS = "big_bump_temperature"
FILL_STENCILS = ((-2, -1, 1, 2), (-2, -1, 2, 3), (-3, -2, 1, 2))

# The fill error is measured in the bands the cycle test used.
NIRCAM_WIDE = ("F090W", "F115W", "F150W", "F200W", "F277W", "F356W", "F444W")
BAND_REDSHIFTS = (5, 6, 7, 8)

# Internal run-set name -> shipped name, as it appears in a grid name after "_cloudy-"
SHIPPED_RUN_SET_NAMES = {"blr": "standard_blr", "dense": "gas_enshrouded"}

OUTPUT_EXTENSIONS = ("cont", "emergent_elin", "intrinsic_elin", "out", "ovr", "in", "yaml")

BUILD_NOTE = (
    "Built by smock-recipes agn_grid/build_agn_grids.py with syncretize's "
    "create_synthesizer_grid.py from an incident-grid copy with WeightVariable 'None', so "
    "nothing was rescaled during the build and the grid is per erg/s of L_bol; "
    "WeightVariable was then set to bolometric_luminosities, so synthesizer rescales to "
    "per Lsun once, on load. An unmodified build is 3.83e33 too bright under an Lsun "
    "default unit (17 Oct 2026)."
)


def shipped_name(name):
    """A grid's shipped name: its run-set name swapped by SHIPPED_RUN_SET_NAMES."""
    for internal, shipped in SHIPPED_RUN_SET_NAMES.items():
        if f"_cloudy-{internal}_" in name:
            return name.replace(f"_cloudy-{internal}_", f"_cloudy-{shipped}_")
    raise ValueError(f"no shipped name for {name}")


def rename_to_shipped(path, internal_name):
    """Step 7: rename the finished HDF5 and record its internal name in it."""
    shipped = path.with_name(shipped_name(internal_name) + ".hdf5")
    with h5py.File(path, "r+") as f:
        f.attrs["smock_internal_name"] = internal_name
    path.rename(shipped)
    return shipped


def grids_to_build():
    """(grid name, incident, base YAML, extra YAML or None) for every production grid."""
    return [(syncretize_grid_name(inc, base, extra), inc, base, extra)
            for run_set in GRID_RUN_SETS for inc, base, extra in RUN_SETS[run_set]]


def family(base):
    """'blr' or 'dense', from the base YAML's name."""
    return base.split("_")[0]


def key(row):
    """A parsed run's (grid, incident index, index)."""
    return row["grid"], int(row["incident_index"]), int(row["index"])


def read_rerun_stage(reruns, run_list, csv_path):
    """[(row, rerun directory, model directory)] for one rerun stage.

    The rerun directory is the run's location under reruns/ above <grid>/<incident>
    (e.g. "fix_test_2/nudge"). Paths are taken relative to reruns/, so the work
    directory can move.
    """
    lines = [line.split() for line in (reruns / run_list).read_text().splitlines() if line.strip()]
    with open(reruns / csv_path) as f:
        rows = list(csv.DictReader(f))
    if len(rows) != len(lines):
        sys.exit(f"{csv_path} has {len(rows)} rows for {len(lines)} runs in {run_list}")
    stage = []
    for row, (model_dir, index) in zip(rows, lines):
        rel = Path(model_dir.split("/reruns/", 1)[1])
        if (rel.parent.name, int(rel.name), int(index)) != key(row):
            sys.exit(f"{csv_path} and {run_list} disagree at {model_dir} {index}")
        stage.append((row, str(rel.parent.parent), reruns / rel))
    return stage


def choose_sources(work):
    """Each model's chosen run: {(grid, incident, index): (row, source, model directory)}.

    source is "production" or a rerun directory. Rerun stages not parsed yet are
    reported and skipped.
    """
    with open(work / "production_runs.csv") as f:
        chosen = {key(r): (r, "production", work / "grids" / r["grid"] / r["incident_index"])
                  for r in csv.DictReader(f)}
    reruns = work / "reruns"
    for run_list, csv_path in RERUN_SOURCES:
        if not (reruns / csv_path).exists():
            print(f"warning: {reruns / csv_path} not found (not parsed yet?); its reruns are not used")
            continue
        for row, source, model_dir in read_rerun_stage(reruns, run_list, csv_path):
            if row["ok"] == "True" and (source, *key(row)) not in IGNORE:
                chosen[key(row)] = (row, source, model_dir)
    return chosen


def read_limit_cycles(work):
    """{(grid, incident, index): cycle half-amplitude}, the largest over every NIRCam wide
    band at z = 5-8 (fractional; the window of iterations 11-30)."""
    with open(work / "reruns" / LIMIT_CYCLES) as f:
        rows = list(csv.DictReader(f))
    cycles = {}
    for r in rows:
        grid_incident, index = r["key"].split()
        grid, incident = grid_incident.split("/")
        cycles[(grid, int(incident), int(index))] = max(float(r["half_amp_nonlya"]),
                                                       float(r["half_amp_lya"]))
    return cycles


def is_full_column_unconverged(row):
    """True for a run that ended normally at its stop column but didn't converge."""
    return (row["exited_ok"] == "True" and row["aborted"] == "False"
            and row["stop_reason"] == "H column dens reached." and row["converged"] == "False")


def incident_axes(path):
    """The incident grid's axis names and values, as stored."""
    with h5py.File(path) as f:
        axes = [a for a in f.attrs["axes"]]
        return axes, {a: f[f"axes/{a}"][:] for a in axes}


def kept_indices(axes, values, keep):
    """Per axis, the indices of the kept values."""
    return [np.flatnonzero([keep.get(a, lambda v: True)(v) for v in values[a]]) for a in axes]


def model_grid(axes, values):
    """syncretize's model order: (value tuples, grid-index tuples), one per model.

    The index tuples are where syncretize puts each model in the grid arrays. They are not
    a C-order unravelling of the model number, so always use them.
    """
    if not axes:
        return [()], [()]
    *_, model_list, index_list = get_grid_props_cloudy(axes, values, verbose=False)
    return [tuple(m) for m in model_list], [tuple(int(k) for k in ix) for ix in index_list]


def grid_index(plan, i, j):
    """Model (incident i, photoionisation j)'s index into the grid arrays."""
    return plan["inc_index"][i] + plan["photo_index"][j]


def shift(idx, axis, offset):
    """idx moved by offset along axis."""
    return idx[:axis] + (idx[axis] + offset,) + idx[axis + 1:]


def photoionisation_axes(param_dir, base, extra):
    """The varied photoionisation axes, as syncretize reads them (sorted names)."""
    _, variable = get_cloudy_params(base, param_dir=param_dir)
    if extra:
        variable |= get_cloudy_params(extra, param_dir=param_dir)[1]
    axes = sorted(variable)
    return axes, {a: variable[a] for a in axes}


def write_cut_incident(src, dst, keep_idx):
    """Copy the incident grid with only the kept axis values, and WeightVariable "None".

    Every dataset whose leading dimensions are the grid shape is sliced; the rest
    (wavelengths) are copied. See the module docstring on WeightVariable.
    """
    with h5py.File(src) as f, h5py.File(dst, "w") as g:
        shape = tuple(len(f[f"axes/{a}"]) for a in f.attrs["axes"])

        def copy(name, obj):
            if isinstance(obj, h5py.Group):
                g.require_group(name).attrs.update(obj.attrs)
                return
            data = obj[()]
            if name.startswith("axes/"):
                axis = list(f.attrs["axes"]).index(name.split("/", 1)[1])
                data = data[keep_idx[axis]]
            elif data.shape[:len(shape)] == shape:
                data = data[np.ix_(*keep_idx)]
            g.create_dataset(name, data=data).attrs.update(obj.attrs)

        f.visititems(copy)
        g.attrs.update(f.attrs)
        g.attrs["WeightVariable"] = "None"
        g.attrs["smock_build_note"] = ("Build-time copy for syncretize: kept axis values "
                                       "only, WeightVariable 'None'.")


def write_cut_yaml(src, dst, axes, kept_values):
    """Copy a parameter YAML with each varied axis trimmed to its kept values."""
    params = yaml.safe_load(src.read_text())
    for axis in axes:
        if axis in params and isinstance(params[axis], list):
            params[axis] = [v for v in params[axis] if float(v) in kept_values[axis]]
        elif "." in axis:
            raise NotImplementedError(f"cutting the nested axis {axis}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(yaml.safe_dump(params, sort_keys=False))


def link(dst, target):
    """A relative symlink dst -> target, so the view survives moving the work directory."""
    dst.symlink_to(os.path.relpath(target.resolve(), dst.parent))


def plan_grid(work, name, incident, base, extra, chosen, cycles):
    """Everything needed to build one grid; nothing is written."""
    keep = KEEP[family(base)]
    inc_axes, inc_values = incident_axes(work / "incident_grids" / f"{incident}.hdf5")
    inc_keep = kept_indices(inc_axes, inc_values, keep)
    photo_axes, photo_values = photoionisation_axes(PARAMS, base, extra)
    photo_keep = kept_indices(photo_axes, photo_values, keep)

    full_inc, _ = model_grid(inc_axes, inc_values)
    full_photo, _ = model_grid(photo_axes, photo_values)
    cut_inc, inc_index = model_grid(inc_axes, {a: inc_values[a][k] for a, k in zip(inc_axes, inc_keep)})
    cut_photo, photo_index = model_grid(photo_axes,
                                        {a: photo_values[a][k] for a, k in zip(photo_axes, photo_keep)})
    old_inc = {v: i for i, v in enumerate(full_inc)}
    old_photo = {v: i for i, v in enumerate(full_photo)}

    models = []  # (new incident, new index, old incident, old index, model dir, source, flags)
    cycle_amp = {}  # (new incident, new index): half-amplitude, for the limit cycles kept
    for i, inc_v in enumerate(cut_inc):
        for j, photo_v in enumerate(cut_photo):
            oi, oj = old_inc[inc_v], old_photo[photo_v]
            row, source, model_dir = chosen[(name, oi, oj)]
            limit_cycle = row["ok"] != "True" and (name, oi, oj) in cycles
            if row["ok"] != "True" and not limit_cycle:  # no usable run: keep production's .out only
                models.append((i, j, oi, oj, model_dir, "hole", FLAGS["hole"]))
                continue
            flags = 0 if source == "production" else FLAGS["rerun"]
            if limit_cycle:
                if source != "production" or not is_full_column_unconverged(row):
                    sys.exit(f"limit cycle {name} {oi}/{oj}: expected an unconverged full-column "
                             f"production run, got {source} ({row['stop_reason']!r})")
                flags |= FLAGS["limit_cycle"]
                cycle_amp[(i, j)] = cycles[(name, oi, oj)]
            if "nudge" in source:
                flags |= FLAGS["nudged"]
            if row["thick_es"] == "True":
                flags |= FLAGS["thick_es"]
            if "masers off" in (model_dir / f"{oj}.in").read_text():
                flags |= FLAGS["masers_off"]
            models.append((i, j, oi, oj, model_dir, source, flags))
    plan = {
        "name": name, "incident": incident, "base": base, "extra": extra,
        "inc_axes": inc_axes, "inc_keep": inc_keep, "photo_axes": photo_axes,
        "photo_values": photo_values, "photo_keep": photo_keep,
        "shape": tuple(len(k) for k in inc_keep) + tuple(len(k) for k in photo_keep),
        "full_shape": tuple(len(inc_values[a]) for a in inc_axes)
                      + tuple(len(photo_values[a]) for a in photo_axes),
        "inc_index": inc_index, "photo_index": photo_index,
        "models": models, "cycle_amp": cycle_amp,
    }
    plan["fills"] = plan_fills(plan)
    return plan


def flag_array(plan):
    """The models' flags as a grid-shaped array."""
    flags = np.zeros(plan["shape"], dtype=np.int32)
    for i, j, *_, f in plan["models"]:
        flags[grid_index(plan, i, j)] = f
    return flags


def fill_stencil(flags, idx, axis):
    """The first of FILL_STENCILS usable for the hole at idx, or None.

    Every stencil point must be a converged model (no hole, no limit cycle), and a gap
    stencil is used only when the point it skips is itself a hole.
    """
    bad = FLAGS["hole"] | FLAGS["limit_cycle"]
    n = flags.shape[axis]
    for offs in FILL_STENCILS:
        skipped = set(range(min(offs), max(offs) + 1)) - set(offs) - {0}
        if (all(0 <= idx[axis] + o < n and not flags[shift(idx, axis, o)] & bad for o in offs)
                and all(flags[shift(idx, axis, o)] & FLAGS["hole"] for o in skipped)):
            return offs
    return None


def plan_fills(plan):
    """{(i, j): stencil offsets} for every hole; stops if a hole can't be filled."""
    if FILL_AXIS not in plan["inc_axes"]:
        return {}
    axis = plan["inc_axes"].index(FILL_AXIS)
    flags = flag_array(plan)
    fills = {}
    for i, j, oi, oj, *_, f in plan["models"]:
        if f & FLAGS["hole"]:
            fills[(i, j)] = fill_stencil(flags, grid_index(plan, i, j), axis)
            if fills[(i, j)] is None:
                sys.exit(f"hole {plan['name']} {oi}/{oj}: no usable stencil along {FILL_AXIS}")
    return fills


def describe(plan):
    """One paragraph per grid for the dry run."""
    m = plan["models"]
    sources = {}
    for *_, source, _ in m:
        sources[source] = sources.get(source, 0) + 1
    flags = {k: sum(1 for *_, f in m if f & bit) for k, bit in FLAGS.items()}
    print(f"{plan['name']}\n    shape {plan['full_shape']} -> {plan['shape']} ({len(m)} models)"
          f"\n    sources {sources}\n    flags {flags}")
    orig = {(i, j): (oi, oj) for i, j, oi, oj, *_ in m}
    for (i, j), offs in plan["fills"].items():
        print(f"    fill {orig[(i, j)][0]}/{orig[(i, j)][1]} along {FILL_AXIS}: T{offs}")


def write_view(work, build, plan):
    """The cut incident grid, the cut YAMLs and the final-view symlink tree."""
    name = plan["name"]
    (build / "grids").mkdir(parents=True)
    write_cut_incident(work / "incident_grids" / f"{plan['incident']}.hdf5",
                       build / "grids" / f"{plan['incident']}.hdf5", plan["inc_keep"])
    kept = {a: set(map(float, plan["photo_values"][a][k]))
            for a, k in zip(plan["photo_axes"], plan["photo_keep"])}
    for yaml_name in filter(None, (plan["base"], plan["extra"])):
        write_cut_yaml(PARAMS / yaml_name, build / "params" / yaml_name, plan["photo_axes"], kept)

    for i, j, oi, oj, src, source, flags in plan["models"]:
        dst = build / "final" / name / str(i)
        if not dst.exists():
            dst.mkdir(parents=True)
            # syncretize reads the line list from each incident directory
            link(dst / "linelist.dat", work / "grids" / name / str(oi) / "linelist.dat")
        for ext in ("out",) if flags & FLAGS["hole"] else OUTPUT_EXTENSIONS:
            target = src / f"{oj}.{ext}"
            if target.exists():
                link(dst / f"{j}.{ext}", target)


def run_syncretize(build, plan):
    """syncretize's create_synthesizer_grid.py, unchanged, on the final view."""
    cmd = [sys.executable, str(SYNCRETIZE_CLOUDY / "create_synthesizer_grid.py"),
           "--grid-dir", str(build / "grids"), "--cloudy-output-dir", str(build / "final"),
           "--incident-grid", plan["incident"],
           "--cloudy-paramfile", plan["base"].removesuffix(".yaml"), "--include-spectra"]
    if plan["extra"]:
        cmd += ["--cloudy-paramfile-extra", plan["extra"].removesuffix(".yaml")]
    with open(build / "syncretize.log", "w") as log:
        subprocess.run(cmd, cwd=build, stdout=log, stderr=subprocess.STDOUT, check=True)
    return build / "grids" / f"{plan['name']}.hdf5"


def lagrange_weights(x, xs):
    """Weights of the polynomial through the points xs, evaluated at x."""
    return np.array([np.prod([(x - xj) / (xi - xj) for xj in xs if xj != xi]) for xi in xs])


def interpolate_log(stack, weights):
    """The polynomial in log10 through the rows of stack (one per stencil point).

    Zero wherever any stencil point is zero: a log floor given a negative weight can blow
    the fill up by many dex (archive/2026-10-21_hole_fill_components).
    """
    positive = (stack > 0).all(axis=0)
    logs = np.log10(np.where(positive, stack, 1.0))
    return np.where(positive, 10 ** np.tensordot(weights, logs, axes=1), 0.0)


def fill_model(f, idx, axis, offs):
    """Filled transmitted, nebular and linecont spectra, line luminosities and
    normalisation for the model at idx, from its stencil along axis.

    Everything is filled per ionising photon: CLOUDY's incident Q is the same along T at
    fixed U and n_H (to 2e-4 dex), so a stored quantity divided by the grid's Q
    (log10_specific_ionising_luminosity/HI) is the raw CLOUDY output up to one constant,
    which is what the leave-one-out tests interpolated. The nebular is filled as nebular
    continuum + linecont, so the continuum can't go negative. Every part is clipped to its
    two bracketing neighbours (author, 22 Oct 2026): where a part is heavily absorbed or
    rises steeply in T, the cubic overshoots (25x in one test; linecont 5-9x at 1658 A in
    the dense core's adjacent pair).
    """
    points = [shift(idx, axis, o) for o in offs]
    log_t = np.log10(f[f"axes/{FILL_AXIS}"][:])
    weights = lagrange_weights(log_t[idx[axis]], log_t[[p[axis] for p in points]])
    q = {p: 10 ** f["log10_specific_ionising_luminosity/HI"][p] for p in points + [idx]}
    below = offs.index(max(o for o in offs if o < 0))
    above = offs.index(min(o for o in offs if o > 0))

    def stack(name):
        return np.array([f[name][p] / q[p] for p in points])

    def fill_clipped(stack):
        bracket = stack[[below, above]]
        return np.clip(interpolate_log(stack, weights), bracket.min(0), bracket.max(0))

    transmitted, nebular, linecont = (stack(f"spectra/{n}") for n in ("transmitted", "nebular", "linecont"))
    out = {
        "spectra/transmitted": fill_clipped(transmitted),
        "spectra/linecont": fill_clipped(linecont),
        "lines/luminosity": fill_clipped(stack("lines/luminosity")),
    }
    out["spectra/nebular"] = fill_clipped(nebular - linecont) + out["spectra/linecont"]
    out = {k: v * q[idx] for k, v in out.items()}
    # normalisation = grid Q / Q of CLOUDY's incident, and the latter is fixed along T
    log_q_cloudy = np.mean([np.log10(q[p] / f["spectra/normalisation"][p]) for p in points])
    out["spectra/normalisation"] = q[idx] / 10 ** log_q_cloudy
    return out


def line_continua(spectra, lam, line_lam):
    """syncretize's line continuum datasets: the spectra interpolated at the lines."""
    transmitted = spectra["spectra/transmitted"]
    nebular_continuum = spectra["spectra/nebular"] - spectra["spectra/linecont"]
    continua = {"incident": spectra["spectra/incident"], "transmitted": transmitted,
                "nebular_continuum": nebular_continuum,
                "total_continuum": nebular_continuum + transmitted}
    return {f"lines/{k}": np.interp(line_lam, lam, v) for k, v in continua.items()}


def band_means(lam, lnu, filters):
    """Photon-weighted mean L_nu in each NIRCam wide band at z = 5-8 (lam rest-frame A)."""
    nu = 2.998e18 / lam
    means = []
    for z in BAND_REDSHIFTS:
        for flt in filters:
            w = np.interp(lam * (1 + z), flt.lam.value, flt.t, left=0, right=0) / nu
            means.append(np.trapezoid(lnu * w, nu) / np.trapezoid(w, nu))
    return np.array(means)


def fill_error(f, flags, idx, axis, offs, n_inc, filters):
    """Leave-one-out error of the fill next to the hole at idx.

    Each converged neighbour along the photoionisation axes whose own stencil (the same
    offsets) is converged is refilled and compared with its CLOUDY spectrum. The error is
    the largest |filled / true - 1| of transmitted + nebular (what synthesizer adds) over
    the bands where that is >= 1% of the incident (where the fill is reliable, as in the
    tests), the largest over neighbours. NaN if no neighbour qualifies.
    """
    bad = FLAGS["hole"] | FLAGS["limit_cycle"] | FLAGS["filled"]
    lam = f["spectra/wavelength"][:]
    worst = np.nan
    for ax in range(n_inc, flags.ndim):
        for d in (-1, 1):
            c = shift(idx, ax, d)
            if not 0 <= c[ax] < flags.shape[ax] or flags[c] & bad:
                continue
            if any(not 0 <= c[axis] + o < flags.shape[axis] or flags[shift(c, axis, o)] & bad
                   for o in offs):
                continue
            pred = fill_model(f, c, axis, offs)
            true = f["spectra/transmitted"][c] + f["spectra/nebular"][c]
            true_b = band_means(lam, true, filters)
            pred_b = band_means(lam, pred["spectra/transmitted"] + pred["spectra/nebular"], filters)
            use = true_b >= 0.01 * band_means(lam, f["spectra/incident"][c], filters)
            if use.any():
                worst = np.nanmax([worst, np.max(np.abs(pred_b[use] / true_b[use] - 1))])
    return worst


def fill_holes(path, plan):
    """Fill every hole of the built grid along T (plan["fills"]), in place.

    Writes the filled spectra, lines, line continua and normalisation, sets failures to 0,
    and marks each model filled instead of hole in plan["models"]. The incident is copied
    from a model at the same incident index (it doesn't depend on the gas). Returns
    {(i, j): leave-one-out error}.
    """
    if not plan["fills"]:
        return {}
    from synthesizer.instruments import FilterCollection
    filters = FilterCollection(filter_codes=[f"JWST/NIRCam.{b}" for b in NIRCAM_WIDE])
    axis = plan["inc_axes"].index(FILL_AXIS)
    flags = flag_array(plan)
    errors = {}
    with h5py.File(path, "r+") as f:
        if set(f.attrs["spec_names"]) != {"incident", "transmitted", "nebular", "linecont"}:
            sys.exit(f"{path}: unexpected spectra {list(f.attrs['spec_names'])}")
        lam, line_lam = f["spectra/wavelength"][:], f["lines/wavelength"][:]
        for (i, j), offs in plan["fills"].items():
            idx = grid_index(plan, i, j)
            filled = fill_model(f, idx, axis, offs)
            donor = next(grid_index(plan, i, k) for k in range(len(plan["photo_index"]))
                         if not flags[grid_index(plan, i, k)] & FLAGS["hole"])
            filled["spectra/incident"] = f["spectra/incident"][donor]
            filled |= line_continua(filled, lam, line_lam)
            for name, value in filled.items():
                f[name][idx] = value
            f["failures"][idx] = 0
            errors[(i, j)] = fill_error(f, flags, idx, axis, offs, len(plan["inc_axes"]), filters)
    plan["models"] = [(i, j, oi, oj, d, "filled", (fl & ~FLAGS["hole"]) | FLAGS["filled"])
                      if (i, j) in errors else (i, j, oi, oj, d, source, fl)
                      for i, j, oi, oj, d, source, fl in plan["models"]]
    return errors


def finish_grid(path, plan, build, fill_errors):
    """Restore WeightVariable, add smock_flags, smock_cycle_half_amplitude,
    smock_fill_error and provenance attributes, and write MANIFEST.tsv."""
    flags = flag_array(plan)
    cycle_amp = np.where(flags & FLAGS["hole"], np.nan, 0.0)
    fill_error = np.where(flags & FLAGS["hole"], np.nan, 0.0)
    for (i, j), amp in plan["cycle_amp"].items():
        cycle_amp[grid_index(plan, i, j)] = amp
    for (i, j), err in fill_errors.items():
        fill_error[grid_index(plan, i, j)] = err
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                            text=True, cwd=Path(__file__).parent).stdout.strip()
    with h5py.File(path, "r+") as f:
        f.attrs["WeightVariable"] = "bolometric_luminosities"
        f.attrs["smock_build_note"] = BUILD_NOTE
        f.attrs["smock_build"] = f"{datetime.date.today()} smock-recipes {commit}; manifest {build / 'MANIFEST.tsv'}"
        d = f.create_dataset("smock_flags", data=flags)
        d.attrs["Description"] = ("Per-model bit flags: " + ", ".join(f"{v} {k}" for k, v in FLAGS.items())
                                  + ". A hole has zero spectra and failures = 1.")
        d = f.create_dataset("smock_cycle_half_amplitude", data=cycle_amp)
        d.attrs["Description"] = (
            "Per-model fractional error of a limit_cycle model's spectrum: (max - min) / "
            "(2 mean) of the band luminosity over CLOUDY iterations 11-30, the largest over "
            "the 7 NIRCam wide bands at z = 5-8. 0 for converged models, NaN for holes.")
        d = f.create_dataset("smock_fill_error", data=fill_error)
        d.attrs["Description"] = (
            "Per-model fractional error of a filled model's spectrum (a cubic in log T): the "
            "largest leave-one-out |filled / true - 1| of transmitted + nebular among its "
            "converged neighbours along the photoionisation axes, over the 7 NIRCam wide "
            "bands at z = 5-8 where that is >= 1% of the incident (heavily absorbed bands are "
            "not reliable). 0 for unfilled models, NaN for holes.")

    with open(build / "MANIFEST.tsv", "w") as f:
        f.write("incident_index\tindex\torig_incident_index\torig_index\tsource\tflags\n")
        for i, j, oi, oj, _, source, fl in plan["models"]:
            f.write(f"{i}\t{j}\t{oi}\t{oj}\t{source}\t{fl}\n")


def check_normalisation(path, plan, incident_dir):
    """Known answer: loaded, the grid's incident spectrum equals the incident grid's own.

    Both go through synthesizer's Grid, so both get its per-Lsun conversion once.
    Returns the largest |ratio - 1| over 1000-10000 A among non-hole models.
    """
    from synthesizer.grid import Grid
    ref = Grid(plan["incident"], grid_dir=str(incident_dir), ignore_lines=True)
    out = Grid(path.stem, grid_dir=str(path.parent), ignore_lines=True)
    lam_r, lam_o = ref.lam.value, out.lam.value
    window = (lam_r > 1000) & (lam_r < 10000)
    ref_spec = ref.spectra["incident"][np.ix_(*plan["inc_keep"])]
    out_spec = out.spectra["incident"]
    n_inc = len(plan["inc_keep"])
    worst = 0.0
    for i, j, *_, f in plan["models"]:
        if f & FLAGS["hole"]:
            continue
        idx = grid_index(plan, i, j)
        ratio = np.interp(lam_r[window], lam_o, out_spec[idx]) / ref_spec[idx[:n_inc]][window]
        worst = max(worst, abs(float(np.median(ratio)) - 1))
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--work-dir", default="cloudy_runs")
    # Builds go to <work>/build unless given elsewhere, e.g. to rebuild from a kept copy of
    # the runs without writing there.
    parser.add_argument("--build-dir", help="where the build directories go (default: <work-dir>/build)")
    parser.add_argument("--grids", nargs="+", help="substrings of grid names (default: all)")
    parser.add_argument("--rebuild", action="store_true", help="replace existing build directories")
    parser.add_argument("--apply", action="store_true", help="build (default: dry run)")
    args = parser.parse_args()

    work = Path(args.work_dir).resolve()
    build_root = Path(args.build_dir).resolve() if args.build_dir else work / "build"
    chosen = choose_sources(work)
    cycles = read_limit_cycles(work)
    grids = [g for g in grids_to_build() if not args.grids or any(s in g[0] for s in args.grids)]
    plans = [plan_grid(work, *g, chosen, cycles) for g in grids]
    # every limit cycle should land in a grid; one that doesn't means a key mismatch
    kept = {(p["name"], oi, oj) for p in plans
            for _, _, oi, oj, _, _, flags in p["models"] if flags & FLAGS["limit_cycle"]}
    if not args.grids and kept != set(cycles):
        sys.exit(f"{len(set(cycles) - kept)} limit cycles in {LIMIT_CYCLES} are in no grid")
    for plan in plans:
        describe(plan)
        print(f"    ships as {shipped_name(plan['name'])}.hdf5")
    if not args.apply:
        print("Dry run: nothing written. Re-run with --apply.")
        return

    for plan in plans:
        build = build_root / plan["name"]
        if build.exists():
            if not args.rebuild:
                sys.exit(f"{build} exists; use --rebuild to replace it.")
            shutil.rmtree(build)
        write_view(work, build, plan)
        path = run_syncretize(build, plan)
        fill_errors = fill_holes(path, plan)
        finish_grid(path, plan, build, fill_errors)
        for (i, j), err in fill_errors.items():
            print(f"    filled {i}/{j}: leave-one-out error {err:.1%}")
        worst = check_normalisation(path, plan, work / "incident_grids")
        status = "ok" if worst < 0.01 else "FAILED"
        print(f"{path}\n    normalisation check {status}: max |median ratio - 1| = {worst:.2e}")
        if status != "ok":
            sys.exit("normalisation check failed")
        print(f"    shipped as {rename_to_shipped(path, plan['name']).name}")


if __name__ == "__main__":
    main()

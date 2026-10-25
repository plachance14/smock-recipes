"""Stage edited copies of CLOUDY inputs for a rerun, leaving the originals untouched.

Reads a run list ("<model_dir> <index>" per line, as make_cloudy_inputs.py writes) and
copies each run's <index>.in and <index>.yaml to <out_dir>/<grid>/<incident_index>/.
The shared files (input.sed, input.npy, linelist.dat) are symlinked from the original
model directory, by relative links so the work directory can move. Keeping the
<grid>/<incident_index> layout lets parse_cloudy_runs.py read the reruns like the
originals. Edits:
    --add-command CMD         append a CLOUDY command, e.g. "set pressure ionize 30000"
                              (the decks have no blank line, so CLOUDY reads to the end).
                              CLOUDY only matches keywords at the start of a word
                              (parser.h nMatch1), so a run-together "set presioniz"
                              fails: write each keyword as its own word. "{index}" is
                              replaced by the run's index, for per-run save files, e.g.
                              'save convergence reason "{index}.convreason"'
    --omit-zero-turbulence    drop "turbulence 0.0", which CLOUDY refuses
                              (make_cloudy_inputs.omit_zero_turbulence)
    --max-iterations N        "iterate to convergence" -> "iterate to convergence N"
                              (CLOUDY's default limit is 10; the first number on the
                              line is the limit, parse_commands.cpp ParseIterations)
    --replace-line OLD NEW    replace the deck line OLD (exact, must occur once) with NEW,
                              e.g. "hden 11.5 log" "hden 11.501 log" to nudge a crashed run
With no edit, the input is copied unchanged (a plain rerun).

The staged runs are appended to <out_dir>/runs.txt, for run_cloudy.slurm (RUN_LIST=).
Existing staged inputs are never overwritten.

Usage (from the top of smock-recipes, which holds cloudy_runs/; see README.md):
    python agn_grid/make_rerun_inputs.py RUNS.txt cloudy_runs/reruns/fix_test \
        --add-command "set pressure ionize 30000"      # dry run
    python agn_grid/make_rerun_inputs.py ... --apply
"""

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

from make_cloudy_inputs import omit_zero_turbulence

SHARED = ("input.sed", "input.npy", "linelist.dat")


def edit_deck(text, add_commands, omit_turbulence, max_iterations=None, replace_lines=(), index=None):
    """The deck with the requested edits applied ("{index}" in added commands -> index)."""
    for old, new in replace_lines:
        lines = text.split("\n")
        if lines.count(old) != 1:
            sys.exit(f"expected the line {old!r} exactly once")
        text = "\n".join(new if line == old else line for line in lines)
    if omit_turbulence:
        text = omit_zero_turbulence(text)
    if max_iterations:
        text, n = re.subn(r"^iterate to convergence\s*$", f"iterate to convergence {max_iterations}",
                          text, flags=re.MULTILINE)
        if n != 1:
            sys.exit("expected exactly one bare 'iterate to convergence' line")
    for cmd in add_commands:
        text += cmd.replace("{index}", str(index)) + "\n"
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("run_list")
    parser.add_argument("out_dir")
    parser.add_argument("--add-command", action="append", default=[])
    parser.add_argument("--omit-zero-turbulence", action="store_true")
    parser.add_argument("--max-iterations", type=int)
    parser.add_argument("--replace-line", nargs=2, action="append", default=[], metavar=("OLD", "NEW"))
    parser.add_argument("--apply", action="store_true", help="write (default: dry run)")
    args = parser.parse_args()

    out = Path(args.out_dir).resolve()
    runs = [line.split() for line in Path(args.run_list).read_text().splitlines() if line.strip()]
    staged = []
    for model_dir, index in runs:
        src = Path(model_dir).resolve()
        dst = out / src.parent.name / src.name
        old = (src / f"{index}.in").read_text()
        new = edit_deck(old, args.add_command, args.omit_zero_turbulence, args.max_iterations,
                        args.replace_line, index)
        if (dst / f"{index}.in").exists():
            sys.exit(f"{dst / f'{index}.in'} exists; inputs are not overwritten.")
        if args.omit_zero_turbulence and new == old:
            print(f"warning: no 'turbulence 0' line in {src / f'{index}.in'}")
        staged.append((src, dst, index, new))
    print(f"{len(staged)} runs -> {out}"
          + "".join(f"\n    + {c}" for c in args.add_command)
          + ("\n    - turbulence 0" if args.omit_zero_turbulence else "")
          + (f"\n    iterate to convergence {args.max_iterations}" if args.max_iterations else "")
          + "".join(f"\n    {o!r} -> {n!r}" for o, n in args.replace_line))
    if not args.apply:
        print("Dry run: nothing written. Re-run with --apply.")
        return

    for src, dst, index, new in staged:
        dst.mkdir(parents=True, exist_ok=True)
        for name in SHARED:
            if (src / name).exists() and not (dst / name).exists():
                (dst / name).symlink_to(os.path.relpath((src / name).resolve(), dst.resolve()))
        (dst / f"{index}.in").write_text(new)
        shutil.copyfile(src / f"{index}.yaml", dst / f"{index}.yaml")
    with open(out / "runs.txt", "a") as f:
        f.writelines(f"{dst} {index}\n" for _, dst, index, _ in staged)
    print(f"appended to {out / 'runs.txt'}")


if __name__ == "__main__":
    main()

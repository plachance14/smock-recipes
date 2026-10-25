"""Cut the test groups out of an ASTRID PIG into a small BigFile with the same layout.

    python test_cutout/make_test_cutout.py              # dry run: show what would be written
    python test_cutout/make_test_cutout.py --apply      # write the cutout, then validate it
    python test_cutout/make_test_cutout.py --validate   # re-check an existing cutout

The groups are chosen by subhalo index: the test set picked on 15 Oct 2026 (ROADMAP,
Phase 1). Each subhalo's FoF group is the one whose particle ranges contain the
subhalo's (SubhaloGroupNr is not always the FoF row). The cutout holds:
  - FOFGroups/* and SubGroups/*: every block, for those groups and all their subhalos;
  - 0/, 4/, 5/: the gas, star and BH fields in FIELDS, for every particle of those groups;
  - no dark matter. There is no 1/ block, but the catalogue's type-1 lengths and masses
    are kept as they were, and the type-1 offsets are renumbered like the others.

Offsets are renumbered so the cutout reads like a full PIG. Group starts are the
cumulative sum of LengthByType (as in the full PIG), and subhalo offsets move with their
group. GroupFirstSub and SubhaloGroupNr point at the new rows; the original rows are kept
in FOFGroups/OriginalIndex and SubGroups/OriginalIndex.

Validation reads every subhalo and every group from the full PIG and from the cutout,
block by block, and asserts the arrays are identical. It then checks units: summed
particle masses against SubhaloMassType, and the stellar half-mass radius from the
particles against SubhaloHalfmassRadType.
"""
import argparse
import datetime
from pathlib import Path

import bigfile
import numpy as np

from smock import config

SNAP = "047"                                           # z = 6
TEST_SUBHALOS = (303604, 54437, 24900, 4283857, 42595)  # typical, LRD-like, multi, small, massive
FIELDS = {   # what smock.data reads, plus velocities (and BH Mass for the mass check)
    0: ["Position", "Velocity", "Mass", "Metallicity", "SmoothingLength", "StarFormationRate"],
    4: ["Position", "Velocity", "Mass", "TotalMassReturned", "Metallicity", "StarFormationTime"],
    5: ["ID", "Position", "Velocity", "Mass", "BlackholeMass", "BlackholeAccretionRate"],
}
# Catalogue blocks whose values change in the cutout; every other block is copied as is.
RENUMBERED = {"FOFGroups/OffsetByType", "FOFGroups/GroupFirstSub",
              "SubGroups/SubhaloOffsetType", "SubGroups/SubhaloGroupNr"}


def rows(f, block, i0, i1):
    return np.asarray(f[block][i0:i1])


# ---- finding the groups and their subhalos ---------------------------------------

def group_ranges(f, n_star_needed=0, n_groups_needed=0):
    """Particle starts (exclusive cumsum of LengthByType) and lengths per FoF group, for at
    least n_groups_needed groups and enough to cover the first n_star_needed star particles.
    Shapes (n_groups, 6)."""
    chunks, n, step = [], 0, 1_000_000
    while (not chunks or n < n_groups_needed
           or np.concatenate(chunks)[:, 4].sum() < n_star_needed):
        chunks.append(rows(f, "FOFGroups/LengthByType", n, n + step).astype(np.int64))
        n += step
    length = np.concatenate(chunks)
    start = np.zeros_like(length)
    start[1:] = np.cumsum(length, axis=0)[:-1]
    return start, length


def contained(sub_off, sub_len, g_start, g_len):
    """True for each subhalo whose non-empty particle ranges all lie inside the group's."""
    ok = (sub_len == 0) | ((sub_off >= g_start) & (sub_off + sub_len <= g_start + g_len))
    return ok.all(axis=-1) & (sub_len.sum(axis=-1) > 0)


def group_of(f, sidx, g_start, g_len):
    """Index of the FoF group containing subhalo sidx (searched by its star offset)."""
    off = rows(f, "SubGroups/SubhaloOffsetType", sidx, sidx + 1)[0].astype(np.int64)
    ln = rows(f, "SubGroups/SubhaloLenType", sidx, sidx + 1)[0].astype(np.int64)
    g = int(np.searchsorted(g_start[:, 4], off[4], side="right") - 1)
    assert contained(off, ln, g_start[g], g_len[g]), f"subhalo {sidx} not inside group {g}"
    return g


def subhalos_of(f, sidx, g_start, g_len, window=2000):
    """(first, stop) rows of the contiguous run of subhalos inside the group holding sidx."""
    nsub = f["SubGroups/SubhaloLenType"].size
    while True:
        lo, hi = max(0, sidx - window), min(nsub, sidx + window)
        off = rows(f, "SubGroups/SubhaloOffsetType", lo, hi).astype(np.int64)
        ln = rows(f, "SubGroups/SubhaloLenType", lo, hi).astype(np.int64)
        inside = contained(off, ln, g_start, g_len)
        k = sidx - lo
        first = k - np.argmin(inside[k::-1]) + 1 if not inside[:k + 1].all() else 0
        stop = k + np.argmin(inside[k:]) if not inside[k:].all() else hi - lo
        if (first > 0 or lo == 0) and (stop < hi - lo or hi == nsub):
            return lo + first, lo + stop
        window *= 4


def plan(f, subhalos):
    """Groups (ascending) and their subhalo row ranges, for the chosen subhalos."""
    star_off = [int(rows(f, "SubGroups/SubhaloOffsetType", s, s + 1)[0, 4]) for s in subhalos]
    g_start, g_len = group_ranges(f, max(star_off) + 1)
    groups = sorted({group_of(f, s, g_start, g_len) for s in subhalos})
    sub_runs = []
    for g in groups:
        member = next(s for s in subhalos if group_of(f, s, g_start, g_len) == g)
        sub_runs.append(subhalos_of(f, member, g_start[g], g_len[g]))
    return np.array(groups), sub_runs, g_start, g_len


def describe(f, groups, sub_runs, g_len):
    """Print what the cutout would hold, with its size."""
    nbytes = {t: sum(f[f"{t}/{k}"].dtype.itemsize for k in ks) for t, ks in FIELDS.items()}
    nsubs = rows_of(f, "FOFGroups/GroupNsubs", groups)
    total = 0
    print(f"{'group':>7s} {'subhalos':>17s} {'GroupNsubs':>10s} {'gas':>8s} {'stars':>7s} {'BHs':>4s} {'MB':>6s}")
    for g, (a, b), n in zip(groups, sub_runs, nsubs):
        mb = sum(g_len[g, t] * nbytes[t] for t in FIELDS) / 1e6
        total += mb
        print(f"{g:7d} {a:8d}-{b - 1:8d} {n:10d} {g_len[g, 0]:8d} {g_len[g, 4]:7d} {g_len[g, 5]:4d} {mb:6.1f}")
    print(f"particle data ~{total:.1f} MB (gas {nbytes[0]} B, star {nbytes[4]} B, BH {nbytes[5]} B per particle)")


def rows_of(f, block, idx):
    """Rows idx of a catalogue block, in the given order, reading each run of
    consecutive indices in one go."""
    idx = np.asarray(idx)
    breaks = np.nonzero(np.diff(idx) != 1)[0] + 1
    return np.concatenate([rows(f, block, run[0], run[-1] + 1) for run in np.split(idx, breaks)])


# ---- writing ---------------------------------------------------------------------

def catalogue_blocks(f):
    return sorted(b for b in f.list_blocks() if b.startswith(("FOFGroups/", "SubGroups/")))


def write_cutout(src, out_dir, subhalos):
    groups, sub_runs, g_start, g_len = plan(src, subhalos)
    sub_rows = np.concatenate([np.arange(a, b) for a, b in sub_runs])
    n_sub = np.array([b - a for a, b in sub_runs])
    new_g_len = g_len[groups]
    new_g_start = np.zeros_like(new_g_len)
    new_g_start[1:] = np.cumsum(new_g_len, axis=0)[:-1]
    sub_group = np.repeat(np.arange(groups.size), n_sub)   # new group row of each subhalo
    shift = (new_g_start - g_start[groups])[sub_group]      # per subhalo and particle type

    with bigfile.File(str(out_dir), create=True) as dst:
        for block in catalogue_blocks(src):
            if block.startswith("FOFGroups/"):
                data = rows_of(src, block, groups)
            else:
                data = np.concatenate([rows(src, block, a, b) for a, b in sub_runs])
            if block == "FOFGroups/OffsetByType":
                data = new_g_start.astype(data.dtype)
            elif block == "FOFGroups/GroupFirstSub":
                data = np.r_[0, np.cumsum(n_sub)[:-1]].astype(data.dtype)
            elif block == "SubGroups/SubhaloOffsetType":
                data = (data.astype(np.int64) + shift).astype(data.dtype)
            elif block == "SubGroups/SubhaloGroupNr":
                data = sub_group.astype(data.dtype)
            dst.create_from_array(block, data, Nfile=1)
        dst.create_from_array("FOFGroups/OriginalIndex", groups.astype(np.int64), Nfile=1)
        dst.create_from_array("SubGroups/OriginalIndex", sub_rows.astype(np.int64), Nfile=1)

        for t, fields in FIELDS.items():
            for field in fields:
                block = f"{t}/{field}"
                data = np.concatenate([rows(src, block, g_start[g, t], g_start[g, t] + g_len[g, t])
                                       for g in groups])
                dst.create_from_array(block, data, Nfile=1)

        # String attributes don't read back with bigfile 0.1.52 + numpy 2, so text goes in README.txt.
        with dst.create("Header") as h:
            for key in src["Header"].attrs.keys():
                h.attrs[key] = src["Header"].attrs[key]
            n_part = new_g_len.sum(axis=0)
            n_part[[1, 2, 3]] = 0                           # only types 0, 4, 5 are stored
            h.attrs["NumPartInGroupTotal"] = n_part.astype(np.int64)
            h.attrs["NumFOFGroupsTotal"] = np.array([groups.size], dtype=np.int64)
            h.attrs["CutoutSubhalos"] = np.array(subhalos, dtype=np.int64)
            h.attrs["CutoutGroups"] = groups.astype(np.int64)
    return groups, sub_runs


def write_readme(out_dir, src_path, subhalos, groups):
    (out_dir / "README.txt").write_text(
        f"ASTRID test cutout, written {datetime.date.today()} by smock-recipes/test_cutout/make_test_cutout.py\n"
        f"Source: {src_path}\n"
        f"Chosen subhalos (original indices): {list(subhalos)}\n"
        f"Their FoF groups (original indices): {groups.tolist()}\n"
        "Layout: as a PIG (Header, FOFGroups/, SubGroups/, 0/, 4/, 5/); no dark matter (1/).\n"
        "Offsets and GroupFirstSub/SubhaloGroupNr are renumbered for the cutout; original rows are\n"
        "in FOFGroups/OriginalIndex and SubGroups/OriginalIndex. See the script's docstring.\n"
        "ASTRID: https://astrid.psc.edu/\n")


# ---- validation ------------------------------------------------------------------

def check_identical(src, cut):
    """Every group and subhalo, read from both files by their own offsets, must match."""
    g_orig = np.asarray(cut["FOFGroups/OriginalIndex"][:])
    s_orig = np.asarray(cut["SubGroups/OriginalIndex"][:])
    g_start_src, g_len_src = group_ranges(src, n_groups_needed=g_orig.max() + 1)
    g_len_cut = np.asarray(cut["FOFGroups/LengthByType"][:]).astype(np.int64)
    g_off_cut = np.asarray(cut["FOFGroups/OffsetByType"][:]).astype(np.int64)
    assert np.array_equal(g_off_cut[1:], np.cumsum(g_len_cut, axis=0)[:-1]), "cutout offsets != cumsum"
    s_off_cut = np.asarray(cut["SubGroups/SubhaloOffsetType"][:]).astype(np.int64)
    s_len_cut = np.asarray(cut["SubGroups/SubhaloLenType"][:]).astype(np.int64)
    s_off_src = rows_of(src, "SubGroups/SubhaloOffsetType", s_orig).astype(np.int64)

    n_arrays = 0
    for t, fields in FIELDS.items():
        for field in fields:
            block = f"{t}/{field}"
            full_cut = np.asarray(cut[block][:])
            for i, g in enumerate(g_orig):
                a, ln = g_start_src[g, t], g_len_src[g, t]
                assert np.array_equal(rows(src, block, a, a + ln), full_cut[g_off_cut[i, t]:g_off_cut[i, t] + ln]), \
                    f"group {g} {block}"
                n_arrays += 1
            for j in range(s_orig.size):
                ln = s_len_cut[j, t]
                if ln == 0:
                    continue
                a = s_off_src[j, t]
                assert np.array_equal(rows(src, block, a, a + ln), full_cut[s_off_cut[j, t]:s_off_cut[j, t] + ln]), \
                    f"subhalo {s_orig[j]} {block}"
                n_arrays += 1
    for block in catalogue_blocks(src):
        if block in RENUMBERED:
            continue
        orig = rows_of(src, block, g_orig if block.startswith("FOFGroups/") else s_orig)
        assert np.array_equal(orig, np.asarray(cut[block][:]), equal_nan=True), block
    assert np.array_equal(np.asarray(cut["FOFGroups/GroupNsubs"][:]),
                          np.bincount(np.asarray(cut["SubGroups/SubhaloGroupNr"][:]), minlength=g_orig.size)), \
        "GroupNsubs != subhalos found inside each group"
    print(f"identical: {n_arrays} particle arrays ({g_orig.size} groups, {s_orig.size} subhalos), "
          f"and every catalogue block not renumbered")


def half_mass_radius(pos, mass, centre, box):
    """Radius enclosing half the mass, about centre, with periodic wrapping (input units)."""
    d = (pos - centre + box / 2) % box - box / 2
    r = np.sqrt((d**2).sum(axis=1))
    o = np.argsort(r)
    cum = np.cumsum(mass[o])
    return np.interp(0.5 * cum[-1], cum, r[o])


def check_units(cut, min_stars=100):
    """Summed particle masses vs SubhaloMassType; particle r_half,* vs SubhaloHalfmassRadType."""
    box = float(cut["Header"].attrs["BoxSize"][0])
    off = np.asarray(cut["SubGroups/SubhaloOffsetType"][:]).astype(np.int64)
    ln = np.asarray(cut["SubGroups/SubhaloLenType"][:]).astype(np.int64)
    mtype = np.asarray(cut["SubGroups/SubhaloMassType"][:]).astype(np.float64)
    rhalf = np.asarray(cut["SubGroups/SubhaloHalfmassRadType"][:]).astype(np.float64)
    spos = np.asarray(cut["SubGroups/SubhaloPos"][:]).astype(np.float64)
    mass = {t: np.asarray(cut[f"{t}/Mass"][:]).astype(np.float64) for t in FIELDS}
    star_pos = np.asarray(cut["4/Position"][:])
    for t, name in [(0, "gas"), (4, "stars"), (5, "BHs")]:
        sel = ln[:, t] > 0
        summed = np.array([mass[t][o:o + n].sum() for o, n in zip(off[sel, t], ln[sel, t])])
        rel = summed / mtype[sel, t] - 1
        print(f"  {name:5s}: summed Mass / SubhaloMassType - 1 over {sel.sum()} subhalos: "
              f"max |.| = {np.abs(rel).max():.2e}")
    sel = np.nonzero(ln[:, 4] >= min_stars)[0]
    ratio = np.array([rhalf[j, 4] / half_mass_radius(star_pos[off[j, 4]:off[j, 4] + ln[j, 4]],
                                                     mass[4][off[j, 4]:off[j, 4] + ln[j, 4]], spos[j], box)
                      for j in sel])
    print(f"  r_half,*: SubhaloHalfmassRadType / particle value (both about SubhaloPos, Position units) "
          f"over {sel.size} subhalos with >= {min_stars} stars: median {np.median(ratio):.4f}, "
          f"range {ratio.min():.4f}-{ratio.max():.4f}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--pig-dir", type=Path, default=config.EXTERNAL["astrid_pig_dir"],
                   help="directory holding PIG_<snap>_subfind (read only)")
    p.add_argument("--snap", default=SNAP)
    p.add_argument("--out-dir", type=Path, default=config.TEST_CUTOUT_PIG_DIR,
                   help="the cutout is written to <out-dir>/PIG_<snap>_subfind")
    p.add_argument("--apply", action="store_true", help="write the cutout (default: dry run)")
    p.add_argument("--validate", action="store_true", help="only validate an existing cutout")
    args = p.parse_args()

    src_path = args.pig_dir / f"PIG_{args.snap}_subfind"
    out_path = args.out_dir / f"PIG_{args.snap}_subfind"
    if args.pig_dir.resolve() in out_path.resolve().parents:
        raise SystemExit(f"refusing to write inside the source directory {args.pig_dir}")
    src = bigfile.File(str(src_path))

    if args.validate:
        cut = bigfile.File(str(out_path))
        check_identical(src, cut)
        check_units(cut)
        return
    groups, sub_runs, g_start, g_len = plan(src, TEST_SUBHALOS)
    print(f"source {src_path}\ncutout {out_path}")
    describe(src, groups, sub_runs, g_len)
    if not args.apply:
        print("dry run: nothing written (use --apply)")
        return
    if out_path.exists():
        raise SystemExit(f"{out_path} exists; delete it by hand to rewrite it")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    groups, _ = write_cutout(src, out_path, TEST_SUBHALOS)
    write_readme(out_path, src_path, TEST_SUBHALOS, groups)
    print("written; validating")
    cut = bigfile.File(str(out_path))
    check_identical(src, cut)
    check_units(cut)


if __name__ == "__main__":
    main()

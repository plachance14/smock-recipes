"""Build synthesizer incident grids of CLOUDY's ``agn`` continuum shape.

The spectra feed syncretize's incident-grid route (``create_cloudy_input_grid.py``
turns each one into a ``table SED``). Each output file has a single axis, the
big-bump temperature, because syncretize runs every incident point against every
photoionisation point. A separate file is written for each global (alpha_ox,
alpha_uv) pair, so the alpha variations stay one-at-a-time.

Functional form, ported from the CLOUDY c23.01 source (``parse_agn.cpp`` and
``ffun1`` in ``cont_ffun.cpp``), with x = h nu / 1 Ryd:

    F_nu = x^a_uv exp(-x/x_bb) exp(-0.01/x)              big bump, IR cutoff at 0.01 Ryd
         + a x^a_x exp(-1/x)            0.1 < x < 7350   X-ray power law (zero below 0.1 Ryd)
         + a 7350^a_x (x/7350)^-2       x >= 7350        above 100 keV

with x_bb = kT/1 Ryd. The constant a is set by alpha_ox using the two components
separately, not the total spectrum: the X-ray term alone (without its exp(-1/x)
factor) at 147 Ryd (2 keV) over the big bump alone (without the IR cutoff and the
X-ray term) at 0.3645 Ryd (2500 A) equals 403.3^alpha_ox.

Each spectrum is normalised to L_bol = 1 erg/s over CLOUDY's energy mesh
(10 MHz to 100 MeV), so the grid is per unit bolometric luminosity
(WeightVariable = "bolometric_luminosities").

Usage (sharable_mock_obs env):
    python agn_grid/make_incident_grid.py --set core --dry-run
    python agn_grid/make_incident_grid.py --set all --apply
    python agn_grid/make_incident_grid.py --check-con path/to/temp_5_0.con
"""

import argparse
from pathlib import Path

import numpy as np
from unyt import Angstrom, Hz, K, c, dimensionless, erg, h, s

# CLOUDY energy mesh limits (mesh.h: 10 MHz to 100 MeV), in Ryd.
RYD_HZ = 3.289841960e15  # 1 Ryd in Hz (infinite-mass Rydberg frequency)
E_MIN_RYD = 1e7 / RYD_HZ
E_MAX_RYD = 1e8 / 13.605693122994  # 100 MeV in Ryd
TE1RYD = 157887.512  # 1 Ryd in K (CLOUDY physconst TE1RYD = h c R_inf / k)
SEXP_LIMIT = 84.0  # CLOUDY's sexp() returns 0 beyond this

# Points per decade of the tabulated SED. 0.005 dex resolves the exp(-x/x_bb)
# cutoff at 10^4 K to well under 1% between points.
DEX_STEP = 0.005

# The 20-point big-bump temperature axis (cloudy_grid_plan.txt).
LOG_T_AXIS = [
    4.00, 4.05, 4.10, 4.15, 4.20, 4.25, 4.30,
    4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00,
    5.25, 5.50, 5.75, 6.00, 6.25, 6.50,
]
ALPHA_X = -1.0  # CLOUDY default; fixed for all grids
FID_AOX, FID_AUV = -1.4, -0.5

# Each set is a list of (name tag, log T values, alpha_ox, alpha_uv).
INCIDENT_SETS = {
    "core": [("T20", LOG_T_AXIS, FID_AOX, FID_AUV)],
    "variations": (
        [("T20", LOG_T_AXIS, aox, FID_AUV) for aox in (-1.0, -1.8, -2.2, -3.0, -5.0)]
        + [("T20", LOG_T_AXIS, FID_AOX, auv) for auv in (-1.0, -0.1, 1 / 3)]
    ),
    # Pilot: T extremes for the corner runs, T = 10^5 K for the fiducials.
    "pilot": [
        ("pilotTcorners", [4.0, 6.5], FID_AOX, FID_AUV),
        ("pilotT5", [5.0], FID_AOX, FID_AUV),
        ("pilotT5", [5.0], -5.0, FID_AUV),
    ],
}


def sexp(x):
    """CLOUDY's safe exp(-x): exp(-x) for x < 84, else exactly 0."""
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    ok = x < SEXP_LIMIT
    out[ok] = np.exp(-x[ok])
    return out


def agn_photon_flux(x, log_T, alpha_ox, alpha_uv, alpha_x=ALPHA_X):
    """CLOUDY's ``agn`` continuum as ffun1 returns it (photons per unit energy).

    Args:
        x: photon energy in Ryd (array).
        log_T: log10 of the big-bump temperature in K.
        alpha_ox, alpha_uv, alpha_x: slopes as on the CLOUDY ``agn`` command.

    Returns:
        Photon flux per unit energy, arbitrary normalisation.
    """
    x = np.asarray(x, dtype=float)
    x_bb = 10.0**log_T / TE1RYD  # big-bump cutoff energy in Ryd

    # Normalisation of the X-ray term from alpha_ox (ParseAgn).
    ratio = 403.3 ** (alpha_ox - 1.0)
    big_bump_2500 = 0.3645 ** (alpha_uv - 1.0) * sexp(0.3645 / x_bb)
    xrays_2kev = 147.0 ** (alpha_x - 1.0)
    a_x = ratio / (xrays_2kev / big_bump_2500)

    f = x ** (alpha_uv - 1.0) * sexp(x / x_bb) * sexp(0.01 / x)
    mid = (x > 0.1) & (x < 7350.0)
    high = x >= 7350.0
    f[mid] += x[mid] ** (alpha_x - 1.0) * a_x * sexp(1.0 / x[mid])
    f[high] += 7350.0 ** (alpha_x - 1.0) * a_x / (x[high] / 7350.0) ** 3
    return f


def wavelength_grid():
    """Log-uniform wavelengths in Angstrom spanning CLOUDY's mesh, ascending."""
    x = 10 ** np.arange(np.log10(E_MIN_RYD), np.log10(E_MAX_RYD), DEX_STEP)
    nu = x * RYD_HZ * Hz
    return np.sort((c / nu).to("Angstrom").value) * Angstrom


def agn_lnu_per_lbol(lam, log_T, alpha_ox, alpha_uv, alpha_x=ALPHA_X):
    """L_nu of the ``agn`` shape normalised to L_bol = 1 erg/s.

    Args:
        lam: wavelengths (unyt, ascending).
        log_T, alpha_ox, alpha_uv, alpha_x: as in agn_photon_flux.

    Returns:
        L_nu in erg/s/Hz per erg/s of bolometric luminosity, on ``lam``.
    """
    nu = (c / lam).to("Hz").value
    x = nu / RYD_HZ
    lnu = x * agn_photon_flux(x, log_T, alpha_ox, alpha_uv, alpha_x)  # F_nu = E * N_E
    order = np.argsort(nu)
    lbol = np.trapezoid(lnu[order], nu[order])
    return lnu / lbol


def grid_name(tag, alpha_ox, alpha_uv):
    """File stem, e.g. agn-cloudy_T20_aox-1.40_auv-0.50."""
    return f"agn-cloudy_{tag}_aox{alpha_ox:+.2f}_auv{alpha_uv:+.2f}".replace("+", "")


def write_incident_grid(path, log_T, alpha_ox, alpha_uv):
    """Write one incident grid with syncretize's GridFile.

    Args:
        path: output .hdf5 path.
        log_T: list of log10 big-bump temperatures in K (the single axis).
        alpha_ox, alpha_uv: global slopes for this file.
    """
    # Imported here so --check-con and --dry-run work without syncretize.
    from synthesizer_grids.grid_io import GridFile

    lam = wavelength_grid()
    spec = np.array([agn_lnu_per_lbol(lam, lt, alpha_ox, alpha_uv) for lt in log_T])

    model = {
        "name": path.stem,
        "type": "agn",
        "family": "cloudy_agn",
        "alpha_ox": alpha_ox,
        "alpha_uv": alpha_uv,
        "alpha_x": ALPHA_X,
        "source": "CLOUDY c23.01 agn command (parse_agn.cpp, cont_ffun.cpp)",
    }
    out = GridFile(str(path))
    out.write_grid_common(
        model=model,
        axes={"big_bump_temperature": 10.0 ** np.array(log_T) * K},
        descriptions={
            "big_bump_temperature": "Big-bump temperature of the CLOUDY agn continuum"
        },
        wavelength=lam,
        log_on_read={"big_bump_temperature": True},
        spectra={"incident": spec * erg / s / Hz},
        weight="bolometric_luminosities",
    )
    # Q(HI), Q(HeII) per unit L_bol; create_synthesizer_grid uses it to
    # rescale CLOUDY's per-area output back to per unit L_bol.
    out.add_specific_ionising_lum()


def check_against_con(con_path, log_T=5.0, alpha_ox=-1.4, alpha_uv=-0.5,
                      alpha_x=-1.0, band=(0.01, 100.0)):
    """Compare our shape with the incident column of an old CLOUDY .con file.

    The .con must be saved with ``units nm`` (as the old decks were). Both
    spectra are normalised by their nu F_nu integral over ``band`` (Ryd).

    Returns:
        (x, ratio) over the band: energies in Ryd and ours / CLOUDY.
    """
    lam_nm, inc = np.loadtxt(con_path, usecols=(0, 1), unpack=True)
    x = 91.1267 / lam_nm  # Ryd: 1 Ryd = 91.1267 nm (infinite-mass)
    keep = (x >= band[0]) & (x <= band[1]) & (inc > 0)
    x, nufnu_cloudy = x[keep], inc[keep]
    nufnu_ours = x**2 * agn_photon_flux(x, log_T, alpha_ox, alpha_uv, alpha_x)

    order = np.argsort(x)
    lnx = np.log(x[order])
    norm_c = np.trapezoid(nufnu_cloudy[order], lnx)
    norm_o = np.trapezoid(nufnu_ours[order], lnx)
    ratio = (nufnu_ours / norm_o) / (nufnu_cloudy / norm_c)
    return x[order], ratio[order]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--set", choices=[*INCIDENT_SETS, "all"], default="all")
    parser.add_argument("--grid-dir", default="cloudy_runs/incident_grids",
                        help="output directory (default is git-ignored)")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", default=True,
                      help="list what would be written (default)")
    mode.add_argument("--apply", action="store_true", help="write the files")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--check-con", metavar="CON_FILE",
                        help="compare with an old .con (agn 5.0 -1.40 -0.5 -1.0) and exit")
    args = parser.parse_args()

    if args.check_con:
        x, ratio = check_against_con(args.check_con)
        dev = np.abs(ratio - 1)
        print(f"{args.check_con}: {x.size} points in 0.01-100 Ryd")
        print(f"  ours/CLOUDY: median {np.median(ratio):.4f}, max |dev| {dev.max():.4f} "
              f"at {x[dev.argmax()]:.3g} Ryd, 95th pct {np.percentile(dev, 95):.4f}")
        for lo, hi in [(0.01, 0.1), (0.1, 1), (1, 10), (10, 100)]:
            sel = (x >= lo) & (x < hi)
            print(f"  {lo:>5}-{hi:<4} Ryd: max |dev| {dev[sel].max():.4f}")
        return

    sets = list(INCIDENT_SETS) if args.set == "all" else [args.set]
    grid_dir = Path(args.grid_dir)
    for set_name in sets:
        for tag, log_T, aox, auv in INCIDENT_SETS[set_name]:
            path = grid_dir / f"{grid_name(tag, aox, auv)}.hdf5"
            exists = path.exists()
            action = "skip (exists)" if exists and not args.overwrite else "write"
            print(f"[{set_name}] {action}: {path}  ({len(log_T)} T values)")
            if args.apply and action == "write":
                grid_dir.mkdir(parents=True, exist_ok=True)
                if exists:
                    path.unlink()
                write_incident_grid(path, log_T, aox, auv)
    if not args.apply:
        print("Dry run: nothing written. Re-run with --apply.")


if __name__ == "__main__":
    main()

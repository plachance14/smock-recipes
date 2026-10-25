"""Make the NIRCam PSFs for the 0.03" image grid with STPSF, one FITS file per filter.

Run with the Python that has stpsf (it needs a different numpy from the main env, so it
lives in a separate virtualenv, with smock installed into it; see README.md), with
STPSF_PATH pointing at its data:

    STPSF_PATH=<stpsf-data> <stpsf-venv>/bin/python psfs/make_psfs.py --out DIR
    ... --overwrite        replace existing files

Most users never run this: the PSFs are synthesizer data, downloaded with the test data.

Each file, NIRCam_<filter>_0.03arcsec.fits, holds two of STPSF's own products, centred on
the central pixel and each normalised to sum 1 over the array:
  DET_SAMP  the PSF on the 0.03" grid, integrated over each pixel (STPSF's detector-sampled
            product). Used directly as synthesizer's PSF, with psf_resample_factor = 1. It is
            exact for a point source at a pixel centre and within 2% of the exact answer for
            compact clumps (see the recipe README).
  OVERSAMP  STPSF's 5x oversampled product (0.006" pixels), for drawing images on a finer grid.
The headers keep STPSF's own (instrument, detector, OPD, wavelengths), plus our settings,
the fraction of the PSF inside the array before normalising, and the encircled energy in an
r = 0.1" aperture (from OVERSAMP), which the depth-to-noise conversion needs.

Choices (author, 14 Oct 2026): 0.03" pixels like the CEERS and JADES mosaics; the
undistorted products, because drizzling removes the distortion; STPSF's default OPD,
detector and source spectrum.
"""
import argparse
import datetime
import os
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

from smock import config  # installed into the STPSF venv with --no-deps (config has no dependencies)

FILTERS = ("F090W", "F115W", "F150W", "F200W", "F277W", "F356W", "F444W")
PIXEL_ARCSEC = 0.03
FOV_PIXELS = 201        # odd, so the PSF centre is the centre of the middle pixel; 6.03"
OVERSAMPLE = 5          # odd, so the PSF centre is the centre of a sub-pixel too
EE_RADIUS_ARCSEC = 0.1  # the survey depth aperture (smock.imaging)


def encircled_energy_fine(oversamp, pixel_arcsec, radius_arcsec):
    """Fraction of the array's flux within radius_arcsec of the centre, on the fine grid."""
    m = oversamp.shape[0]
    xy = (np.arange(m) - (m - 1) / 2) * pixel_arcsec
    r2 = xy[:, None] ** 2 + xy[None, :] ** 2
    return float(oversamp[r2 < radius_arcsec**2].sum() / oversamp.sum())


def make_psf(stpsf, filt):
    nrc = stpsf.NIRCam()
    nrc.filter = filt                     # also picks the SW or LW channel and its default detector
    nrc.pixelscale = PIXEL_ARCSEC
    hdul = nrc.calc_psf(fov_pixels=FOV_PIXELS, oversample=OVERSAMPLE)
    over = hdul["OVERSAMP"].data.astype(np.float64)
    det = hdul["DET_SAMP"].data.astype(np.float64)
    assert det.shape == (FOV_PIXELS, FOV_PIXELS)

    header = hdul["DET_SAMP"].header.copy()
    header.remove("EXTNAME", ignore_missing=True)      # each extension gets its own name below
    for key, value, comment in (
            ("SMOCKPIX", PIXEL_ARCSEC, "[arcsec] pixel scale of DET_SAMP"),
            ("SMOCKOVR", OVERSAMPLE, "oversampling of OVERSAMP"),
            ("FOVSUM", float(det.sum()), "PSF fraction inside the array, before normalising"),
            ("EE_R0P1", encircled_energy_fine(over, PIXEL_ARCSEC / OVERSAMPLE, EE_RADIUS_ARCSEC),
             "encircled energy in r=0.1\" (of the array's flux)"),
            ("STPSFVER", stpsf.__version__, "stpsf version"),
            ("DATE_MK", datetime.date.today().isoformat(), "made by smock-recipes/psfs/make_psfs.py")):
        header[key] = (value, comment)

    out = fits.HDUList([fits.PrimaryHDU(header=header)])
    for name, arr in (("DET_SAMP", det), ("OVERSAMP", over)):
        out.append(fits.ImageHDU((arr / arr.sum()).astype(np.float32), header=header, name=name))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=config.PSF_DIR,
                   help="output directory (default: config.PSF_DIR = %(default)s)")
    p.add_argument("--filters", nargs="+", default=FILTERS)
    p.add_argument("--overwrite", action="store_true", help="replace existing files")
    args = p.parse_args()

    if "STPSF_PATH" not in os.environ:
        sys.exit("Set STPSF_PATH to the STPSF data directory (see the docstring).")
    import stpsf

    args.out.mkdir(parents=True, exist_ok=True)
    for filt in args.filters:
        path = args.out / f"NIRCam_{filt}_{PIXEL_ARCSEC}arcsec.fits"
        if path.exists() and not args.overwrite:
            print(f"{path} exists; skipped (--overwrite to replace)")
            continue
        hdul = make_psf(stpsf, filt)
        hdul.writeto(path, overwrite=args.overwrite)
        h = hdul[0].header
        print(f"{filt}: {path.name}  detector {h.get('DET_NAME')}, in-array fraction {h['FOVSUM']:.4f}, "
              f"EE(0.1\") {h['EE_R0P1']:.3f}", flush=True)


if __name__ == "__main__":
    main()

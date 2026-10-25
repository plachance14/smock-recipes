# Recipe: NIRCam PSFs at 0.03″

**Status:** done. Product made 16 Oct 2026; reproduced bit for bit on 23 Oct 2026 (§5).

## 1. Product
Seven files, `NIRCam_<filter>_0.03arcsec.fits` for F090W, F115W, F150W, F200W, F277W, F356W and
F444W, each with STPSF's `DET_SAMP` (detector-sampled, undistorted, on the 0.03″ grid) and
`OVERSAMP` (5× oversampled) products. They go to `SYNTHESIZER_DIR/instrument_cache/psfs/`
(`smock.config.PSF_DIR`; 30 MB), and smock reads them with `smock.imaging.load_stpsf_psf`.
Users download them with the test data (smock's Google Drive folder), so they never need STPSF.
No raw data kept.

## 2. Requirements
- **STPSF 2.2.0 in its own virtualenv.** stpsf 2.2.0 pins numpy 2.3.5, which would downgrade
  the main env's 2.5.4 and risk synthesizer's compiled extensions. The venv is built from the
  `sharable_mock_obs` Python:

  ```bash
  python -m venv <stpsf-venv>                         # with sharable_mock_obs active
  <stpsf-venv>/bin/python -m pip install stpsf==2.2.0
  <stpsf-venv>/bin/python -m pip install --no-deps -e <path-to>/smock   # for smock.config
  ```

  The last line replaces the old `sys.path` insert (author, 17 Oct 2026): `make_psfs.py`
  imports only `smock.config`, which has no dependencies, so `--no-deps` is enough.
- **STPSF's data package, 2.2.0** (92 MB download), unpacked into `<stpsf-data>`; set
  `STPSF_PATH` to it.
- The command in §3 expects both in `external/` at the top of this repository
  (`external/stpsf-venv/`, `external/stpsf-data/`; not committed). Virtualenvs hard-code their
  absolute path, so build the venv where it will stay rather than moving it.
- Compute: a few minutes on the login node.

## 3. Steps
From the top of this repository:

```bash
STPSF_PATH=external/stpsf-data external/stpsf-venv/bin/python \
    psfs/make_psfs.py --out DIR          # add --overwrite to replace existing files
```

`--out` defaults to `smock.config.PSF_DIR`.

## 4. Problems solved
- **Sampling:** `DET_SAMP` is used directly, with
  `psf_resample_factor = 1`. It is exact for a point source at a pixel centre and within 2% of
  the exact fine-grid answer for compact clumps. synthesizer's own supersampling is worse (up
  to −24% at the peak in F090W). A first proposal, pixel-centre sampling, was withdrawn after
  testing against the exact answer (`archive/2026-10-16_psf_sampling_test/`).
- **numpy conflict:** the separate venv (above).
- **Import of smock.config inside the venv:** smock installed with `--no-deps` (above).

## 5. Verification
Each file's header records `FOVSUM` (the PSF fraction inside the array before normalising)
and `EE_R0P1` (encircled energy in r = 0.1″, used by the depth-to-noise conversion).

**Reproduction (23 Oct 2026):** all 7 filters regenerated with §3's command, in a venv rebuilt
with the commands above (its unpinned dependencies had drifted, e.g. synphot 1.7.0 →
1.8.0), are bit-identical to the installed files in both products, `DET_SAMP` and `OVERSAMP`,
with the same `FOVSUM` and `EE_R0P1` (`archive/2026-10-23_psf_reproduction/compare_psfs.txt`).

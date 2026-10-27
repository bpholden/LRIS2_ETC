# CLAUDE.md

Exposure time calculator (ETC) for LRIS-2, the planned successor to the LRIS spectrograph on Keck I. It computes per-wavelength SNR for a slit spectrum of a template source scaled to a given magnitude. It also compares LRIS-2 against the current LRIS.

## Installing and running

The code is the `lris2_etc` package, installed from [pyproject.toml](pyproject.toml). All data files ship inside the package and are found through `lris2_etc.DATA_DIR`, so the working directory doesn't matter.

```sh
pip install -e .                          # editable install for development
lris2-etc                                 # LRIS-2 blue+red SNR plot, defaults: r'=22, 1200 s, starb1
lris2-etc -m 20 -t 600 -z 0.5 -T elliptical --red_grism R725p2 --blue_grism B1250p2
lris2-etc --flux_plots                    # also show intermediate flux plots at each pipeline stage
lris2-etc -x 1.5                          # extract 1.5 x seeing FWHM along the slit instead of the full slit
lris2-compare                             # LRIS-2 vs. LRIS SNR / flux / sky
lris2-moffat-beta 0.8 0.904               # Moffat beta from FWHM and side of the 50%-energy square
python plot_throughput_ratio.py           # LRIS-2 vs. LRIS throughput and their ratio
```

`lris2-etc`, `lris2-compare` and `lris2-moffat-beta` are console scripts for `lris2_etc.run_lris2:main`, `lris2_etc.compare_old_new:main` and `lris2_etc.moffat_beta:main`; `python -m lris2_etc.run_lris2` also works. The top-level `plot_*.py` scripts are not part of the package and import it with `from lris2_etc import ...`.

Dependencies: numpy ≥ 2.0 (code uses `numpy.trapezoid`, not `trapz`), astropy, scipy, matplotlib. There is no test suite. All output is interactive `plt.show()` windows; nothing is written to disk.

Inside the package, modules import each other relatively (`from . import Sky`). The module names keep their original capitalization (`lris2_etc.ExpCalc`, `lris2_etc.Instrument`, …).

## Architecture

- [ExpCalc.py](lris2_etc/ExpCalc.py): `ExpCalc` is the core pipeline. `compute_spectrum(time, slit_length, slit_width, aperture=None, beta=3)` runs these steps in order. The extraction length is the slit length, or `aperture` × seeing FWHM (capped at the slit length) when `aperture` is given:
  1. Fraction of a Moffat PSF (`Moffat.moffat_frac`, β from `beta`, seeing used as FWHM) inside slit width × extraction length, by exact 2D integration.
  2. Read the template (FITS ext 1, `WAVELENGTH`/`FLUX` columns in erg/s/cm²/Å) and redshift its wavelengths.
  3. Scale to `app_mag` by computing the template's AB mag through the filter (`Mag.compute_ABmag`).
  4. Multiply by exposure time, slit fraction and telescope area, then convert to photons.
  5. Apply atmospheric extinction (`Transmission`) and instrument throughput.
  6. Interpolate the sky onto the template grid and scale it by extraction length and throughput.
  7. Convert source and sky from per Å to per dispersion pixel (× `Ang_per_pix`·`bind`). Compute variance = source + sky + npix·dark·t + npix·RN², where npix is the number of spatial pixels in the extraction length.
  8. If `per_angstrom=True`, divide signal, sky and variance by Å/pix, so SNR per Å = SNR per pixel / √(Å/pix).

  Results land in attributes over `waves[in_band]`: `signal`, `sky_flux`, `noise` (σ), `snr`, all per pixel or per Å as requested. `flux` remains the full-grid source spectrum in photons/Å. The scripts expose the option as `--per_angstrom`/`-A`; `aperture` as `--aperture`/`-x`; `run_lris2` also exposes `beta` as `--beta`/`-b`.
- [Instrument.py](lris2_etc/Instrument.py): `Instrument` is a parameter bag. Its configurator methods each take a `Telescope`:
  - `lris2_red` (R400/R725p1/R725p2) and `lris2_blue` (B600/B1250p1/B1250p2) for LRIS-2. The p1/p2 gratings each cover part of the arm, so they set narrower `BLUE_CUTOFF`/`RED_CUTOFF` values. The cutoffs are set by physical limits other than the efficiency curves. Their old design names, which still appear in the PDR spreadsheets, were R750→R725p1, R700→R725p2, B1300→B1250p1 and B1200→B1250p2.
  - `lris_red` and `lris_blue` for the current LRIS.
  - `deimos`, which is incomplete: it sets no `name`, so `read_throughput` raises.

  Each configurator sets the band cutoffs (`BLUE_CUTOFF`/`RED_CUTOFF`, in Å), the pixel scale, read noise and dark current. `read_throughput` dispatches on `name`:
  - LRIS-2 (`read_lris2_throughput`) builds the total efficiency on a 1 nm grid from `BLUE_CUTOFF` to `RED_CUTOFF`. It multiplies the grating file `<grating>_grating.csv`, the arm's `common_path_file` and `detector_qe_file` (set in `lris2_red`/`lris2_blue`), and the telescope's `throughput_file` (set for Keck I only). The individual curves are kept in `throughput_components`. Each curve goes through `spline_efficiency`: a cubic spline inside the data, continued linearly from the spline's end slope beyond it, and clipped to [0, 1]. To use a measured curve, point the matching attribute at the new file.
  - LRIS and DEIMOS read XIDL `sens_<name>_<grating>.fits.gz` files.
- [Telescope.py](lris2_etc/Telescope.py): collecting area (cm²) and plate scale ("/mm). Build it with `Telescope(name='keck1')`, or with `Telescope()` and then `.keckone()`.
- [Mag.py](lris2_etc/Mag.py): reads a 2-column filter curve (Å, weight) from `data/filters/`, computes λ_eff and the AB magnitude of a spectrum.
- [Sky.py](lris2_etc/Sky.py): sky spectrum from `data/sky/sky.eps_pang_parcsec_onemicron.fits`, a linear WCS in e⁻/s/Å/arcsec². `rescale(inst)` multiplies by slit width and Å/pix.
- [Transmission.py](lris2_etc/Transmission.py): hard-coded extinction curve (mag/airmass), returned as a transmission at the given airmass.
- [Moffat.py](lris2_etc/Moffat.py): Moffat PSF helpers. Only `moffat_frac` is used by the pipeline. `moffat_beta` inverts it for β, given a FWHM and the side of the centred square holding 50% of the energy. `moffat_snr*` and `make_1d_mod` (which prints) are exploratory.
- [depreciated/](depreciated/): old code kept for reference. It is not maintained and still uses the pre-package imports.
- [data/sky/combine_blue_red.py](lris2_etc/data/sky/combine_blue_red.py): one-off script that generated the combined sky FITS. Run it from inside `lris2_etc/data/sky/`. It is not installed with the package.

## Data

All data lives under `lris2_etc/data/`, which `pyproject.toml` ships as package data. A new file in a new subdirectory, or with an extension the `package-data` globs don't cover, also has to be added there, or it will be missing from non-editable installs.

- `templates/`: Kinney-Calzetti galaxy, QSO, SN Ia, PN and Pickles stellar templates (FITS tables).
- `filters/`: SDSS u'g'r'i'z', LSST ugrizy and Bessell V. The LSST files have `#` header lines, which astropy skips.
- `throughput/`: all LRIS-2 curves have wavelength in nm.
  - `<grating>_tot_eff.csv`: the older total-efficiency curves (`wavelength`, `throughput`, including Keck). The code no longer reads them; they are kept for comparison.
  - Component curves split out of the PDR efficiency spreadsheets, to be replaced with measured values over time. The code multiplies these to get the LRIS-2 throughput:
    - `Keck1_telescope.csv`: TELE-M1, TELE-M2 and ADC, 310–1000 nm.
    - `LRIS2_{blue,red}_common_path.csv`: Window, COLL-M1, COLL-M2, DICHROIC, FOLD and CAMERA. The blue file runs to 570 nm. Rows 561–570 were added by hand:
      - DICHROIC = 1 − (red DICHROIC), since in the 550–560 nm overlap the two sum to 1.000 within 0.002.
      - Window, COLL-M1 and COLL-M2 are taken from the red file, since the two arms agree exactly from 550 to 560 nm.
      - FOLD and CAMERA repeat their 560 nm values.
    - `<grating>_grating.csv`: `DE` (the grating efficiency at field centre, the PDR "DE (0')" column) and `PRISMS`.
    - `STA5100_predicted_blue_qe.dat` and `STA4850_ID280_red_qe.dat`: detector QE.

    These cover only 310–570 nm (blue) and 550–1000 nm (red), so cutoffs outside that range rely on extrapolation: B1250p1 starts at 307 nm, and the B600 and B1250p2 grating DE curves stop at 560 nm while those cutoffs run to 570 nm. Multiplying them reproduces the PDR totals.
  - `sens_*.fits.gz`: XIDL sensitivity files for the current LRIS.

## Gotchas

- **Name lookup is a substring match.** Templates and filters are found by `if name in filename` over `os.listdir`, and the first hit wins. `listdir` order is arbitrary, so use unambiguous fragments, for example `sdss_rprime` rather than `r`.
- **In-place mutation.** `compute_spectrum` mutates `self.flux` and `self.sky.spec` in place, so calling it twice on the same `ExpCalc` compounds the scaling. Create a new `ExpCalc` per run. Each `ExpCalc` builds its own `Sky()`.
- **Dark current** (`Instrument.dark`) is in e⁻/pix/s; values are written as `<e-/pix/hr>/3600.`. The LRIS-2 values are the intended ones. The legacy values are not accurate: LRIS red/blue (4 e⁻/pix/hr) are guesses, and DEIMOS (3.5 e⁻/pix/hr) is an average over its 8 detectors.
- **Units.** Wavelengths are in Å throughout, except in the LRIS-2 throughput CSVs. Slit width and length and seeing are in arcsec.
- The SNR is an unweighted (boxcar) sum over the extraction length: the whole slit by default, or `aperture` × FWHM. There is no optimal (profile-weighted) extraction.
- `Sky.rescale` multiplies by `Ang_per_pix` by default (`per_pix=True`). `ExpCalc` calls it with `per_pix=False` and does the per-pixel conversion itself.

## Conventions

Classes are configured by calling methods after construction, not through constructor arguments, as with `Instrument().lris2_red(tel)`. Docstrings mostly just repeat the signature. CLI scripts use `argparse` with short flags (`-m -t -z -a -s -w -l -f -T`) and share one set of defaults: r'=22, 1200 s, airmass 1.2, seeing 1.0", slit 0.7"×8", template starb1.

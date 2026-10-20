# CLAUDE.md

Exposure time calculator (ETC) for LRIS-2, the planned successor to the LRIS spectrograph on Keck I. It computes per-wavelength SNR for a slit spectrum of a template source scaled to a given magnitude. It also compares LRIS-2 against the current LRIS.

## Running

There is no package, build, test suite or requirements file. The modules are flat scripts imported by name, so **run everything from the repo root**: data paths such as `data/templates`, `data/filters` and `data/throughput` are relative to the working directory.

```sh
python run_lris2.py                      # LRIS-2 blue+red SNR plot, defaults: r'=22, 1200 s, starb1
python run_lris2.py -m 20 -t 600 -z 0.5 -T elliptical --red_grism R700 --blue_grism B1200
python run_lris2.py --flux_plots         # also show intermediate flux plots at each pipeline stage
python compare_old_new.py                # LRIS-2 vs. LRIS SNR / flux / sky
python plot_throughput_ratio.py          # LRIS-2 vs. LRIS throughput and their ratio
```

Dependencies: numpy ≥ 2.0 (code uses `numpy.trapezoid`, not `trapz`), astropy, scipy, matplotlib. All output is interactive `plt.show()` windows; nothing is written to disk.

## Architecture

- [ExpCalc.py](ExpCalc.py): `ExpCalc` is the core pipeline. `compute_spectrum(time, slit_length, slit_width)` runs these steps in order:
  1. Slit-loss fraction from a Moffat PSF (`Moffat.moffat_frac`, β=3, seeing used as FWHM).
  2. Read the template (FITS ext 1, `WAVELENGTH`/`FLUX` columns in erg/s/cm²/Å) and redshift its wavelengths.
  3. Scale to `app_mag` by computing the template's AB mag through the filter (`Mag.compute_ABmag`).
  4. Multiply by exposure time, slit fraction and telescope area, then convert to photons.
  5. Apply atmospheric extinction (`Transmission`) and instrument throughput.
  6. Interpolate the sky onto the template grid and scale it by slit length and throughput.
  7. Convert source and sky from per Å to per dispersion pixel (× `Ang_per_pix`·`bind`). Compute variance = source + sky + npix·dark·t + npix·RN², where npix is the number of spatial pixels along the slit.
  8. If `per_angstrom=True`, divide signal, sky and variance by Å/pix, so SNR per Å = SNR per pixel / √(Å/pix).

  Results land in attributes over `waves[in_band]`: `signal`, `sky_flux`, `noise` (σ), `snr`, all per pixel or per Å as requested. `flux` remains the full-grid source spectrum in photons/Å. The scripts expose the option as `--per_angstrom`/`-A`.
- [Instrument.py](Instrument.py): `Instrument` is a parameter bag. Its configurator methods each take a `Telescope`:
  - `lris2_red` (R400/R700/R750) and `lris2_blue` (B600/B1200/B1300) for LRIS-2.
  - `lris_red` and `lris_blue` for the current LRIS.
  - `deimos`, which is incomplete: it sets no `name`, so `read_throughput` raises.

  Each configurator sets the band cutoffs (`BLUE_CUTOFF`/`RED_CUTOFF`, in Å), the pixel scale, read noise and dark current. `read_throughput` dispatches on `name`:
  - LRIS-2 reads `data/throughput/<grating>_tot_eff.csv`, which is in **nm** and converted to Å on load.
  - LRIS and DEIMOS read XIDL `sens_<name>_<grating>.fits.gz` files.
- [Telescope.py](Telescope.py): collecting area (cm²) and plate scale ("/mm). Build it with `Telescope(name='keck1')`, or with `Telescope()` and then `.keckone()`.
- [Mag.py](Mag.py): reads a 2-column filter curve (Å, weight) from `data/filters/`, computes λ_eff and the AB magnitude of a spectrum.
- [Sky.py](Sky.py): sky spectrum from `data/sky/sky.eps_pang_parcsec_onemicron.fits`, a linear WCS in e⁻/s/Å/arcsec². `rescale(inst)` multiplies by slit width and Å/pix.
- [Transmission.py](Transmission.py): hard-coded extinction curve (mag/airmass), returned as a transmission at the given airmass.
- [Moffat.py](Moffat.py): Moffat PSF helpers. Only `moffat_frac` is used by the pipeline. `moffat_snr*` and `make_1d_mod` (which prints) are exploratory.
- [plot_galspec.py](plot_galspec.py): an older, standalone, procedural version of the pipeline that predates `ExpCalc`. It has its own dichroic split at 5500 Å and no slit-loss term. It is not kept in sync with `ExpCalc`.
- [data/sky/combine_blue_red.py](data/sky/combine_blue_red.py): one-off script that generated the combined sky FITS. Run it from inside `data/sky/`.

## Data

- `data/templates/`: Kinney-Calzetti galaxy, QSO, SN Ia, PN and Pickles stellar templates (FITS tables).
- `data/filters/`: SDSS u'g'r'i'z', LSST ugrizy and Bessell V. The LSST files have `#` header lines, which astropy skips.
- `data/throughput/`: LRIS-2 total efficiency CSVs (`wavelength` in nm, `throughput`). The LRIS `sens_*.fits.gz` files exist locally but are **not tracked**, because `.gitignore` excludes `*.gz`. A fresh clone can't run `compare_old_new.py` or `plot_throughput_ratio.py` without them.

## Gotchas

- **Name lookup is a substring match.** Templates and filters are found by `if name in filename` over `os.listdir`, and the first hit wins. `listdir` order is arbitrary, so use unambiguous fragments, for example `sdss_rprime` rather than `r`.
- **Hard-coded sky path.** `Sky` loads from `$HOME/src/LRIS2_ETC/data/sky`, not from a relative path. Running from a checkout anywhere else breaks it.
- **In-place mutation.** `compute_spectrum` mutates `self.flux` and `self.sky.spec` in place, so calling it twice on the same `ExpCalc` compounds the scaling. Create a new `ExpCalc` per run. Each `ExpCalc` builds its own `Sky()`.
- **Dark current** (`Instrument.dark`) is in e⁻/pix/s; values are written as `<e-/pix/hr>/3600.`. The LRIS-2 values are the intended ones. The legacy values are not accurate: LRIS red/blue (4 e⁻/pix/hr) are guesses, and DEIMOS (3.5 e⁻/pix/hr) is an average over its 8 detectors.
- **Units.** Wavelengths are in Å throughout, except in the LRIS-2 throughput CSVs. Slit width and length and seeing are in arcsec.
- The SNR is computed for the whole slit length (all spatial pixels in the slit), not for an optimal or seeing-sized extraction aperture.
- `Sky.rescale` multiplies by `Ang_per_pix` by default (`per_pix=True`). `ExpCalc` calls it with `per_pix=False` and does the per-pixel conversion itself.
- In `lris_blue`, the grating check `grating not in ("600_4000_D560")` tests against a string, not a tuple, so it is a substring check.
- `combine_blue_red.py` reads `bfn` for both the blue and the red spectrum, which looks like a bug.

## Conventions

Classes are configured by calling methods after construction, not through constructor arguments, as with `Instrument().lris2_red(tel)`. Docstrings mostly just repeat the signature. CLI scripts use `argparse` with short flags (`-m -t -z -a -s -w -l -f -T`) and share one set of defaults: r'=22, 1200 s, airmass 1.2, seeing 1.0", slit 0.7"×8", template starb1.

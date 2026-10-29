This will be an ETC for the successor to LRIS

## Installing

```sh
pip install .        # or: pip install -e .  for development
```

## Running

```sh
lris2-etc --help
lris2-etc -m 22 -t 1200 -T starb1 --red_grism R400 --blue_grism B600
lris2-etc -c config.ini 
lris2-compare        # LRIS-2 compared with the current LRIS
```

There is an example config file called `example_config.ini` 

## lris2-etc runtime options

In a config file, use the long option name without the dashes as the key, for example `slit_width = 1.0`. An option given on the command line overrides the value in the config file.

| Option | Default | Description |
|---|---|---|
| `--mag`, `-m` | `22` | AB magnitude of the source in the `--filter` passband. The template is scaled to this. |
| `--filter`, `-f` | `sdss_rprime.dat` | Filter used to scale the template, chosen from `lris2_etc/data/filters/`: SDSS u'g'r'i'z' (`sdss_uprime`, …), LSST ugrizy (`LSST_LSST.z`, …) or `Bessell_V`. |
| `--template`, `-T` | `starb1` | Source spectrum, chosen from `lris2_etc/data/templates/`: galaxies (`starb1`–`starb6`, `elliptical`, `s0`, `sa`, `sb`, `sc`, `bulge`), `qso`, `sn1a`, `pn`, `orion`, and Pickles stars (`A0V`, `G5V`, …). |
| `--redshift`, `-z` | `0` | Redshift applied to the template. |
| `--time`, `-t` | `1200` | Exposure time in seconds. |
| `--airmass`, `-a` | `1.2` | Airmass, used for atmospheric extinction. |
| `--seeing`, `-s` | `1.0` | Seeing FWHM in arcsec. |
| `--beta`, `-b` | `4` | Moffat β of the seeing profile. Lower values put more light in the wings. `lris2-moffat-beta FWHM SIDE` gives β from a FWHM and the side of the square holding 50% of the light. |
| `--slit_width`, `-w` | `0.7` | Slit width in arcsec. |
| `--slit_length`, `-l` | `8` | Slit length in arcsec. |
| `--aperture`, `-x` | none (whole slit) | Extraction length along the slit, in units of the seeing FWHM, for example `1.5`. It is capped at the slit length. If left out, the whole slit length is summed. |
| `--blue_grism` | `B600` | Blue grating: `B600`, `B1250p1` or `B1250p2`. |
| `--red_grism` | `R400` | Red grating: `R400`, `R725p1` or `R725p2`. |
| `--per_angstrom`, `-A` | off | Report signal, noise and SNR per Å rather than per dispersion pixel. |
| `--flux_plots` | off | Also show plots of the source spectrum at the intermediate steps. |
| `--config`, `-c` | none | INI file with a `[lris2-etc]` section that sets any of the options above. |

The `--filter` and `--template` names match any file whose name contains them, and the first match wins. Use a fragment that matches only one file, such as `sdss_rprime` rather than `r`.

As well as the SNR plot, `lris2-etc` prints the median SNR for each arm over the half-maximum range of `--filter` and within ±10 Å of its effective wavelength.


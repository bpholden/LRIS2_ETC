import argparse
import os

import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants

import Telescope
import Instrument
import Sky
import Mag
import Transmission

DICHROIC = 5500
BLUE_CUTOFF = 3100
RED_CUTOFF = 9500

def parse_args():
    parser = argparse.ArgumentParser(description='Plot a galaxy spectrum')
    parser.add_argument('--redshift', '-z', type=float, default=0, help='Redshift')
    parser.add_argument('--time', '-t', type=float, default=1200, help='exposure time (s)')
    parser.add_argument('--airmass', '-a', type=float, default=1.2, help='airmass')
    parser.add_argument('--seeing', '-s', type=float, default=1.0, help='seeing (arcsec)')
    parser.add_argument('--slit_width', '-w', type=float, default=0.7, help='slit width (arcsec)')
    parser.add_argument('--slit_length', '-l', type=float, default=8, help='slit length (arcsec)')
    parser.add_argument('--mag', '-m', type=float, default=22, help='filter magnitude')
    parser.add_argument('--filter', '-f', type=str, default='sdss_rprime.dat', help='filter')
    parser.add_argument('--flux_plots', action='store_true', help='make flux plots')
    parser.add_argument('--template', '-T', type=str, default='starb1_template.fits', help='template file')
    return parser.parse_args()


def photons(wavelength, flux):

    c = astropy.constants.c.value * 1e10 # convert to Angstroms
    flux *= wavelength /(6.626e-27 * c)  # h is in ergs s, c is in Angstroms s^-1

    return flux

def read_template(fn, redshift=0):

    hdus = fits.open(fn)
    dat = hdus[1].data

    waves = dat['WAVELENGTH']
    flux = dat['FLUX']
    waves *= (1 + redshift)
    return waves, flux


def compute_extinction(wave, spec, airmass=1.0):

    trans = Transmission.Transmission(airmass=airmass)
    extinct_interp = np.interp(wave, trans.wave, trans.extinct)

    return spec * extinct_interp

def main():

    keck_1 = Telescope.Telescope(name='keck1')
    args = parse_args()

    lrisred = Instrument.Instrument()
    lrisred.lris2_red(keck_1)
    lrisblue = Instrument.Instrument()
    lrisblue.lris2_blue(keck_1)

    lrisred.swidth = args.slit_width
    lrisblue.swidth = args.slit_width
    lrisred.sheight = args.slit_length
    lrisblue.sheight = args.slit_length

    filen = 'data/templates'
    filen = os.path.join(filen, args.template)
    waves, flux = read_template(filen, redshift=args.redshift)

    mmag = Mag.Mag(args.filter)
    abs_mag = mmag.compute_ABmag(waves, flux)
    print('abs mag = ', abs_mag)
    scale = 10**(0.4*(abs_mag - args.mag))
    flux *= scale

    bsky = Sky.Sky()
    rsky = Sky.Sky()

    flux *= args.time
    bsky.spec *= args.time
    rsky.spec *= args.time

    bsky.rescale(lrisblue)
    rsky.rescale(lrisred)

    if args.flux_plots:
        plt.plot(waves, flux, 'k-')

        plt.xlabel(r'Wavelength ($\AA$)')
        plt.ylabel(r'Intensity ($ergs\ \AA^{-1}\ cm^{-2}$)')
        plt.show()

    flux = flux * keck_1.area
    if args.flux_plots:
        plt.plot(waves, flux, 'k-')
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.ylabel(r'($ergs\ \AA^{-1}$)')
        plt.show()

    flux = photons(waves, flux)
    flux = compute_extinction(waves, flux, airmass=args.airmass)

    flux[waves <= DICHROIC] *= lrisblue.Ang_per_pix
    flux[waves > DICHROIC] *= lrisred.Ang_per_pix

    red_npix = int(args.seeing / lrisred.scale_perp)
    blue_npix = int(args.seeing / lrisblue.scale_perp)

    in_band = (waves > BLUE_CUTOFF) & (waves < RED_CUTOFF)
    rsky_red_band = (rsky.wave > DICHROIC) & (rsky.wave < RED_CUTOFF)
    bsky_blue_band = (bsky.wave > BLUE_CUTOFF) & (bsky.wave <= DICHROIC)

    blue_waves = (waves > BLUE_CUTOFF) & (waves <= DICHROIC)
    red_waves = (waves > DICHROIC) & (waves < RED_CUTOFF)

    sky_bflux = np.interp(waves[blue_waves], bsky.wave[bsky_blue_band], bsky.spec[bsky_blue_band])
    sky_rflux = np.interp(waves[red_waves], rsky.wave[rsky_red_band], rsky.spec[rsky_red_band])

    sky_bflux *= blue_npix
    sky_rflux *= red_npix

    snr = np.zeros_like(waves)
    blue_noise = np.sqrt(flux[blue_waves] + sky_bflux + blue_npix*lrisblue.readnoise**2)
    snr[blue_waves] = flux[blue_waves]/blue_noise

    red_noise = np.sqrt(flux[red_waves] + sky_rflux + red_npix*lrisred.readnoise**2)
    snr[red_waves] = flux[red_waves]/red_noise

    if args.flux_plots:
        plt.plot(waves[in_band], flux[in_band], 'k-')
        plt.plot(waves[blue_waves], sky_bflux, 'b-')
        plt.plot(waves[red_waves], sky_rflux, 'r-')
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.ylabel(r'($\gamma\ pix^{-1}$)')
        plt.show()


    plt.plot(waves[in_band], snr[in_band], 'k-')
    plt.plot(waves[blue_waves], snr[blue_waves], 'b-')
    plt.plot(waves[red_waves], snr[red_waves], 'r-')
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel(r'SNR')
    plt.show()


if __name__ == "__main__":
    main()

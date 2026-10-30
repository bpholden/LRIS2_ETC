import argparse

import matplotlib.pyplot as plt

from . import ExpCalc
from . import Instrument
from . import Telescope


def parse_args():
    '''
    parse_args()
    '''
    parser = argparse.ArgumentParser(description='Plot a galaxy spectrum')
    parser.add_argument('--redshift', '-z', type=float, default=0, help='Redshift (0.0)')
    parser.add_argument('--time', '-t', type=float, default=1200, help='exposure time (1200 s)')
    parser.add_argument('--airmass', '-a', type=float, default=1.2, help='Airmass (1.2)')
    parser.add_argument('--seeing', '-s', type=float, default=1.0, help='seeing (1.0 arcsec)')
    parser.add_argument('--slit_width', '-w', type=float, default=0.7, \
                        help='slit width (0.7 arcsec)')
    parser.add_argument('--slit_length', '-l', type=float, default=8, \
                        help='slit length (8 arcsec)')
    parser.add_argument('--mag', '-m', type=float, default=22, help='filter magnitude')
    parser.add_argument('--filter', '-f', type=str, default='sdss_rprime.dat', \
                        help='Filter file (SDSS r)')
    parser.add_argument('--template', '-T', type=str, default='starb1', \
                        help='Template file (starb1)')
    parser.add_argument('--flux_plots', action='store_true', help='make flux plots')
    parser.add_argument('--per_angstrom', '-A', action='store_true', \
                        help='signal and noise per Angstrom instead of per pixel')
    parser.add_argument('--aperture', '-x', type=float, default=None, \
                        help='extraction length in units of the seeing FWHM (full slit)')
    parser.add_argument('--optimal', '-O', action='store_true', \
                        help='Horne optimal (variance weighted) extraction instead of a sum')
    return parser.parse_args()

def main():

    args = parse_args()

    keck_1 = Telescope.Telescope()
    keck_1.keckone()
    lris2_blue = Instrument.Instrument()
    lris2_blue.lris2_blue(keck_1)
    lris2_red = Instrument.Instrument()
    lris2_red.lris2_red(keck_1)

    lris_blue = Instrument.Instrument()
    lris_blue.lris_blue(keck_1)
    lris_red = Instrument.Instrument()
    lris_red.lris_red(keck_1)

    blue2_exp_calc = ExpCalc.ExpCalc(args.mag, args.filter, args.seeing, args.airmass, \
                                    args.redshift, args.template, lris2_blue, keck_1)
    red2_exp_calc = ExpCalc.ExpCalc(args.mag, args.filter, args.seeing, args.airmass, \
                                   args.redshift, args.template, lris2_red, keck_1)
    blue_exp_calc = ExpCalc.ExpCalc(args.mag, args.filter, args.seeing, args.airmass, \
                                    args.redshift, args.template, lris_blue, keck_1)
    red_exp_calc = ExpCalc.ExpCalc(args.mag, args.filter, args.seeing, args.airmass, \
                                      args.redshift, args.template, lris_red, keck_1)

    if args.flux_plots:
        blue2_exp_calc.flux_plots = True
        red2_exp_calc.flux_plots = True
        blue_exp_calc.flux_plots = True
        red_exp_calc.flux_plots = True

    blue2_exp_calc.compute_spectrum(args.time, args.slit_length,\
                                                args.slit_width, per_angstrom=args.per_angstrom,
                                                aperture=args.aperture,
                                                optimal=args.optimal)
    red2_exp_calc.compute_spectrum(args.time, args.slit_length, \
                                                args.slit_width, per_angstrom=args.per_angstrom,
                                                aperture=args.aperture,
                                                optimal=args.optimal)
    blue_exp_calc.compute_spectrum(args.time, args.slit_length,\
                                                args.slit_width, per_angstrom=args.per_angstrom,
                                                aperture=args.aperture,
                                                optimal=args.optimal)
    red_exp_calc.compute_spectrum(args.time, args.slit_length, \
                                                args.slit_width, per_angstrom=args.per_angstrom,
                                                aperture=args.aperture,
                                                optimal=args.optimal)

    unit = r'\AA^{-1}' if args.per_angstrom else r'pix^{-1}'
    _, ax = plt.subplots(3, 1, layout='constrained', figsize=(12, 6))

    ax[0].plot(blue2_exp_calc.waves[blue2_exp_calc.in_band], blue2_exp_calc.snr, 'b-', label='Blue2')
    ax[0].plot(red2_exp_calc.waves[red2_exp_calc.in_band], red2_exp_calc.snr, 'r-', label='Red2')
    ax[0].plot(blue_exp_calc.waves[blue_exp_calc.in_band], blue_exp_calc.snr, 'b--', label='Blue')
    ax[0].plot(red_exp_calc.waves[red_exp_calc.in_band], red_exp_calc.snr, 'r--', label='Red')
    ax[0].set_xlabel(r'Wavelength ($\AA$)')
    ax[0].set_ylabel(('SNR per Angstrom' if args.per_angstrom else 'SNR per pixel')
                     + (', optimal' if args.optimal else ''))

    ax[1].plot(blue2_exp_calc.waves[blue2_exp_calc.in_band], blue2_exp_calc.signal, 'b-', label='Blue2')
    ax[1].plot(red2_exp_calc.waves[red2_exp_calc.in_band], red2_exp_calc.signal, 'r-', label='Red2')
    ax[1].plot(blue_exp_calc.waves[blue_exp_calc.in_band], blue_exp_calc.signal, 'b--', label='Blue')
    ax[1].plot(red_exp_calc.waves[red_exp_calc.in_band], red_exp_calc.signal, 'r--', label='Red')
    ax[1].set_xlabel(r'Wavelength ($\AA$)')
    ax[1].set_ylabel(rf'Photons ($\gamma\ {unit}$)')
    #ax[1].set_ylim(0, 1.5e3)
    ax[1].set_yscale('log')

    ax[2].plot(blue2_exp_calc.waves[blue2_exp_calc.in_band], blue2_exp_calc.sky_flux, 'b-', label='Blue2')
    ax[2].plot(red2_exp_calc.waves[red2_exp_calc.in_band], red2_exp_calc.sky_flux, 'r-', label='Red2')
    ax[2].plot(blue_exp_calc.waves[blue_exp_calc.in_band], blue_exp_calc.sky_flux, 'b--', label='Blue')
    ax[2].plot(red_exp_calc.waves[red_exp_calc.in_band], red_exp_calc.sky_flux, 'r--', label='Red')
    ax[2].set_xlabel(r'Wavelength ($\AA$)')
    ax[2].set_ylabel(rf'Sky Photons ($\gamma\ {unit}$)')
    #ax[2].set_ylim(0, 1.5e3)
    ax[2].set_yscale('log')

    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()

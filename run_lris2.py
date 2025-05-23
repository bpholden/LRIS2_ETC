'''Run the exposure time calculator'''
import argparse

import matplotlib.pyplot as plt

import ExpCalc
import Instrument
import Telescope

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
    parser.add_argument('--red_grism', type=str, default='R400', help='Red Grism (R400)')
    parser.add_argument('--blue_grism', type=str, default='B600', help='Blue Grism (B600)')
    parser.add_argument('--template', '-T', type=str, default='starb1', \
                        help='Template file (starb1)')
    parser.add_argument('--flux_plots', action='store_true', help='make flux plots')
    return parser.parse_args()

def main(args):
    '''
    Run the exposure time calculator.
    '''
    keck_1 = Telescope.Telescope()
    keck_1.keckone()
    lris2_blue = Instrument.Instrument()
    lris2_blue.lris2_blue(keck_1, grating=args.blue_grism)
    lris2_red = Instrument.Instrument()
    lris2_red.lris2_red(keck_1, grating=args.red_grism)

    blue_exp_calc = ExpCalc.ExpCalc(args.mag, args.filter, args.seeing, args.airmass, \
                                    args.redshift, args.template, lris2_blue, keck_1)
    red_exp_calc = ExpCalc.ExpCalc(args.mag, args.filter, args.seeing, args.airmass, \
                                   args.redshift, args.template, lris2_red, keck_1)

    if args.flux_plots:
        blue_exp_calc.flux_plots = True
        red_exp_calc.flux_plots = True

    blue_exp_calc.compute_spectrum(args.time, args.slit_length,\
                                                          args.slit_width)
    red_exp_calc.compute_spectrum(args.time, args.slit_length, \
                                                      args.slit_width)

    _, _ = plt.subplots(figsize=(12, 6))
    plt.plot(blue_exp_calc.waves[blue_exp_calc.in_band], blue_exp_calc.snr, 'b-', label='Blue')
    plt.plot(red_exp_calc.waves[red_exp_calc.in_band], red_exp_calc.snr, 'r-', label='Red')
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('SNR')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    args = parse_args()
    main(args)

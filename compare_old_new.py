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
    parser.add_argument('--template', '-T', type=str, default='starb1_template.fits', \
                        help='Template file')
    parser.add_argument('--flux_plots', action='store_true', help='make flux plots')
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

    blue2_snr, blue2_good = blue2_exp_calc.compute_spectrum(args.time, args.slit_length,\
                                                            args.slit_width)
    red2_snr, red2_good = red2_exp_calc.compute_spectrum(args.time, args.slit_length, \
                                                        args.slit_width)
    blue_snr, blue_good = blue_exp_calc.compute_spectrum(args.time, args.slit_length,\
                                                            args.slit_width)
    red_snr, red_good = red_exp_calc.compute_spectrum(args.time, args.slit_length, \
                                                        args.slit_width)

    _, _ = plt.subplots(figsize=(12, 6))

    plt.plot(blue2_exp_calc.waves[blue2_good], blue2_snr, 'b-', label='Blue2')
    plt.plot(red2_exp_calc.waves[red2_good], red2_snr, 'r-', label='Red2')
    plt.plot(blue_exp_calc.waves[blue_good], blue_snr, 'b.-', label='Blue')
    plt.plot(red_exp_calc.waves[red_good], red_snr, 'r.-', label='Red')
    plt.legend()
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('SNR')
    plt.show()

if __name__ == '__main__':
    main()

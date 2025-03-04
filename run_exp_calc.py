'''Run the exposure time calculator'''
import argparse

import ExpCalc
import Instrument
import Telescope

def parse_args():

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
    parser.add_argument('--template', '-T', type=str, default='starb1_template.fits', \
                        help='Template file')
    return parser.parse_args()

def main():

    args = parse_args()

    keck_1 = Telescope.Telescope()
    keck_1.keckone()
    lris2_blue = Instrument.Instrument()
    lris2_blue.lris2_blue(grating=args.blue_grism)
    lris2_red = Instrument.Instrument()
    lris2_red.lris2_red(grating=args.red_grism)

'''Moffat beta from a FWHM and the side of the square holding 50% of the energy'''
import argparse

from . import Moffat

def parse_args():
    '''
    parse_args()
    '''
    parser = argparse.ArgumentParser(description='Moffat beta from FWHM and 50% energy square')
    parser.add_argument('fwhm', type=float, help='FWHM (arcsec)')
    parser.add_argument('side', type=float, \
                        help='side of the centred square holding 50%% of the energy (arcsec)')
    return parser

def main():
    '''
    Print the Moffat beta.
    '''
    parser = parse_args()
    args = parser.parse_args()
    try:
        beta = Moffat.moffat_beta(args.fwhm, args.side)
    except ValueError as err:
        parser.error(str(err))
    print(f'beta = {beta:.4f}')

if __name__ == '__main__':
    main()

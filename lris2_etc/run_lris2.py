'''Run the exposure time calculator'''
import argparse
import configparser

import matplotlib.pyplot as plt
import numpy as np

from . import ExpCalc
from . import Instrument
from . import Telescope

CONFIG_SECTION = 'lris2-etc'

def read_config(parser, filename):
    '''
    read_config(parser, filename)

    Reads the [lris2-etc] section of an INI file. Keys are the long option
    names, and values are converted with the option's own type.
    '''
    config = configparser.ConfigParser()
    if not config.read(filename):
        parser.error(f'cannot read config file {filename}')
    if CONFIG_SECTION not in config:
        parser.error(f'config file {filename} has no [{CONFIG_SECTION}] section')
    actions = {action.dest: action for action in parser._actions}
    values = {}
    for key, value in config[CONFIG_SECTION].items():
        action = actions.get(key)
        if action is None or key in ('help', 'config'):
            parser.error(f'unknown option {key} in config file {filename}')
        try:
            if action.nargs == 0: # store_true flags
                values[key] = config[CONFIG_SECTION].getboolean(key)
            else:
                values[key] = action.type(value) if action.type else value
        except ValueError:
            parser.error(f'bad value {value!r} for {key} in config file {filename}')
    return values

def parse_args():
    '''
    parse_args()
    '''
    parser = argparse.ArgumentParser(description='Plot a galaxy spectrum')
    parser.add_argument('--redshift', '-z', type=float, default=0, help='Redshift (0.0)')
    parser.add_argument('--time', '-t', type=float, default=1200, help='exposure time (1200 s)')
    parser.add_argument('--airmass', '-a', type=float, default=1.2, help='Airmass (1.2)')
    parser.add_argument('--seeing', '-s', type=float, default=1.0, help='seeing (1.0 arcsec)')
    parser.add_argument('--beta', '-b', type=float, default=3, help='Moffat PSF beta (3)')
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
    parser.add_argument('--per_angstrom', '-A', action='store_true', \
                        help='signal and noise per Angstrom instead of per pixel')
    parser.add_argument('--aperture', '-x', type=float, default=None, \
                        help='extraction length in units of the seeing FWHM (full slit)')
    parser.add_argument('--config', '-c', type=str, default=None, \
                        help=f'INI file with a [{CONFIG_SECTION}] section; '
                        'command line options override it')

    # config file values replace the defaults, so the command line still wins
    args, _ = parser.parse_known_args()
    config = {}
    if args.config:
        config = read_config(parser, args.config)
        parser.set_defaults(**config)
    args = parser.parse_args()
    for key, value in config.items():
        if getattr(args, key) != value:
            print(f'Note: command line {key} = {getattr(args, key)} '
                  f'overrides {value} from {args.config}')
    return args

def median_snr(exp_calc, wave_min, wave_max):
    '''
    median_snr(exp_calc, wave_min, wave_max)

    Median SNR between wave_min and wave_max (Angstroms), or None if the
    arm has no points there.
    '''
    waves = exp_calc.waves[exp_calc.in_band]
    sel = (waves >= wave_min) & (waves <= wave_max)
    if not sel.any():
        return None
    return np.median(exp_calc.snr[sel])

def print_median_snr(label, exp_calcs, wave_min, wave_max):
    '''
    print_median_snr(label, exp_calcs, wave_min, wave_max)
    '''
    out = []
    for arm, exp_calc in exp_calcs.items():
        snr = median_snr(exp_calc, wave_min, wave_max)
        out.append(f'{arm} ' + ('n/a' if snr is None else f'{snr:.2f}'))
    print(f'{label}: ' + ', '.join(out))

def main(args=None):
    '''
    Run the exposure time calculator.
    '''
    if args is None:
        args = parse_args()
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
                                                          args.slit_width, per_angstrom=args.per_angstrom,
                                                          aperture=args.aperture, beta=args.beta)
    red_exp_calc.compute_spectrum(args.time, args.slit_length, \
                                                      args.slit_width, per_angstrom=args.per_angstrom,
                                                      aperture=args.aperture, beta=args.beta)

    # the filter used to scale the template: its half-maximum range and its centre
    mag = red_exp_calc.mag
    in_filter = mag.weight >= 0.5*mag.weight.max()
    unit = 'per Angstrom' if args.per_angstrom else 'per pixel'
    exp_calcs = {'Blue': blue_exp_calc, 'Red': red_exp_calc}
    print_median_snr(f'Median SNR {unit} over {args.filter} '
                     f'({mag.wave[in_filter].min():.0f}-{mag.wave[in_filter].max():.0f} A)',
                     exp_calcs, mag.wave[in_filter].min(), mag.wave[in_filter].max())
    print_median_snr(f'Median SNR {unit} at {args.filter} centre ({mag.lambda_eff:.0f} +/- 10 A)',
                     exp_calcs, mag.lambda_eff - 10, mag.lambda_eff + 10)

    _, _ = plt.subplots(figsize=(12, 6))
    plt.plot(blue_exp_calc.waves[blue_exp_calc.in_band], blue_exp_calc.snr, 'b-', label='Blue')
    plt.plot(red_exp_calc.waves[red_exp_calc.in_band], red_exp_calc.snr, 'r-', label='Red')
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('SNR per Angstrom' if args.per_angstrom else 'SNR per pixel')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()

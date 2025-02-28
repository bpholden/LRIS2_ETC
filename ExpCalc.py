import numpy as np
import astropy.constants

import Telescope
import Instrument
import Sky
import Mag
import Transmission


class ExpCalc():
    """
    ExpCalc()
    """
    def __init__(self):
        self.telescope = Telescope.Telescope()
        self.instrument = Instrument.Instrument()
        self.sky = Sky.Sky()
        self.mag = Mag.Mag()
        self.transmission = Transmission.Transmission()
        self.flux = None
        self.wave = None

    def photons(self):
        """
        photons(wavelength, flux)
        """
        c = astropy.constants.c.value * 1e10 # convert to Angstroms
        self.flux *= self.wavelength /(6.626e-27 * c)  # h is in ergs s, c is in Angstroms s^-1

        return 
    
    def compute_extinction(self, airmass=1.0):
        """
        compute_extinction(airmass=1.0)
        """

        trans = Transmission.Transmission(airmass=airmass)
        extinct_interp = np.interp(self.wave, trans.wave, trans.extinct)

        return self.flux * extinct_interp
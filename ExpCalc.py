import numpy as np
import astropy.constants
import astropy.io.fits as fits

import Telescope
import Instrument
import Sky
import Mag
import Transmission


class ExpCalc():
    """
    ExpCalc()
    """
    def __init__(self, app_mag, mfilter, seeing, airmass, redshift, template_filename, inst=None, telescope=None):
        self.telescope = telescope
        self.instrument = inst
        self.sky = Sky.Sky()
        self.mag = None
        self.app_mag = app_mag
        self.abs_mag = None
        self.transmission = Transmission.Transmission()
        self.flux = None
        self.wave = None
        self.airmass = airmass
        self.seeing = seeing
        self.filter = mfilter
        self.redshift = redshift
        self.template_filename = template_filename



    def photons(self):
        """
        photons(self)
        """
        c = astropy.constants.c.value * 1e10 # convert to Angstroms
        self.flux *= self.wave /(6.626e-27 * c)  # h is in ergs s, c is in Angstroms s^-1

        return

    def compute_extinction(self):
        """
        compute_extinction(self)
        """

        trans = Transmission.Transmission(airmass=self.airmass)
        extinct_interp = np.interp(self.wave, trans.wave, trans.extinct)
        self.flux *= extinct_interp
        return

    def read_template(self):
        """
        read_template(self)
        """
        hdus = fits.open(self.template_filename)
        dat = hdus[1].data

        self.waves = dat['WAVELENGTH']
        self.flux = dat['FLUX']
        self.waves *= (1 + self.redshift)
        return

    def compute_spectrum(self, time, slit_length, slit_width, inst=None, telescope=None):
        '''
        compute_spectrum(self, time, slit_length, slit_width, inst=None, telescope=None)
        '''

        if telescope:
            self.telescope = telescope
        if inst:
            self.instrument = inst
        self.instrument.slit_length = slit_length
        self.instrument.slit_width = slit_width

        self.mag = Mag.Mag(self.filter)
        self.read_template()
        self.abs_mag = self.mag.compute_ABmag(self.waves, self.flux)
        self.flux *= 10**(0.4*(self.abs_mag - self.mag))

        self.flux *= time
        self.sky.spec *= time

        self.sky.wave, self.sky.spec = Sky.rescale(self.instrument, \
                                                     self.sky.wave, self.sky.spec)
        self.flux = self.flux * self.telescope.area

        self.photons()
        self.compute_extinction()

        self.flux[self.waves] *= self.instrument.Ang_per_pix

        npix = int(self.seeing / self.instrument.scale_perp)

        in_band = (self.waves > self.instrument.BLUE_CUTOFF) &\
              (self.waves < self.instrument.RED_CUTOFF)
        sky_band = (self.sky.wave > self.instrument.BLUE_CUTOFF) &\
              (self.sky.wave < self.instrument.RED_CUTOFF)

        sky_flux = np.interp(self.waves[in_band], self.sky.wave[sky_band], self.sky.spec[sky_band])
        sky_flux *= npix

        snr = np.zeros_like(self.waves)
        snr[in_band] = self.flux[in_band]/np.sqrt(self.flux[in_band] + sky_flux)

        return snr[in_band]
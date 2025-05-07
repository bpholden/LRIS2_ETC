import os
import numpy as np
import astropy.constants
import astropy.io.fits as fits
import matplotlib.pyplot as plt

import Sky
import Mag
import Transmission
import Moffat

class ExpCalc():
    """
    ExpCalc()
    """
    def __init__(self, app_mag, mfilter, seeing, airmass, redshift, template_filename, inst=None, telescope=None):
        self.telescope = telescope
        self.instrument = inst
        self.sky = Sky.Sky()
        self.airmass = airmass
        self.seeing = seeing
        self.filter = mfilter
        self.redshift = redshift
        self.template_filename = template_filename
        self.app_mag = app_mag

        self.mag = None
        self.abs_mag = None
        self.transmission = Transmission.Transmission()
        self.flux = None
        self.waves = None
        self.good_waves = None
        self.throughput_interp = None
        self.noise = None
        self.snr = None
        self.sky_flux = None
        self.in_band = None
        self.flux_plots = False


    def photons(self):
        """
        photons(self)
        """
        c = astropy.constants.c.value * 1e10 # convert to Angstroms
        self.flux *= self.waves /(6.626e-27 * c)  # h is in ergs s, c is in Angstroms s^-1

        return

    def compute_extinction(self):
        """
        compute_extinction(self)
        """

        trans = Transmission.Transmission(airmass=self.airmass)
        extinct_interp = np.interp(self.waves, trans.wave, trans.extinct)
        self.flux *= extinct_interp
        return

    def compute_throughput(self):
        """
        compute_throughput(self)
        """
        self.instrument.read_throughput()
        self.good_waves = (self.waves > self.instrument.BLUE_CUTOFF) &\
                        (self.waves < self.instrument.RED_CUTOFF)
        self.throughput_interp = np.interp(self.waves[self.good_waves], self.instrument.throughput['wavelength'],\
                                       self.instrument.throughput['throughput'])
        self.flux[self.good_waves] *= self.throughput_interp
        return

    def scale_sky_by_throughput(self, sky_waves, sky_photons):
        """
        scale_sky_by_throughput(self, sky_waves, sky_photons)
        """
        good_waves = (sky_waves > self.instrument.BLUE_CUTOFF) &\
                        (sky_waves < self.instrument.RED_CUTOFF)
        throughput_interp = np.interp(sky_waves[good_waves], self.instrument.throughput['wavelength'],\
                                       self.instrument.throughput['throughput'])
        sky_photons *= throughput_interp
        return sky_photons

    def read_template(self):
        """
        read_template(self)
        """
        dirn = 'data/templates'
        template_list = os.listdir(dirn)
        for template_name in template_list:
            if self.template_filename in template_name:
                self.template_filename = template_name
                break
        if self.template_filename == '':
            print('No template found')
            return
        filen = os.path.join(dirn, self.template_filename)
        hdus = fits.open(filen)
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
        self.instrument.sheight = slit_length
        self.instrument.swidth = slit_width
        frac = Moffat.moffat_frac(self.seeing, slit_width, slit_length, \
                                  pix_size=self.instrument.scale_perp)

        self.mag = Mag.Mag(self.filter)
        self.read_template()
        self.abs_mag = self.mag.compute_ABmag(self.waves, self.flux)
        self.flux *= 10**(0.4*(self.abs_mag - self.app_mag))

        self.flux *= time
        self.flux *= frac
        self.sky.spec *= time

        self.sky.rescale(self.instrument)
        self.flux = self.flux * self.telescope.area

        if self.flux_plots:
            plt.plot(self.waves, self.flux, 'k-')

            plt.xlabel(r'Wavelength ($\AA$)')
            plt.ylabel(r'Intensity ($ergs\ \AA^{-1}\ cm^{-2}$)')
            plt.title(f'Flux {self.instrument.name}')
            plt.show()

        #self.flux *= self.instrument.Ang_per_pix
        self.photons()
        if self.flux_plots:
            plt.plot(self.waves, self.flux, 'k-')
            plt.xlabel(r'Wavelength ($\AA$)')
            plt.ylabel(r'($\gamma\ \AA^{-1}$)')
            plt.title(f'Photons {self.instrument.name}')
            plt.show()
        self.compute_extinction()
        self.compute_throughput()
        if self.flux_plots:
            plt.plot(self.waves, self.flux, 'k-')
            plt.xlabel(r'Wavelength ($\AA$)')
            plt.ylabel(r'($\gamma\ \AA^{-1}$)')
            plt.title(f'Photons after throughput and extinction {self.instrument.name}')
            plt.show()

        npix = int(slit_length )
        if npix < 2:
            npix = 2

        self.in_band = (self.waves > self.instrument.BLUE_CUTOFF) &\
              (self.waves < self.instrument.RED_CUTOFF)
        sky_band = (self.sky.wave > self.instrument.BLUE_CUTOFF) &\
              (self.sky.wave < self.instrument.RED_CUTOFF)

        self.sky_flux = np.interp(self.waves[self.in_band], self.sky.wave[sky_band], \
                                  self.sky.spec[sky_band])
        self.sky_flux *= slit_length

        self.sky_flux = self.scale_sky_by_throughput(self.waves[self.in_band], self.sky_flux)

        self.noise = self.flux[self.in_band]
        self.noise += self.sky_flux
        self.noise += npix*self.instrument.dark*time
        self.noise += npix*self.instrument.readnoise**2

        self.snr = np.zeros_like(self.waves[self.in_band])
        self.snr = self.flux[self.in_band]/np.sqrt(self.noise)

        return 

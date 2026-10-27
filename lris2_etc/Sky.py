import os
import numpy as np
import astropy.constants
import astropy.io.fits

from . import DATA_DIR
from . import Mag
from . import Telescope

# Mauna Kea dark sky, AB mag/arcsec^2. LSST z, not SDSS z', which runs past
# the end of the sky spectrum
SKY_MAGS = {'sdss_uprime': 22.9, 'sdss_gprime': 22.0, 'sdss_rprime': 21.4,
            'sdss_iprime': 20.4, 'LSST_LSST.z': 19.4}

class Sky:

    def __init__(self):

        self.sky_dir = os.path.join(DATA_DIR, 'sky')
        #self.bfn = os.path.join(self.sky_dir,'bsky.eps_pang_parcsec.fits')
        #self.rfn = os.path.join(self.sky_dir,'rsky.eps_pang_parcsec_onemicron.fits')
        self.fn = os.path.join(self.sky_dir,'sky.eps_pang_parcsec_onemicron.fits')
        # note this in e-/s/Ang/arcsec^2
    
        #self.bpix, self.bwave, self.bspec = self.read_spec(self.bfn)
        self.pix, self.wave, self.spec = self.read_spec(self.fn)
        self.scale_to_mags()

        #self.wave = np.concatenate((self.bwave, self.rwave))
        #self.spec = np.concatenate((self.bspec, self.rspec))

    def read_spec(self, fn):

        spec_hdu = astropy.io.fits.open(fn)
        hdr = spec_hdu[0].header
        spec = spec_hdu[0].data
        pix = np.arange(0,len(spec))
        wave =  hdr['CRVAL1'] + pix*hdr['CDELT1']

        return pix, wave, spec


    def scale_to_mags(self):
        '''
        scale_to_mags(self)

        The sky was observed with LRIS on Keck and may still include the
        telescope + instrument throughput. Scale it to SKY_MAGS, with the
        correction interpolated in wavelength between the filters.
        '''
        h = astropy.constants.h.cgs.value
        c = astropy.constants.c.value * 1e10 # Angstroms s^-1
        tel = Telescope.Telescope(name='keck1')
        flam = self.spec * h * c / self.wave / tel.area # ergs/s/cm^2/Ang/arcsec^2

        lam_eff = []
        corr = []
        for filt, target in SKY_MAGS.items():
            mag = Mag.Mag(filt)
            weight = np.interp(self.wave, mag.wave, mag.weight, left=0, right=0)
            # photon-weighted f_nu over the part of the filter the sky covers
            fnu = np.trapezoid(flam * self.wave * weight, x=self.wave)
            fnu /= np.trapezoid(c / self.wave * weight, x=self.wave)
            sky_mag = -2.5 * np.log10(fnu) - 48.6
            lam_eff.append(np.trapezoid(self.wave * weight, x=self.wave) /
                           np.trapezoid(weight, x=self.wave))
            corr.append(10**(-0.4 * (target - sky_mag)))

        self.spec *= np.interp(self.wave, lam_eff, corr)
        return

    def rescale(self, instrument, per_pix=True):

        #instrument.Ang_per_pix *= instrument.pixel_size # dwave is for a 1" slit
        #sq_arcsec = instrument.swidth * instrument.sheight
        sq_arcsec = instrument.swidth #* instrument.scale_perp
        self.spec *= sq_arcsec
        # this should be e-/s/Ang
        if per_pix:
            self.spec *= instrument.Ang_per_pix
            # this should be e-/pix/s

        return

def compute_dwave(wave, R):

    dwave = wave / R

    return dwave

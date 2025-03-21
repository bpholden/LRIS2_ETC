import os
import numpy as np
import astropy.io.fits

class Sky:

    def __init__(self):

        self.sky_dir = 'data/sky'
        #self.bfn = os.path.join(self.sky_dir,'bsky.eps_pang_parcsec.fits')
        #self.rfn = os.path.join(self.sky_dir,'rsky.eps_pang_parcsec_onemicron.fits')
        self.fn = os.path.join(self.sky_dir,'sky.eps_pang_parcsec_onemicron.fits')

        #self.bpix, self.bwave, self.bspec = self.read_spec(self.bfn)
        self.pix, self.wave, self.spec = self.read_spec(self.fn)

        #self.wave = np.concatenate((self.bwave, self.rwave))
        #self.spec = np.concatenate((self.bspec, self.rspec))

    def read_spec(self, fn):

        spec_hdu = astropy.io.fits.open(fn)
        hdr = spec_hdu[0].header
        spec = spec_hdu[0].data
        pix = np.arange(0,len(spec))
        wave =  hdr['CRVAL1'] + pix*hdr['CDELT1']

        return pix, wave, spec


def rescale(instrument, wave, spec):

    #instrument.Ang_per_pix *= instrument.pixel_size # dwave is for a 1" slit
    #sq_arcsec = instrument.swidth * instrument.sheight
    sq_arcsec = instrument.swidth * instrument.scale_perp
    spec *= sq_arcsec
    # this should be e-/s/Ang
    #spec *= instrument.Ang_per_pix
    # this should be e-/pix/s

    return wave, spec

def compute_dwave(wave, R):

    dwave = wave / R

    return dwave

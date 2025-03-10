import os
import numpy
import astropy.constants
import astropy.io.ascii

class Mag:
    def __init__(self, fn=None):
        self.fn = fn
        self.wave = None
        self.weight = None
        self.lambda_eff = None
        self.mag = None
        self.filter_dir = os.path.join('data', 'filters')
        if self.fn is not None:
            self.read_filter(fn)
            self.lambda_eff = self.compute_lambda_eff()

    def __repr__(self):
        return f'<Mag {self.fn} {self.lambda_eff:0.2f}>'
    
    def __str__(self):
        return f'<Mag {self.fn} {self.lambda_eff:0.2f}>'

    def read_filter(self, fn):
        '''
        read in a filter file

        the filter file should have two columns:
        wavelength (should be in the Angstroms)
        weight (the weight of the filter at that wavelength)
        '''
        self.fn = fn
        fn = os.path.join(self.filter_dir, fn)
        data = astropy.io.ascii.read(fn)
        self.wave = data['col1']
        self.weight = data['col2']
        self.lambda_eff = self.compute_lambda_eff()

    def compute_lambda_eff(self):
        '''
        compute the effective wavelength of a filter
        
        returns the effective wavelength in the units 
        of the input filter file, but the units should
        be Angstroms
        '''
        wv= numpy.trapz(self.wave, self.wave * self.weight)
        wv = wv / numpy.trapz(self.wave, self.weight)
        return wv

    def compute_ABmag(self, spec_wave, spec):
        '''
        compute the magnitude of a spectrum

        spec_wave is in Angstroms
        spec is in ergs/s/cm^2/Angstrom
        '''
        weights_interp = numpy.interp(spec_wave, self.wave, self.weight)
        tot_filt = numpy.trapz(spec_wave, weights_interp)
        tot_flux = numpy.trapz(spec_wave, weights_interp * spec)
        tot_flux = tot_flux / tot_filt
        tot_flux = tot_flux / astropy.constants.c.value
        tot_flux = tot_flux * 1e-10 # convert c to Angstroms
        tot_flux = tot_flux * self.lambda_eff**2
        self.mag = -2.5 * numpy.log10(tot_flux) - 48.6
        return self.mag

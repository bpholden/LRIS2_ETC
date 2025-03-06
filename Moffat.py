# Moffat
import astropy.modeling
import numpy
from scipy import integrate


def prep_values(scale_size, dim_size, pix_size):
    '''
    prep_values
    '''

    sc_size = scale_size / pix_size

    nsamps = int(dim_size / pix_size)
    if nsamps % 2 == 0:
        nsamps += 1
    samps = numpy.linspace(-dim_size/(2*pix_size),  dim_size/(2*pix_size), nsamps)

    xs, ys = numpy.meshgrid(samps, samps)

    return xs, ys, sc_size

def make_1d_mod(size, beta, amp=1):
    '''
    make_mod(size, beta, amp=1)

    Makes the corrected 1D Moffat function.
    The size is assuemd to be a FWHM, turns that into a Moffat
    scale size
    '''
    moffat_size = size / numpy.sqrt(2**(1/beta) - 1)
    moffat_size /= 2
    print(f"size scale = {moffat_size:0.3f}")
    amp_norm = (beta - 1)/  (numpy.pi*moffat_size **2)
    amp *= amp_norm
    print(f"scale of amp = {amp_norm:0.4f}")
    moffat_mod = astropy.modeling.functional_models.Moffat1D(amplitude=amp,
                                                         gamma=moffat_size, alpha=beta)

    return moffat_mod

def make_2d_mod(size, beta, amp=1):
    '''
    make_mod(size, beta, amp=1)
    Makes a 2D Moffat
    '''
    moffat_size = size / numpy.sqrt(2**(1/beta) - 1)
    moffat_size /= 2
    amp_norm = (beta - 1)/  (numpy.pi * moffat_size **2)
    amp *= amp_norm

    moffat_mod = astropy.modeling.functional_models.Moffat2D(amplitude=amp,
                                                         gamma=moffat_size, alpha=beta)

    return moffat_mod

def moffat_snr(flux, size, beta=3, width=0.75, height=8., pix_size=0.15):
    '''
    moffat_snr
    '''
    xs, ys, sc_size = prep_values(size, height, pix_size)
    width /= pix_size
    height /= pix_size

    moffat_mod = make_2d_mod(sc_size, beta, flux)

    background_flux = flux / (1.0/pix_size)**2
    background_flux *= 10

    background = astropy.modeling.models.Const2D(amplitude=background_flux)

    source_flux = integrate.dblquad(moffat_mod, -width/2, width/2, -height/2, height/2)
    sum_background_flux = integrate.dblquad(background, -width/2, width/2, -height/2, height/2)
    snr = source_flux[0] / numpy.sqrt(sum_background_flux[0])
    return snr

def moffat_snr_optimal(flux, size, beta=3, width=0.75, height=8., pix_size=0.15):
    '''
    moffat_snr_optimal
    '''

    xs, ys, sc_size = prep_values(size, height, pix_size)
    width /= pix_size
    height /= pix_size

    slit = astropy.modeling.functional_models.Box2D(1, x_width=width, y_width=height)
    moffat_profile = make_2d_mod(sc_size, beta, 1)*slit
    background_flux = 100*flux / (1.0/pix_size)**2
    background_mod = astropy.modeling.models.Const2D(amplitude=background_flux)*slit

    moffat_mod_eval = moffat_profile(xs, ys)
    background_mod_eval = background_mod(xs, ys)

    profile = numpy.trapz(moffat_mod_eval, axis=1)
    background_profile = numpy.trapz(background_mod_eval, axis=1)

    weighted_spectrum = flux*profile*profile  / background_profile
    weighted_noise = profile*profile / background_profile

    numerator = numpy.sum(weighted_spectrum)
    denominator = numpy.sum(weighted_noise)

    source = numerator/denominator
    var = 1/denominator
    snr = source / numpy.sqrt(var)
    return snr

def moffat_frac(size, width, height, beta=3, pix_size=0.15):
    '''
    moffat_frac
    '''
    xs, ys, sc_size = prep_values(size, height, pix_size)
    width /= pix_size
    height /= pix_size

    slit = astropy.modeling.functional_models.Box2D(1, x_width=width, y_width=height)
    moffat_profile = make_2d_mod(sc_size, beta, 1)*slit
    moffat_mod_eval = moffat_profile(xs, ys)
    frac = numpy.sum(moffat_mod_eval)
    return frac
  
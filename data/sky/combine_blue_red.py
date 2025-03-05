import numpy as np
import astropy.io.fits

def read_file(fn):
    spec_hdu = astropy.io.fits.open(fn)
    spec = spec_hdu[0].data
    hdr = spec_hdu[0].header
    pix = np.arange(0,len(spec))
    wave =  hdr['CRVAL1'] + pix*hdr['CDELT1']

    return pix, wave, spec, hdr['CDELT1']

def compute_finer_wave(cur_wave, new_da):
    npix = (cur_wave.max() - cur_wave.min())
    npix /= new_da
    npix = int(npix)

    fin_wave = np.arange(0,npix)
    fin_wave *= new_da
    fin_wave += cur_wave.min()

    return fin_wave


bfn = 'bsky.eps_pang_parcsec.fits'
rfn = 'rsky.eps_pang_parcsec_onemicron.fits'

bpix, bwave, bspec, bda = read_file(bfn)
rpix, rwave, rspec, rda = read_file(bfn)

blue_max = bwave.max()
red_min = rwave.min()

fin_da = bda
if rda < bda:
    fin_da = rda

full_npix = rwave.max() - bwave.min()
full_npix /= fin_da
full_npix = int(full_npix)

fin_wave = np.arange(0, full_npix, dtype='float64')
fin_wave *= fin_da
fin_wave += bwave.min()

#blue
inblue = fin_wave < bwave.max()
bfin_spec = np.interp(fin_wave[inblue], bwave, bspec)
#red
inred = fin_wave > rwave.min()
rfin_spec = np.interp(fin_wave[inred], rwave, rspec)
# final
fin_spec = np.zeros_like(fin_wave)
# add
fin_spec[inblue] += bfin_spec
fin_spec[inred] += rfin_spec


first_hdu = astropy.io.fits.PrimaryHDU(data=fin_spec)
hdr = first_hdu.header
hdr['CDELT1'] = fin_da
hdr['CRVAL1'] = fin_wave.min()
hdu = astropy.io.fits.PrimaryHDU(data=fin_spec, header=hdr)
hdu.writeto('sky.eps_pang_parcsec_onemicron.fits', overwrite=True)




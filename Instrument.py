"""Instrument class for the snr calculator"""
import os
import astropy.io.ascii

class Instrument:
    """
    Instrument()
    """
    def __init__(self):
        self.name = ""
        self.sheight = 0.0 # Slit height
        self.swidth = 0.0 # Slit width
        self.dark = 0.0 # dark current
        self.readnoise = 0.0 # read noise
        self.bind = 1 # dispersion binning
        self.bins = 1 # spatial binning
        self.scale_para = 1.0 # pixel scale
        self.scale_perp = 1.0 # pixel scale
        self.mag_para = 1.0 #
        self.mag_perp = 1.0 #
        self.R = 1.0 # lambda/delta-lambda
        self.pixel_size= 1.
        self.mlambda= 0.0
        self.dichroic= ''
        self.grating= ''
        self.cwave= 0. # Central wavelength
        self.wvmnx= [] #
        self.dely= 0.0
        self.throughput = None
        self.BLUE_CUTOFF = 0
        self.RED_CUTOFF = 1e6

    def __repr__(self):
        return f'<Instrument {self.name} {self.grating}>'

    def read_throughput(self):
        '''
        read_throughput(self)
        '''
        grating_filename = os.path.join('data', 'throughput', self.grating + ".csv")
        self.throughput = astropy.io.ascii.read(grating_filename)

    def lris2_red(self, grating="R400"):
        '''
        lris2_red(self, grating="R400")

        Builds the red side of LRIS-2, assumse R400 grating 
        unless otherwise specified.
        '''
        self.name = "LRIS-2 Red"
        if grating not in ("R400", "R700", "R750"):
            raise ValueError(f"Grating {grating} not supported for {self.name}")

        self.grating = grating
        self.sheight = 8
        self.swidth = 0.7
        self.dark = 0.0
        self.readnoise = 3.5
        self.bind = 1
        self.bins = 1
        self.scale_para = 0.15
        self.scale_perp = 0.15
        self.mag_para = 7.4
        self.mag_perp = 7.4
        self.Ang_per_pix = 1.13 #
        self.BLUE_CUTOFF = 5500
        self.RED_CUTOFF = 9500

        self.pixel_size= 15 # 15.0 microns is 0.15 "

    def lris2_blue(self, grating="B600"):
        '''
        lris2_blue(self, grating="B600")

        Builds the blue side of LRIS-2, assumse B600 grating 
        unless otherwise specified.
        '''
        self.name = "LRIS-2 Blue"

        if grating not in ("B600", "B1200", "B1300"):
            raise ValueError(f"Grating {grating} not supported for {self.name}")

        self.grating = grating
        self.sheight = 8
        self.swidth = 0.7
        self.dark = 0.0
        self.readnoise = 3.5
        self.bind = 1
        self.bins = 1
        self.scale_para = 0.15
        self.scale_perp = 0.15
        self.mag_para = 7.4
        self.mag_perp = 7.4
        self.Ang_per_pix = 0.62 #
        self.RED_CUTOFF = 5700
        self.BLUE_CUTOFF = 3100

        self.pixel_size= 15 # 15.0 microns is 0.15 "




def deimos(tel,slit=1.0,grating="600"):
    dei = instrument()
    dei.MAG_PERP = 8.03  # Modified to give observed resolution 0.75" maps to 4.5 pixels
    dei.MAG_PARA = 8.03  # Same
    dei.PIXEL_SIZE = 15.0 # in microns

    dei.SCALE_PERP = tel.PLATE_SCALE*dei.MAG_PERP*(dei.PIXEL_SIZE/1000.) # Arcsec
    dei.SCALE_PARA = tel.PLATE_SCALE*dei.MAG_PARA*(dei.PIXEL_SIZE/1000.) # Arcsec

    dei.grating = grating
    if grating == "600":
        dei.R     = 11538.5    # 1 pixel (native) dispersion
    elif grating == "1200":
        dei.R    = 22727.3    # 1 pixel (native) dispersion
    elif grating == "900":
        dei.R    = 17307.8    # 1 pixel (native) dispersion


    ## Detector
    dei.readno = 2.6
    dei.dark = 4.  # electrons/pix/hr
    dei.bind = 1
    dei.bins = 1

    ## Wavelength range
    dei.wvmnx = [4000., 10000.]

    ## Slit
    dei.swidth = slit
    dei.sheight = 10.0  # arcsec

    # dei.throughput = 

    return dei

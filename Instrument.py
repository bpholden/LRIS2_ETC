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
        self.throughput_dir = 'data/throughput'

    def __repr__(self):
        return f'<Instrument {self.name} {self.grating}>'

    def __str__(self):
        return f'>Instrument {self.name} {self.grating}>'

    def read_xidl_throughput(self):
        '''
        read_lris_throughput(self)
        '''
        short_name = f"sens_{self.name}_{self.grating}.fits.gz"
        grating_filename = os.path.join(self.throughput_dir, short_name)
        hdus = astropy.io.fits.open(grating_filename)
        self.throughput = {}
        self.throughput['wavelength'] = hdus[2].data['WAV'][0]
        self.throughput['throughput'] = hdus[2].data['EFF'][0]

    def read_lris2_throughput(self):
        '''
        read_lris2_throughput(self)
        '''
        grating_filename = os.path.join(self.throughput_dir, self.grating + "_tot_eff.csv")
        self.throughput = astropy.io.ascii.read(grating_filename)
        self.throughput['wavelength'] *= 10 # we work in Angstroms but these tables are in nm

    def read_throughput(self):
        '''
        read_throughput(self)
        '''

        if "LRIS-2" in self.name:
            self.read_lris2_throughput()
        elif "LRIS" in self.name:
            self.read_xidl_throughput()
        elif 'DEIMOS' in self.name:
            self.read_xidl_throughput()
        else:
            raise ValueError(f"Throughput not available for {self.name}")


    def lris2_red(self, tel, grating="R400"):
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
        self.dark = 0.001
        self.readnoise = 3.5
        self.bind = 1
        self.bins = 1
        self.mag_para = 7.3
        self.mag_perp = 7.3
        self.Ang_per_pix = 1.13 #
        self.BLUE_CUTOFF = 5500.0
        self.RED_CUTOFF = 9500.0

        self.pixel_size= 15 # 15.0 microns is 0.15 "
        self.scale_perp = tel.plate_scale*self.mag_perp*(self.pixel_size/1000.) # Arcsec
        self.scale_para = tel.plate_scale*self.mag_para*(self.pixel_size/1000.)

    def lris2_blue(self, tel, grating="B600"):
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
        self.dark = 0.001
        self.readnoise = 3.5
        self.bind = 1
        self.bins = 1
        self.mag_para = 7.3
        self.mag_perp = 7.3
        self.Ang_per_pix = 0.62 #
        self.RED_CUTOFF = 5700.0
        self.BLUE_CUTOFF = 3100.0

        self.pixel_size= 15 # 15.0 microns is 0.15 "
        self.scale_perp = tel.plate_scale*self.mag_perp*(self.pixel_size/1000.) # Arcsec
        self.scale_para = tel.plate_scale*self.mag_para*(self.pixel_size/1000.)

    def lris_red(self, tel, grating="400_8500_D560"):
        '''
        lris_red(self, grating="400_8500_D560")
        '''

        self.name = "LRISr"
        if grating not in ("400_8500_D560", "600_7500_D680", "600_10000_D560"):
            raise ValueError(f"Grating {grating} not supported for {self.name}")
        self.grating = grating
        self.sheight = 8
        self.swidth = 1.2
        self.dark = 0.001
        self.readnoise = 4.5
        self.bind = 1
        self.bins = 1
        self.mag_para = 6.5
        self.mag_perp = 6.5

        self.Ang_per_pix = 1.16
        self.RED_CUTOFF = 10300.0
        self.BLUE_CUTOFF = 5600.0

        self.pixel_size= 15 

        self.scale_perp = tel.plate_scale*self.mag_perp*(self.pixel_size/1000.) # Arcsec
        self.scale_para = tel.plate_scale*self.mag_para*(self.pixel_size/1000.)


    def lris_blue(self, tel, grating="600_4000_D560"):
        '''
        lris_blue(self, grating="600_4000_D560")
        '''

        self.name = "LRISb"
        if grating not in ("600_4000_D560"):
            raise ValueError(f"Grating {grating} not supported for {self.name}")

        self.grating = grating
        self.sheight = 8
        self.swidth = 1.2
        self.dark = 0.001
        self.readnoise = 3.7
        self.bind = 1
        self.bins = 1
        self.mag_para = 6.5
        self.mag_perp = 6.5
        self.Ang_per_pix = 0.63
        self.RED_CUTOFF = 5600.0
        self.BLUE_CUTOFF = 3100.0

        self.pixel_size= 13.5

        self.scale_perp = tel.plate_scale*self.mag_perp*(self.pixel_size/1000.) # Arcsec
        self.scale_para = tel.plate_scale*self.mag_para*(self.pixel_size/1000.)


    def deimos(self, tel, grating="600"):
        '''
        deimos(self, grating="600")
        '''
        self.mag_perp = 8.03  # Modified to give observed resolution 0.75" maps to 4.5 pixels
        self.mag_para = 8.03  # Same
        self.pixel_size = 15.0 # in microns

        self.scale_perp = tel.plate_scale*self.mag_perp*(self.pixel_size/1000.) # Arcsec
        self.scale_para = tel.plate_scale*self.mag_para*(self.pixel_size/1000.)

        self.grating = grating
        if grating == "600":
            self.R     = 11538.5    # 1 pixel (native) dispersion
        elif grating == "1200":
            self.R    = 22727.3    # 1 pixel (native) dispersion
        elif grating == "900":
            self.R    = 17307.8    # 1 pixel (native) dispersion

        ## Detector
        self.readnoise = 2.6
        self.dark = 0.001  # electrons/pix/hr
        self.bind = 1
        self.bins = 1

        ## Wavelength range
        self.RED_CUTOFF = 10000
        self.BLUE_CUTOFF = 4000
        # dei.throughput
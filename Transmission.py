import numpy

class Transmission:
    def __init__(self, airmass=1.0):

        self.airmass = airmass

        self.wave = [3000,3100,3200, 3300, 3400., \
          3500,3600,3700, 3800, 3900, \
          4000,4250,4500, 4750, 5000, \
          5250,5500,5750, 6000, 6500, \
          7000,8000,9000,10000,12000 ]

        self.mag_trans =  [4.90, 1.37, 0.82, 0.57, 0.51, \
          0.42, 0.37, 0.33, 0.30, 0.27, \
          0.25, 0.21, 0.17, 0.14, 0.13, \
          0.12, 0.12, 0.12, 0.11, 0.11, \
          0.10, 0.07, 0.05, 0.04, 0.03]

        self.mag_trans = numpy.array(self.mag_trans)
        # these are magnitudes

        self.extinct = 10**(-0.4*self.mag_trans*self.airmass)

        self.extinct = numpy.array(self.extinct)
        self.wave = numpy.array(self.wave)

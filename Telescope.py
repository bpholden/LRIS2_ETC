class Telescope(object):

    def __init__(self, name=None):
        self.name = '' # KeckI, KeckII, Lick-3m
        self.area= 0.
        self.plate_scale= 0.
        if name:
            if name.lower() == 'keck1':
                self.keckone()
            elif name.lower() == 'keck2':
                self.kecktwo()
            elif name.lower() == 'lick':
                self.lick()

    def __repr__(self):
         return '<Telescope "%s" %f %f>' % (self.name, self.area, self.plate_scale)

    def lick(self):
        self.area  = 63617. # cm^2 -- 3m telescope with 10% central obscuration
        self.plate_scale = 3.86 # "/mm  for Kast
        self.name = 'Lick-3m'


    def keck(self):
        self.area  = 723674. # cm^2 -- 10m telescope with 7.9% central obscuration
        self.plate_scale = 1.379 # "/mm


    def keckone(self):
        self.keck()
        self.name = 'KeckI'


    def kecktwo(self):
        self.keck()
        self.name = 'KeckII'

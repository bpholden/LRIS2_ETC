import matplotlib.pyplot as plt

import ExpCalc
import Instrument
import Telescope
import Sky

keck_1 = Telescope.Telescope()
keck_1.keckone()
lris2_red = Instrument.Instrument()
lris2_red.lris2_red(keck_1)
lris_red = Instrument.Instrument()
lris_red.lris_red(keck_1)
lris_red.swidth = 0.7
lr_s = Sky.Sky()
lr2_s = Sky.Sky()


plt.plot(lr_s.wave, lr_s.spec)
plt.plot(lr2_s.wave, lr2_s.spec)
plt.xlim(5650,5825)
plt.ylim(0, 1)
plt.show()

lr_s.spec = Sky.rescale(lris_red, lr_s.spec)
lr2_s.spec = Sky.rescale(lris2_red, lr2_s.spec)


plt.plot(lr_s.wave, lr_s.spec)
plt.plot(lr2_s.wave, lr2_s.spec)
plt.xlim(5650,5825)
plt.ylim(0, 3)
plt.show()


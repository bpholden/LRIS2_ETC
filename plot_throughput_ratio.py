import matplotlib.pyplot as plt
import numpy as np
import ExpCalc
import Instrument
import Telescope


keck_1 = Telescope.Telescope()
keck_1.keckone()
lris2_red = Instrument.Instrument()
lris2_red.lris2_red(keck_1)
lris_red = Instrument.Instrument()
lris_red.lris_red(keck_1)


lris2_blue = Instrument.Instrument()
lris2_blue.lris2_blue(keck_1)
lris_blue = Instrument.Instrument()
lris_blue.lris_blue(keck_1)

mag = 22
filter_name = 'sdss_rprime.dat'
seeing = 1.0
airmass = 1.2
template = 'starb1'
redshift = 0.0

blue2_exp_calc = ExpCalc.ExpCalc(mag, filter_name, seeing, airmass, \
                                    redshift, template, lris2_blue, keck_1)
red2_exp_calc = ExpCalc.ExpCalc(mag, filter_name, seeing, airmass, \
                                   redshift, template, lris2_red, keck_1)
blue_exp_calc = ExpCalc.ExpCalc(mag, filter_name, seeing, airmass, \
                                    redshift, template, lris_blue, keck_1)
red_exp_calc = ExpCalc.ExpCalc(mag, filter_name, seeing, airmass, \
                                      redshift, template, lris_red, keck_1)

blue2_exp_calc.read_template()
red2_exp_calc.read_template()
blue_exp_calc.read_template()
red_exp_calc.read_template()

blue2_exp_calc.compute_throughput()
red2_exp_calc.compute_throughput()
blue_exp_calc.compute_throughput()
red_exp_calc.compute_throughput()

blue2_waves = blue2_exp_calc.waves[blue2_exp_calc.good_waves]
red2_waves = red2_exp_calc.waves[red2_exp_calc.good_waves]
blue_waves = blue_exp_calc.waves[blue_exp_calc.good_waves]
red_waves = red_exp_calc.waves[red_exp_calc.good_waves]

plt.plot(blue2_waves,
         100*blue2_exp_calc.throughput_interp, 'b-', label='LRIS-2 B600')
plt.plot(red2_waves,
          100*red2_exp_calc.throughput_interp, 'r-', label='LRIS-2 R400')
plt.plot(blue_waves,
          100*blue_exp_calc.throughput_interp, 'b--', label='LRIS 600/4000')
plt.plot(red_waves,
          100*red_exp_calc.throughput_interp, 'r--', label='LRIS 400/8500')
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Throughput')
plt.legend()
plt.ylim(0, 100)
plt.show()

red_np_interp = np.interp(red2_waves, red_exp_calc.instrument.throughput['wavelength'],\
                                       red_exp_calc.instrument.throughput['throughput'])
blue_np_interp = np.interp(blue2_waves, blue_exp_calc.instrument.throughput['wavelength'],\
                                       blue_exp_calc.instrument.throughput['throughput'])

blue_good = blue2_waves < 5550
red_good = red2_waves > 5700
blue_rat = blue2_exp_calc.throughput_interp / blue_np_interp
red_rat = red2_exp_calc.throughput_interp / red_np_interp
plt.plot(blue2_waves[blue_good], blue_rat[blue_good], 'b-', label='Ratio of B600 to 600/4000')
plt.plot(red2_waves[red_good], red_rat[red_good], 'r-', label='Ratio of R400 to 400/8500')
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Throughput Ratio')
#plt.legend()
plt.ylim(1, 3)
plt.show()

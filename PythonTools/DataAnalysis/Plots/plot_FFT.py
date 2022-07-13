
"""
Helping file to plot the FFT of some column in some file faster
"""

import sys
sys.path.append("../")
from utils import *

import numpy as np
import matplotlib.pyplot as plt

folder = "./g1/"
filename = "strain_22_cutoff_low_0.015_at_infinity.dat"

# Spacetime masses:
masses = [
    1,
]

# Files to plot:
files = [
    folder + filename,
]

# Labels for each:
labels = ["EsGB"]

column = 3

# plot_FFT_from_file(files, use_dts=True, max_frequency=0.15, column=column)
# exit()

################################################################################
use_tex()
################################################################################

datas = [np.transpose(np.loadtxt(file)) for file in files]

times = [data[0] for data in datas]

# wheres = [[t > 150 and t < 1400 for t in time] for time in times]
wheres = [[t > 0 and t < 1400 for t in time] for time in times]

strains_re_infty = [data[column] for data in datas]
# strains_im_infty = [data[column+1] for data in datas]
dts = [time[1] - time[0] for time in times]

times = [time[where] for where, time in zip(wheres, times)]
strains_re_infty = [strain[where] for where, strain in zip(wheres, strains_re_infty)]
# strains_im_infty = [strain[where] for where, strain in zip(wheres, strains_im_infty)]

FFTs = [FFT(strain) for strain in strains_re_infty]
# FFTs = [FFT(strain_re + 1j * strain_im) for strain_re, strain_im in zip(strains_re_infty, strains_im_infty)]

freqs = [FFT[0]/dt for dt, FFT in zip(dts, FFTs)]
strains_infty_FFT = [FFT[1] for FFT in FFTs]

max_frequency = 0.15

# PLOT
fig, ax1 = plt.subplots(figsize=(16*resize, 9*resize))

ax1.set_xlabel(r'$f\,M$')
ax1.set_ylabel('Re' + r'$(r \tilde{h}^+_{22}) / M$')

ax1.set_xlim(0, max_frequency)
# ax1.set_ylim(-250, 170)

for mass, label, freq, strain in zip(masses, labels, freqs, strains_infty_FFT):
    # ax1.plot(freq, np.abs(strain) / len(strain), label=label)
    # ax1.plot(freq, np.angle(strain) / len(strain), label=label)
    ax1.plot(freq * mass, np.real(strain) / mass, label=label)

ax1.legend()

# # Inset:
# ax1ins = ax1.inset_axes([0.08, 0.08, 0.9, 0.3])
# for mass, label, freq, strain in zip(masses, labels, freqs, strains_infty_FFT):
#     ax1ins.plot(freq * mass, np.real(strain) / mass, label=label)
# ax1ins.set_xlim(0.051, 0.092)
# ax1ins.set_ylim(-50, 50)
# ax1.indicate_inset_zoom(ax1ins, edgecolor="black")

plt.draw()
out_name = 'FFT_strain'
# plt.savefig(out_name + ".pdf", bbox_inches = 'tight', pad_inches=0.3)
plt.savefig(out_name + ".png", bbox_inches = 'tight', pad_inches=0.3)

plt.show()

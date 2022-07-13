
import sys
sys.path.append("../")
from utils import *

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

filename = "./run_g2_v2_h112/mismatches_vs_GR.dat"

################################################################################
use_tex()
################################################################################

data = np.loadtxt(filename, unpack=True)
masses = data[0]
mismatches = data[1] * 100

mismatch_constant = mismatches[-1]

# chop first points
masses = masses[:-1]
mismatches = mismatches[:-1]

# remember that the last element (mass=0) is the constant PSD, not LIGO PSD

# interpolate in log scale!
mismatches_interp = make_interp_spline(np.log(masses), mismatches, k=3)
masses_new = np.exp(np.linspace(np.log(masses[0]), np.log(masses[-1]), 100))
mismatches_full = mismatches_interp(np.log(masses_new))

# PLOT
fig, ax1 = plt.subplots(figsize=(16*resize, 9*resize))

ax1.set_xlabel(r'$M/M_{\odot}$')
ax1.set_ylabel('mismatch (\%)')

ax1.set_xlim(10,200)
ax1.set_ylim(9,14)
# ax1.set_yscale('log')

ax1.plot([masses_new[0], masses_new[-1]], [mismatch_constant, mismatch_constant], label="Flat PSD")
ax1.plot(masses_new, mismatches_full, label="LIGO PSD")

ax1.legend()

ax1.minorticks_on()
ax1.tick_params('both', length=20, width=2, which='major')
ax1.tick_params('both', length=10, width=1, which='minor')


plt.draw()
out_name = 'mismatches_g2_vs_GR'
plt.savefig(out_name + ".pdf", bbox_inches = 'tight')
plt.savefig(out_name + ".png", bbox_inches = 'tight')

# plt.show()

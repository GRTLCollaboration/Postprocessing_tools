
"""
Analyse GW Convergence
(converge of real part, converge of amplitude, convergence of phase, ...)
"""

import sys
sys.path.append("../Base")
sys.path.append("../")
from utils import *
import SmallDataIOReader, IntegrationMethod, DataIntegration

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import make_interp_spline # to interpolate the amplitude
from scipy.optimize import minimize # to align time shift of GW amplitude
from functools import partial

Ns = np.array([96, 112, 128])

mass = 1

files = [f"./run_g2_v{i+1}_h{Ns[i]:03d}/strain_22_cutoff_low_0.015_at_infinity.dat" for i in range(len(Ns))]
datas = [np.loadtxt(file) for file in files]

times = np.array([data[:,0] / mass for data in datas])
strains = np.array([(data[:,3] + 1j * data[:,4]) / mass for data in datas])

time_min = np.max([time[0] for time in times])
time_max = np.min([time[-1] for time in times])
print("MIN TIME =", time_min)
print("MAX TIME =", time_max)

# This supposedly returns the converge order, but it's super noisy
# def getConvergenceOrder(datas):
#     ratio = abs((datas[0] - datas[1])/(datas[1] - datas[2]))
#     initial_guess = 3.14 # just because (like this I can detect that this is not a coincidence)
#     return fsolve(func, initial_guess, args=(ratio))[0]
# convergenceOrder = [getConvergenceOrder([low, med, high]) for low, med, high in zip(amp_low_offset, amp_med_offset, amp_high_new)]
# print(convergenceOrder)

# Get convergence factor based on 'Ns' (3 grid resolutions)
func = lambda n, c : (np.power(1./Ns[0],n)-np.power(1./Ns[1],n))/\
                        (np.power(1./Ns[1],n)-np.power(1./Ns[2],n)) - c
q3_factor = func(3, 0)
q4_factor = func(4, 0)

fixed_time = 0;
strains_decomposed = [convert_complex_to_amplitude_and_phase(time, strain, fixed_time) for time, strain in zip(times, strains)]
strain_amplitudes = [strain_decomposed[0] for strain_decomposed in strains_decomposed]
strain_phases = [strain_decomposed[1] for strain_decomposed in strains_decomposed]

times_new, strain_amplitudes, strain_phases, offsets = align_waves(strain_amplitudes[-1], strain_amplitudes[:-1], times[-1], times[:-1], [-8, -2], strain_phases[-1], strain_phases[:-1])

t_merger = times_new[np.argmax(strain_amplitudes[-1])]
times_new -= t_merger

amp_low_med = [abs(low - med) for low, med in zip(strain_amplitudes[0], strain_amplitudes[1])]
amp_med_high = [abs(med - high) for med, high in zip(strain_amplitudes[1], strain_amplitudes[2])]
q3_amp_high = [ high * q3_factor for high in amp_med_high]
q4_amp_high = [ high * q4_factor for high in amp_med_high]

phase_low_med = [abs(low - med) for low, med in zip(strain_phases[0], strain_phases[1])]
phase_med_high = [abs(med - high) for med, high in zip(strain_phases[1], strain_phases[2])]
q3_phase_high = [ high * q3_factor for high in phase_med_high]
q4_phase_high = [ high * q4_factor for high in phase_med_high]

# plt.plot(times_new, phase_low_med, label = r"$Low-Med$")
# plt.plot(times_new, phase_med_high, label = r"$Med-High$")
# plt.plot(times_new, q3_phase_high, "--", label = r"$Q3(Med-High)$")
# plt.plot(times_new, q4_phase_high, "--", label = r"$Q4(Med-High)$")
# plt.yscale("log")
# plt.legend()
# plt.show()
# exit()

################################################################################
use_tex(32, 24)
################################################################################


fig, axs = plt.subplots(3, figsize=(9*resize, 16*resize))
[ax1, ax2, ax3] = axs
ax1_twin = ax1.twinx()

fig.subplots_adjust(wspace=0, hspace=0)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

for ax in [ax1, ax1_twin, ax2, ax3]:
    ax.minorticks_on()
    ax.tick_params('both', length=20, width=2, which='major', direction='in')
    ax.tick_params('both', length=10, width=1, which='minor', direction='in')

    ax.axvline(0, linestyle='--', color='black') # vertical lines
    ax.set_xlim(-1160, 60)

# PLOT 1 - amplitude and phase
leg = []
leg += ax1.plot(times_new, strain_amplitudes[-1], color=plt.cm.Set1(1))
leg += ax1_twin.plot(times_new, -strain_phases[-1], color=plt.cm.Set1(2))

leg = fig.legend(leg, [r"$r h^A_{22}/M$", r"$h^\phi_{22}$"]
    , loc="upper left"
    , bbox_to_anchor=(0.135, 0.83)
)
plt.gca().add_artist(leg)

ax1.set_ylabel(r"$r h_{22}^A / M$")
ax1_twin.set_ylabel(r"$h_{22}^\phi$", rotation=270, labelpad=70)

# PLOT 2 - amplitude convergence
pipe = r"$|$"
leg = []
leg += ax2.plot(times_new, amp_low_med, label = pipe + "Low-Med" + pipe)
leg += ax2.plot(times_new, amp_med_high, label = pipe + "Med-High" + pipe)
leg += ax2.plot(times_new, q3_amp_high, "--", label = r"$Q_3(|$" + "Med-High" + r"$|)$")
leg += ax2.plot(times_new, q4_amp_high, "--", label = r"$Q_4(|$" + "Med-High" + r"$|)$")

leg = fig.legend(leg, [l.get_label() for l in leg]
    , loc="upper left"
    , bbox_to_anchor=(0.47, 0.39)
)
leg.get_frame().set_alpha(None) # set transparent legend
leg.get_frame().set_facecolor((1, 1, 1, 0.1))
plt.gca().add_artist(leg)

ax2.set_yscale('log')
ax2.set_ylim(1e-7, 8e-2)
ax2.set_ylabel(r'$\Delta (r h_{22}^A / M)$')

# PLOT 3 - phase convergence
ax3.plot(times_new, phase_low_med, label = pipe + "Low-Med" + pipe)
ax3.plot(times_new, phase_med_high, label = pipe + "Med-High" + pipe)
ax3.plot(times_new, q3_phase_high, "--", label = r"$Q_3(|$" + "Med-High" + r"$|)$")
ax3.plot(times_new, q4_phase_high, "--", label = r"$Q_4(|$" + "Med-High" + r"$|)$")

ax3.set_yscale('log')
ax3.set_ylim(1e-5, 1e1)
ax3.set_xlabel(r'$\Delta t/M$')
ax3.set_ylabel(r'$\Delta h_{22}^\phi$')
# ax3.legend()

plt.draw()
out_name = 'GW_strain_convergence'
plt.savefig(out_name + ".pdf", bbox_inches = 'tight')
plt.savefig(out_name + ".png", bbox_inches = 'tight')

# plt.show()

# ##################################################################
# # PLOT Weyl4 Real Part CONVERGENCE

# ax2.set_ylim(5e-11, 1e-2)
# ax2.set_xlim(-200, 1450)

# ax2.legend(loc="lower left", bbox_to_anchor=(0.0,0.8))

# ax2ins = ax2.inset_axes([0.47, 0.1, 0.5, 0.4])
# ax2ins.plot(times_minus_r, weyl4_re_low_med, label = r"$Low-Med$")
# ax2ins.plot(times_minus_r, weyl4_re_med_high, label = r"$Med-High$")
# ax2ins.plot(times_minus_r, q3_weyl4_re_high, "--", label = r"$Q3(Med-High)$")
# ax2ins.plot(times_minus_r, q4_weyl4_re_high, "--", label = r"$Q4(Med-High)$")
# ax2ins.set_xlim(1125/initial_mass_g2_0_005_v3, 1225/initial_mass_g2_0_005_v3)
# ax2ins.set_ylim(1e-5, 7e-3)
# ax2ins.set_yscale("log")
# ax2ins.minorticks_on()
# ax2ins.tick_params('both', length=20, width=2, which='major')
# ax2ins.tick_params('both', length=10, width=1, which='minor')
# ax2.indicate_inset_zoom(ax2ins, edgecolor="black")



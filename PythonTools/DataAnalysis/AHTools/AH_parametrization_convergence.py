import numpy as np
import matplotlib.pyplot as plt
import math # sqrt, atan

import sys
sys.path.append("../")
from utils import *

#################################################################
# inputs

folders = ["./run_g2_v1_h096/", "./run_g2_v2_h112/", "./run_g2_v3_h128/"]
masses = [1, 1, 1]
center = 512 # L/2

#################################################################
# process data

N = len(folders)

files_AH1 = [folder + "stats_AH1.dat" for folder in folders]
files_AH2 = [folder + "stats_AH2.dat" for folder in folders]

data_AH1 = [np.loadtxt(file) for file in files_AH1]
data_AH2 = [np.loadtxt(file) for file in files_AH2]

# get headers
headers = [0]*N
for d in range(N):
    file = open(folders[d] + "stats_AH1.dat", 'r')
    headers[d] = file.readline()
    file.close()

# minimum/maximum time that appears at least in all files
# assumes AH2 has same times
min_time = max([min(data_AH1[d][:,0]) / masses[d] for d in range(N)])
max_time = min([max(data_AH1[d][:,0]) / masses[d] for d in range(N)])

# create interpolations of x and y positions
from scipy.interpolate import make_interp_spline
spline_x_AH1 = [0]*N
spline_y_AH1 = [0]*N
spline_x_AH2 = [0]*N
spline_y_AH2 = [0]*N

for d in range(N):
    hasNewSpin = ('spin-z-alt' in headers[d])
    x_idx = 13 if hasNewSpin else 8
    spline_x_AH1[d] = make_interp_spline(data_AH1[d][:,0] / masses[d], data_AH1[d][:,x_idx] - center, k=3)   # type: BSpline
    spline_y_AH1[d] = make_interp_spline(data_AH1[d][:,0] / masses[d], data_AH1[d][:,x_idx+1] - center, k=3)  # type: BSpline
    spline_x_AH2[d] = make_interp_spline(data_AH2[d][:,0] / masses[d], data_AH2[d][:,x_idx] - center, k=3)   # type: BSpline
    spline_y_AH2[d] = make_interp_spline(data_AH2[d][:,0] / masses[d], data_AH2[d][:,x_idx+1] - center, k=3)  # type: BSpline

# calculate D's and theta's
Npoints = 300
ts = np.linspace(min_time, max_time, Npoints)


Ds = [0]*N
thetas = [0]*N
for d in range(N):
    x1 = spline_x_AH1[d](ts)
    y1 = spline_y_AH1[d](ts)
    x2 = spline_x_AH2[d](ts)
    y2 = spline_y_AH2[d](ts)

    Ds[d] = np.sqrt((x2-x1)**2 + (y2-y1)**2) / masses[d]
    thetas[d] = np.arctan2(y2-y1, x2-x1)

    # remove discontinuities
    theta_0 = thetas[d][0]
    for p in range(0,len(thetas[d])):
        thetas[d][p] = theta_0 - thetas[d][p] # invert rotation
        if p>0 and abs(thetas[d][p] - thetas[d][p-1]) > math.pi:
            thetas[d][p] = np.nan

######################
Ns = np.array([96, 112, 128])
def func(n, c):
    return (np.power(1./Ns[0],n)-np.power(1./Ns[1],n))/\
           (np.power(1./Ns[1],n)-np.power(1./Ns[2],n)) - c
q3_factor = func(3, 0)
q4_factor = func(4, 0)
#####################

D_diffLowMed = np.abs(Ds[0] - Ds[1])
D_diffMedHigh = np.abs(Ds[1] - Ds[2])
q3_D_high = D_diffMedHigh * q3_factor
q4_D_high = D_diffMedHigh * q4_factor

theta_diffLowMed = np.abs(thetas[0] - thetas[1])
theta_diffMedHigh = np.abs(thetas[1] - thetas[2])
q3_theta_high = theta_diffMedHigh * q3_factor
q4_theta_high = theta_diffMedHigh * q4_factor

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

# PLOT 1 - D and theta
leg = []
leg += ax1.plot(ts, Ds[0+2], color=plt.cm.Set1(1))
leg += ax1_twin.plot(ts, thetas[0+2], color=plt.cm.Set1(2))

leg = fig.legend(leg, [r"$D/M$", r"$\theta$"]
    , loc="upper left"
    , bbox_to_anchor=(0.135, 0.83)
)
plt.gca().add_artist(leg)
leg.get_frame().set_alpha(None) # set transparent legend
leg.get_frame().set_facecolor((1, 1, 1, 0.1))

ax1.set_ylabel(r"$D/M$")
ax1_twin.set_ylabel(r"$\theta$", rotation=270, labelpad=70)
ax1.set_ylim([0, 40])

# PLOT 2 - amplitude convergence
pipe = r"$|$"
leg = []
leg += ax2.plot(ts, D_diffLowMed, label = pipe + "Low-Med" + pipe)
leg += ax2.plot(ts, D_diffMedHigh, label = pipe + "Med-High" + pipe)
leg += ax2.plot(ts, q3_D_high, "--", label = r"$Q_3(|$" + "Med-High" + r"$|)$")
leg += ax2.plot(ts, q4_D_high, "--", label = r"$Q_4(|$" + "Med-High" + r"$|)$")

leg = fig.legend(leg, [l.get_label() for l in leg]
    , loc="upper left"
    , bbox_to_anchor=(0.135, 0.3)
)
leg.get_frame().set_alpha(None) # set transparent legend
leg.get_frame().set_facecolor((1, 1, 1, 0.1))
plt.gca().add_artist(leg)

# ax2.set_yscale('log')
# ax2.set_ylim(1e-7, 8e-2)
ax2.set_ylabel(r'$\Delta D/M$')

# PLOT 3 - phase convergence
ax3.plot(ts, theta_diffLowMed, label = pipe + "Low-Med" + pipe)
ax3.plot(ts, theta_diffMedHigh, label = pipe + "Med-High" + pipe)
ax3.plot(ts, q3_theta_high, "--", label = r"$Q_3(|$" + "Med-High" + r"$|)$")
ax3.plot(ts, q4_theta_high, "--", label = r"$Q_4(|$" + "Med-High" + r"$|)$")

# ax3.set_yscale('log')
# ax3.set_ylim(1e-5, 1e1)
ax3.set_xlabel(r'$t/M$')
ax3.set_ylabel(r'$\Delta\theta$')
# ax3.legend()

plt.draw()
out_name = 'parametrization_convergence'
plt.savefig(out_name + ".pdf", bbox_inches = 'tight')
plt.savefig(out_name + ".png", bbox_inches = 'tight')

# plt.show()

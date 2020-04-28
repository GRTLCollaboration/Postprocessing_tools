# parallel_psi.py
# James Widdicombe
# Last Updated 16/12/2019
# Calculate phi evolution

# Load the modules
import yt
import numpy as np
import time
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Timings
start_time = time.time()

# Enable Parallelism
yt.enable_parallelism()

# Plot parameters
line = 2.5  # Line width
alp = 0.6  # alpha

rcParams.update({"figure.autolayout": True})
rcParams["axes.formatter.limits"] = [-3, 3]
rcParams["font.size"] = 12

# CHANGE ME
# Loading dataset (Load from one folder up)
data_location = "../../outMatterSF_00*"  # Data file location

# Loading dataset
ts = yt.load(data_location)

# Program Parameters
center = ts[0].domain_right_edge / 2.0
adjusted_left = 60
adjusted_right = int(center[0]) * 2 - adjusted_left

# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

for sto, i in ts.piter(storage=storage):

    # All Data
    ad = i.r[
        adjusted_left:adjusted_right,
        adjusted_left:adjusted_right,
        adjusted_left:adjusted_right,
    ]

    maxphi = ad.max("phi")
    minphi = ad.min("phi")

    array = [i.current_time, maxphi, minphi]

    sto.result = array
    sto.result_id = str(i)
if yt.is_root():
    timedata = []
    maxphidata = []
    minphidata = []

    for L in sorted(storage.items()):
        timedata.append(L[1][0])
        maxphidata.append(L[1][1])
        minphidata.append(L[1][2])
    # L2H
    plt.figure(1)
    plt.plot(timedata, maxphidata)
    plt.plot(timedata, minphidata)
    plt.ylabel("$\\phi$ $[M_{pl}]$")
    plt.xlabel("Time $[1/m]$")
    plt.grid()
    plt.savefig("phi.png", bbox_inches="tight")
    plt.close()

    np.savetxt("time.out", timedata)
    np.savetxt("maxphi.out", maxphidata)
    np.savetxt("minphi.out", minphidata)

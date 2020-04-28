# parallel_L.py
# James Widdicombe
# Last Updated 16/12/2019
# Calculate L2M and L2H

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

# Matplotlib Settings
rcParams.update({"figure.autolayout": True})
rcParams["axes.formatter.limits"] = [-3, 3]

# CHANGE ME
# Loading dataset (Load from one folder up)
data_location = "../../plt*"  # Data file location

# Loading dataset
ts = yt.load(data_location)

# Program Parameters
center = ts[0].domain_right_edge / 2.0
adjusted_left = 8
adjusted_right = int(center[0]) * 2 - adjusted_left

# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

# Define H2
def _H2(field, data):
    return data["Ham"] * data["Ham"]


# Define M2
def _M2(field, data):
    return (
        data["Mom1"] * data["Mom1"]
        + data["Mom2"] * data["Mom2"]
        + data["Mom3"] * data["Mom3"]
    )


for sto, i in ts.piter(storage=storage):

    # All Data
    ad = i.r[
        adjusted_left:adjusted_right,
        adjusted_left:adjusted_right,
        adjusted_left:adjusted_right,
    ]

    # Add the M2 and L2 Fields
    i.add_field("H2", _H2, units="")
    i.add_field("M2", _M2, units="")

    # L2H
    meanH2 = ad.mean("H2", weight="cell_volume")
    L2H = np.sqrt(meanH2)

    # L2M
    meanM2 = ad.mean("M2", weight="cell_volume")
    L2M = np.sqrt(meanM2)

    array = [i.current_time, L2H, L2M]

    sto.result = array
    sto.result_id = str(i)
if yt.is_root():
    timedata = []
    L2Hdata = []
    L2Mdata = []

    for L in sorted(storage.items()):
        timedata.append(L[1][0])
        L2Hdata.append(L[1][1])
        L2Mdata.append(L[1][2])
    # L2H
    plt.figure(1)
    plt.plot(timedata, L2Hdata)
    plt.ylabel("$\\mathcal{H}$")
    plt.xlabel("Time $[1/m]$")
    plt.grid()
    plt.savefig("L2H.png")
    plt.close()

    # L2M
    plt.figure(2)
    plt.plot(timedata, L2Mdata)
    plt.ylabel("$\\mathcal{M}$")
    plt.xlabel("Time $[1/m]$")
    plt.grid()
    plt.savefig("L2M.png")
    plt.close()

    np.savetxt("time.out", timedata)
    np.savetxt("L2H.out", L2Hdata)
    np.savetxt("L2M.out", L2Mdata)

# parallel_energy.py
# James Widdicombe
# Last Updated 16/12/2019

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
data_location = "../../outMatterSF_0*"  # Data file location

# Loading dataset
ts = yt.load(data_location)

# Program Parameters
center = ts[0].domain_right_edge / 2.0
adjusted_right = int(center[0]) * 2 - 8


def _V(field, data):
    return 0.5 * data["phi"] ** 2


def _Kinetic(field, data):
    return 0.5 * data["Pi"] ** 2


def _Gradient(field, data):
    return 0.5 * data["phi_gradient_magnitude"] ** 2


# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

for sto, i in ts.piter(storage=storage):

    i.add_gradient_fields(("chombo", u"phi"))
    i.add_field("V", _V, units="")
    i.add_field("Kinetic", _Kinetic, units="")
    i.add_field("Gradient", _Gradient, units="cm**(-2)")

    # All Data
    ad = i.r[
        int(center[0]) - 50 : int(center[0]),
        int(center[0]) - 20 : int(center[0]) + 20,
        int(center[0]) - 20 : int(center[0]) + 20,
    ]

    meanrho = ad.mean("rho", weight="cell_volume")
    Pot = ad.mean("V", weight="cell_volume")
    Kin = float(ad.mean("Kinetic", weight="cell_volume"))
    Grad = float(ad.mean("Gradient", weight="cell_volume"))
    S1 = float(ad.mean("s1", weight="cell_volume"))

    array = [i.current_time, meanrho, Pot, Kin, Grad, S1]

    sto.result = array
    sto.result_id = str(i)
if yt.is_root():
    timedata = []
    rhodata = []
    potdata = []
    kindata = []
    graddata = []
    s1data = []
    vdata = []

    for L in sorted(storage.items()):
        timedata.append(L[1][0])
        rhodata.append(L[1][1])
        potdata.append(L[1][2])
        kindata.append(L[1][3])
        graddata.append(L[1][4])
        s1data.append(L[1][5])
    for i in range(0, len(s1data)):
        ptemp = s1data[i] * (50 * 40 * 40)
        mstar = 0.334311
        v = ptemp / ((mstar * mstar) + (ptemp * ptemp))
        vdata.append(v)
    np.savetxt("time.out", timedata)
    np.savetxt("rho.out", rhodata)
    np.savetxt("pot.out", potdata)
    np.savetxt("kin.out", kindata)
    np.savetxt("grad.out", graddata)
    np.savetxt("s1.out", s1data)
    np.savetxt("v.out", vdata)

    plt.figure(1)
    plt.plot(timedata, rhodata)
    plt.plot(timedata, potdata)
    plt.plot(timedata, kindata)
    plt.plot(timedata, graddata)
    plt.ylabel("Energy")
    plt.xlabel("Time $[1/m]$")
    plt.xlim(500, 1000)
    plt.grid()
    plt.savefig("Energy.png")
    plt.close()

    plt.figure(2)
    plt.plot(timedata, s1data)
    plt.grid()
    plt.savefig("s1.png")
    plt.close()

    plt.figure(3)
    plt.plot(timedata, vdata)
    plt.grid()
    plt.savefig("v.png")
    plt.close()

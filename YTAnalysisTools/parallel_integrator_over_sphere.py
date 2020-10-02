# parallel_spherical_integrator.py
# James Widdicombe
# Last Updated 28/09/2020
# Integrate a field over a sphere

# Import the modules
from matplotlib import rcParams
import matplotlib.pyplot as plt
import yt
import numpy as np
from numpy import pi
import time
import matplotlib

matplotlib.use("Agg")

start_time = time.time()

# Matplotlib Settings
rcParams.update({"figure.autolayout": True})
rcParams["axes.formatter.limits"] = [-3, 3]

# Enable YT Parallelism
yt.enable_parallelism()

# Script Parameters
extraction_radius = 60  # radius of extraction
# Name of variable
RealPart = "chi"
# Imaginary part of the variable you want to integrate over
ImagPart = None
data_location = "../../ScalarField_*.3d.hdf5"  # Data file location

# Loading dataset
ts = yt.load(data_location)

# Define Deltat for power output
time0 = ts[0].current_time
time1 = ts[1].current_time
DeltaT = float(time1 - time0)

# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

# Program Parameters
center = ts[0].domain_right_edge / 2.0

# Definitions for Quadrature scheme
N = 131

coord = np.loadtxt("PointDistFiles/lebedev/lebedev_%03d.txt" % N)
theta = coord[:, 1] * pi / 180
phi = coord[:, 0] * pi / 180 + pi
w = coord[:, 2]

phi_length = len(phi)

# Iterate through dataset
for sto, i in ts.piter(storage=storage):
    # Timings
    L_start = time.time()

    # Init
    Integral = 0 + 1j * 0

    Weyl4data = []

    # k is a counter
    for k in range(phi_length):

        phi_var = phi[k]
        theta_var = theta[k]
        x1 = extraction_radius * \
            np.cos(phi_var) * np.sin(theta_var) + float(center[0])
        y1 = extraction_radius * \
            np.sin(phi_var) * np.sin(theta_var) + float(center[1])
        z1 = extraction_radius * np.cos(theta_var) + float(center[2])
        c = [x1, y1, z1]
        ptn = i.point(c)
        Real = float(ptn[RealPart][0])
        if ImagPart is not None:
            Imag = float(ptn[ImagPart][0])
        else:
            Imag = 0
        Integrand = Real + 1j * Imag

        Integral += (
            4 * pi * w[k] * Integrand * extraction_radius ** 2
        )

        # positive m
    array = [
        i.current_time,
        Integral,
        time.time() - L_start,
    ]
    sto.result = array
    sto.result_id = str(i)
if yt.is_root():

    # Initialize Output arrays
    # l = 2

    timedata = []

    # Diagnostics
    loop_times = []

    Integrated_values = []

    # Swap from storage into arrays
    for L in sorted(storage.items()):
        timedata.append(L[1][0])
        Integrated_values.append(L[1][1])
        loop_times.append(L[1][2])

    All_data = []
    All_data.append(Integrated_values)

    # Output Data
    np.savetxt("time.out", timedata)
    np.savetxt("loop_times.out", loop_times)
    np.savetxt("speed_up.out", [sum(loop_times) / (time.time() - start_time)])
    np.savetxt("Integrated_values.out", Integrated_values)

    # Integrated Plotting
    labels = [
        "integrated_evolution ",
    ]

    # Calculate Retarded Time
    time_retarded = []
    for i in timedata:
        time_retarded.append(float(i) - extraction_radius)
    np.savetxt("timeretarded.out", time_retarded)

    # Individual Modes
    for i in range(0, len(labels)):
        plt.figure(i)
        plt.plot(time_retarded, np.real(All_data[i]))
        plt.xlabel(r"$t_{ret}~[1/m]$")
        plt.ylabel(r"$r\Psi_4$")
        plt.grid()
        plt.savefig(labels[i] + ".png", bbox_inches="tight")
        plt.close()

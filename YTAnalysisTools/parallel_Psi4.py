# parallel_Psi4.py
# James Widdicombe
# Last Updated 16/12/2019
# Decomposition of Psi4 into spin weighted spherical harmonics
# l = 2,3,4

# Import the modules
import yt
import numpy as np
from numpy import pi
import time
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

start_time = time.time()

# Matplotlib Settings
rcParams.update({"figure.autolayout": True})
rcParams["axes.formatter.limits"] = [-3, 3]

# Enable YT Parallelism
yt.enable_parallelism()

# Script Parameters
extraction_radius = 60  # radius of extraction
data_location = "../outMatterSF_*.3d.hdf5"  # Data file location

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

# Spinweighted spherical harmonics s = -2
# =========================================
# L = 2
# =========================================
# m = 0 zero
sY_l2_m0 = np.sqrt(15.0 / (32.0 * np.pi)) * np.sin(theta) ** 2
# positive m values : 1,2
sY_l2_m1 = (
    np.exp(1j * phi)
    * np.sqrt(5.0 / (np.pi))
    * np.cos(theta / 2) ** 3
    * np.sin(theta / 2)
)
sY_l2_m2 = (
    1.0 / 2.0 * np.exp(2 * 1j * phi) * np.sqrt(5.0 / (np.pi)) * np.cos(theta / 2) ** 4
)
# negative m values :-1,-2
sY_l2_m1n = (
    1.0
    / 2.0
    * np.exp(-1j * phi)
    * np.sqrt(5.0 / (np.pi))
    * np.sin(theta / 2) ** 2
    * np.sin(theta)
)
sY_l2_m2n = (
    1.0 / 2.0 * np.exp(-2 * 1j * phi) * np.sqrt(5.0 / (np.pi)) * np.sin(theta / 2) ** 4
)
# =========================================
# L = 3
# =========================================
# m = 0 zero
sY_l3_m0 = (
    1.0 / 4.0 * np.sqrt(105.0 / (2.0 * np.pi)) * np.cos(theta) * np.sin(theta) ** 2
)
# positive m values : 1, 2 ...
sY_l3_m1 = (
    1.0
    / 2.0
    * np.exp(1j * phi)
    * np.sqrt(35.0 / (2.0 * np.pi))
    * np.cos(theta / 2) ** 3
    * (-1 + 3 * np.cos(theta))
    * np.sin(theta / 2)
)
sY_l3_m2 = (
    1.0
    / 2.0
    * np.exp(2 * 1j * phi)
    * np.sqrt(7.0 / (np.pi))
    * np.cos(theta / 2) ** 4
    * (-2 + 3 * np.cos(theta))
)
sY_l3_m3 = (
    -np.exp(3 * 1j * phi)
    * np.sqrt(21.0 / (2 * np.pi))
    * np.cos(theta / 2) ** 5
    * np.sin(theta / 2)
)
# negative m values : -1, -2 ...
sY_l3_m1n = (
    1.0
    / 2.0
    * np.exp(-1j * phi)
    * np.sqrt(35.0 / (2 * np.pi))
    * np.cos(theta / 2)
    * (1 + 3 * np.cos(theta))
    * np.sin(theta / 2) ** 3
)
sY_l3_m2n = (
    1.0
    / 2.0
    * np.exp(-2 * 1j * phi)
    * np.sqrt(7.0 / (np.pi))
    * (2 + 3 * np.cos(theta))
    * np.sin(theta / 2) ** 4
)
sY_l3_m3n = (
    1.0
    / 2.0
    * np.exp(-3 * 1j * phi)
    * np.sqrt(21.0 / (2 * np.pi))
    * np.sin(theta / 2) ** 4
    * np.sin(theta)
)
# =========================================
# L = 4
# =========================================
sY_l4_m0 = (
    3.0
    / 16.0
    * np.sqrt(5.0 / (2.0 * np.pi))
    * (5.0 + 7.0 * np.cos(2.0 * theta))
    * np.sin(theta) ** 2
)
# positive m = 1,2,3,4
sY_l4_m1 = (
    3
    / (2 * np.sqrt(2 * np.pi))
    * np.exp(1j * phi)
    * np.cos(theta / 2) ** 3
    * (6 - 7 * np.cos(theta) + 7 * np.cos(2 * theta))
    * np.sin(theta / 2)
)
sY_l4_m2 = (
    3
    / (4 * np.sqrt(np.pi))
    * np.exp(2.0 * 1j * phi)
    * (9 - 14 * np.cos(theta) + 7 * np.cos(2 * theta))
    * np.cos(theta / 2) ** 4
)
sY_l4_m3 = (
    -3.0
    * np.sqrt(7.0 / (2.0 * np.pi))
    * np.cos(theta / 2) ** 5
    * np.exp(3.0 * 1j * phi)
    * (-1 + 2 * np.cos(theta))
    * np.sin(theta / 2)
)
sY_l4_m4 = (
    3.0
    * np.exp(4 * 1j * phi)
    * np.sqrt(7 / np.pi)
    * np.cos(theta / 2) ** 6
    * np.sin(theta / 2) ** 2
)
# negative m = -1,-2,-3,-4
sY_l4_m1n = (
    3.0
    / (2.0 * np.sqrt(2 * np.pi))
    * np.exp(-1j * phi)
    * np.cos(theta / 2)
    * (6 + 7 * np.cos(theta) + 7 * np.cos(2 * theta))
    * np.sin(theta / 2) ** 3
)
sY_l4_m2n = (
    3.0
    / (4.0 * np.sqrt(np.pi))
    * np.exp(-2 * 1j * phi)
    * (9 + 14 * np.cos(theta) + 7 * np.cos(2 * theta))
    * np.sin(theta / 2) ** 4
)
sY_l4_m3n = (
    3.0
    * np.sqrt(7.0 / (2.0 * np.pi))
    * np.exp(-3 * 1j * phi)
    * (1 + 2 * np.cos(theta))
    * np.sin(theta / 2) ** 5
    * np.cos(theta / 2)
)
sY_l4_m4n = (
    3.0
    / 4.0
    * np.exp(-4 * 1j * phi)
    * np.sqrt(7 / np.pi)
    * np.sin(theta / 2) ** 4
    * np.sin(theta) ** 2
)

# Iterate through dataset
for sto, i in ts.piter(storage=storage):
    # Timings
    L_start = time.time()

    # Initalising
    Sr = 0 + 1j * 0
    # l = 2
    Weyl4_l2_m0 = 0 + 1j * 0
    # positive m
    Weyl4_l2_m1 = 0 + 1j * 0
    Weyl4_l2_m2 = 0 + 1j * 0
    # negative m
    Weyl4_l2_m1n = 0 + 1j * 0
    Weyl4_l2_m2n = 0 + 1j * 0

    # l = 3
    Weyl4_l3_m0 = 0 + 1j * 0
    # positive m
    Weyl4_l3_m1 = 0 + 1j * 0
    Weyl4_l3_m2 = 0 + 1j * 0
    Weyl4_l3_m3 = 0 + 1j * 0

    # negative m
    Weyl4_l3_m1n = 0 + 1j * 0
    Weyl4_l3_m2n = 0 + 1j * 0
    Weyl4_l3_m3n = 0 + 1j * 0

    # l = 4
    Weyl4_l4_m0 = 0 + 1j * 0
    # positive m
    Weyl4_l4_m1 = 0 + 1j * 0
    Weyl4_l4_m2 = 0 + 1j * 0
    Weyl4_l4_m3 = 0 + 1j * 0
    Weyl4_l4_m4 = 0 + 1j * 0
    # negative m
    Weyl4_l4_m1n = 0 + 1j * 0
    Weyl4_l4_m2n = 0 + 1j * 0
    Weyl4_l4_m3n = 0 + 1j * 0
    Weyl4_l4_m4n = 0 + 1j * 0

    Weyl4data = []

    # k is a counter
    for k in range(phi_length):

        phi_var = phi[k]
        theta_var = theta[k]
        x1 = extraction_radius * np.cos(phi_var) * np.sin(theta_var) + float(center[0])
        y1 = extraction_radius * np.sin(phi_var) * np.sin(theta_var) + float(center[1])
        z1 = extraction_radius * np.cos(theta_var) + float(center[2])
        c = [x1, y1, z1]
        ptn = i.point(c)
        ReWeyl = float(ptn["ReWeyl4"][0])
        ImWeyl = float(ptn["ImWeyl4"][0])
        Weyl4 = ReWeyl + 1j * ImWeyl

        Weyl4data.append(Weyl4)

        Weyl4_l2_m0 += (
            4 * pi * w[k] * np.conjugate(sY_l2_m0[k]) * Weyl4 * extraction_radius
        )
        # positive m
        Weyl4_l2_m1 += (
            4 * pi * w[k] * np.conjugate(sY_l2_m1[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l2_m2 += (
            4 * pi * w[k] * np.conjugate(sY_l2_m2[k]) * Weyl4 * extraction_radius
        )
        # negative m
        Weyl4_l2_m1n += (
            4 * pi * w[k] * np.conjugate(sY_l2_m1n[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l2_m2n += (
            4 * pi * w[k] * np.conjugate(sY_l2_m2n[k]) * Weyl4 * extraction_radius
        )
        # l = 3
        Weyl4_l3_m0 += (
            4 * pi * w[k] * np.conjugate(sY_l3_m0[k]) * Weyl4 * extraction_radius
        )
        # positive m
        Weyl4_l3_m1 += (
            4 * pi * w[k] * np.conjugate(sY_l3_m1[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l3_m2 += (
            4 * pi * w[k] * np.conjugate(sY_l3_m2[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l3_m3 += (
            4 * pi * w[k] * np.conjugate(sY_l3_m3[k]) * Weyl4 * extraction_radius
        )
        # negative m
        Weyl4_l3_m1n += (
            4 * pi * w[k] * np.conjugate(sY_l3_m1n[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l3_m2n += (
            4 * pi * w[k] * np.conjugate(sY_l3_m2n[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l3_m3n += (
            4 * pi * w[k] * np.conjugate(sY_l3_m3n[k]) * Weyl4 * extraction_radius
        )
        # l = 4
        Weyl4_l4_m0 += (
            4 * pi * w[k] * np.conjugate(sY_l4_m0[k]) * Weyl4 * extraction_radius
        )
        # positive m
        Weyl4_l4_m1 += (
            4 * pi * w[k] * np.conjugate(sY_l4_m1[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l4_m2 += (
            4 * pi * w[k] * np.conjugate(sY_l4_m2[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l4_m3 += (
            4 * pi * w[k] * np.conjugate(sY_l4_m3[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l4_m4 += (
            4 * pi * w[k] * np.conjugate(sY_l4_m4[k]) * Weyl4 * extraction_radius
        )
        # negative m
        Weyl4_l4_m1n += (
            4 * pi * w[k] * np.conjugate(sY_l4_m1n[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l4_m2n += (
            4 * pi * w[k] * np.conjugate(sY_l4_m2n[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l4_m3n += (
            4 * pi * w[k] * np.conjugate(sY_l4_m3n[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l4_m4n += (
            4 * pi * w[k] * np.conjugate(sY_l4_m4n[k]) * Weyl4 * extraction_radius
        )
    array = [
        i.current_time,
        Weyl4_l2_m0,
        Weyl4_l2_m1,
        Weyl4_l2_m2,
        Weyl4_l2_m1n,
        Weyl4_l2_m2n,
        Weyl4_l3_m0,
        Weyl4_l3_m1,
        Weyl4_l3_m2,
        Weyl4_l3_m3,
        Weyl4_l3_m1n,
        Weyl4_l3_m2n,
        Weyl4_l3_m3n,
        Weyl4_l4_m0,
        Weyl4_l4_m1,
        Weyl4_l4_m2,
        Weyl4_l4_m3,
        Weyl4_l4_m4,
        Weyl4_l4_m1n,
        Weyl4_l4_m2n,
        Weyl4_l4_m3n,
        Weyl4_l4_m4n,
        time.time() - L_start,
        Weyl4data,
    ]
    sto.result = array
    sto.result_id = str(i)
if yt.is_root():

    # Initialize Output arrays
    # l = 2
    Weyl4_l2_m0_data = []
    # positive m
    Weyl4_l2_m1_data = []
    Weyl4_l2_m2_data = []
    # negative m
    Weyl4_l2_m1n_data = []
    Weyl4_l2_m2n_data = []
    # l = 3
    Weyl4_l3_m0_data = []
    # positive m
    Weyl4_l3_m1_data = []
    Weyl4_l3_m2_data = []
    Weyl4_l3_m3_data = []
    # negative m
    Weyl4_l3_m1n_data = []
    Weyl4_l3_m2n_data = []
    Weyl4_l3_m3n_data = []
    # l = 4
    Weyl4_l4_m0_data = []
    # positive m
    Weyl4_l4_m1_data = []
    Weyl4_l4_m2_data = []
    Weyl4_l4_m3_data = []
    Weyl4_l4_m4_data = []
    # negative m
    Weyl4_l4_m1n_data = []
    Weyl4_l4_m2n_data = []
    Weyl4_l4_m3n_data = []
    Weyl4_l4_m4n_data = []

    timedata = []

    # Diagnostics
    loop_times = []

    Weyl4_all_time = []

    # Swap from storage into arrays
    for L in sorted(storage.items()):
        timedata.append(L[1][0])

        Weyl4_l2_m0_data.append(L[1][1])
        # positive m
        Weyl4_l2_m1_data.append(L[1][2])
        Weyl4_l2_m2_data.append(L[1][3])
        # negative m
        Weyl4_l2_m1n_data.append(L[1][4])
        Weyl4_l2_m2n_data.append(L[1][5])
        # l = 3
        Weyl4_l3_m0_data.append(L[1][6])
        # positive m
        Weyl4_l3_m1_data.append(L[1][7])
        Weyl4_l3_m2_data.append(L[1][8])
        Weyl4_l3_m3_data.append(L[1][9])
        # negative m
        Weyl4_l3_m1n_data.append(L[1][10])
        Weyl4_l3_m2n_data.append(L[1][11])
        Weyl4_l3_m3n_data.append(L[1][12])
        # l = 4
        Weyl4_l4_m0_data.append(L[1][13])
        # positive m
        Weyl4_l4_m1_data.append(L[1][14])
        Weyl4_l4_m2_data.append(L[1][15])
        Weyl4_l4_m3_data.append(L[1][16])
        Weyl4_l4_m4_data.append(L[1][17])
        # negative m
        Weyl4_l4_m1n_data.append(L[1][18])
        Weyl4_l4_m2n_data.append(L[1][19])
        Weyl4_l4_m3n_data.append(L[1][20])
        Weyl4_l4_m4n_data.append(L[1][21])

        loop_times.append(L[1][22])

        # Contains all values of Weyl4 for the hole time evolution
        Weyl4_all_time.append(L[1][23])
    # Power calculation
    Int = np.zeros(len(coord), dtype="complex128")
    Spher_data = []
    Energy_data = []
    Energy = 0
    for i in range(len(Weyl4_all_time)):
        Spher = 0
        for k in range(phi_length):
            Int[k] += Weyl4_all_time[i][k] * DeltaT
            Spher += 4 * pi * w[k] * np.absolute(Int[k]) ** 2
        Spher = Spher * (extraction_radius) ** 2 / (16 * np.pi)
        Spher_data.append(Spher)
        Energy += Spher * DeltaT
        Energy_data.append(Energy)
    All_data = []
    All_data.append(Weyl4_l2_m0_data)
    All_data.append(Weyl4_l2_m1_data)
    All_data.append(Weyl4_l2_m2_data)
    All_data.append(Weyl4_l2_m1n_data)
    All_data.append(Weyl4_l2_m2n_data)
    All_data.append(Weyl4_l3_m0_data)
    All_data.append(Weyl4_l3_m1_data)
    All_data.append(Weyl4_l3_m2_data)
    All_data.append(Weyl4_l3_m3_data)
    All_data.append(Weyl4_l3_m1n_data)
    All_data.append(Weyl4_l3_m2n_data)
    All_data.append(Weyl4_l3_m3n_data)
    All_data.append(Weyl4_l4_m0_data)
    All_data.append(Weyl4_l4_m1_data)
    All_data.append(Weyl4_l4_m2_data)
    All_data.append(Weyl4_l4_m3_data)
    All_data.append(Weyl4_l4_m4_data)
    All_data.append(Weyl4_l4_m1n_data)
    All_data.append(Weyl4_l4_m2n_data)
    All_data.append(Weyl4_l4_m3n_data)
    All_data.append(Weyl4_l4_m4n_data)

    # Output Data
    np.savetxt("time.out", timedata)
    np.savetxt("loop_times.out", loop_times)
    np.savetxt("speed_up.out", [sum(loop_times) / (time.time() - start_time)])
    np.savetxt("Weyl4_l2_m0_data.out", Weyl4_l2_m0_data)
    np.savetxt("Weyl4_l2_m1_data.out", Weyl4_l2_m1_data)
    np.savetxt("Weyl4_l2_m2_data.out", Weyl4_l2_m2_data)
    np.savetxt("Weyl4_l2_m1n_data.out", Weyl4_l2_m1n_data)
    np.savetxt("Weyl4_l2_m2n_data.out", Weyl4_l2_m2n_data)
    np.savetxt("Weyl4_l3_m0_data.out", Weyl4_l3_m0_data)
    np.savetxt("Weyl4_l3_m1_data.out", Weyl4_l3_m1_data)
    np.savetxt("Weyl4_l3_m2_data.out", Weyl4_l3_m2_data)
    np.savetxt("Weyl4_l3_m3_data.out", Weyl4_l3_m3_data)
    np.savetxt("Weyl4_l3_m1n_data.out", Weyl4_l3_m1n_data)
    np.savetxt("Weyl4_l3_m2n_data.out", Weyl4_l3_m2n_data)
    np.savetxt("Weyl4_l3_m3n_data.out", Weyl4_l3_m3n_data)
    np.savetxt("Weyl4_l4_m0_data.out", Weyl4_l4_m0_data)
    np.savetxt("Weyl4_l4_m1_data.out", Weyl4_l4_m1_data)
    np.savetxt("Weyl4_l4_m2_data.out", Weyl4_l4_m2_data)
    np.savetxt("Weyl4_l4_m3_data.out", Weyl4_l4_m3_data)
    np.savetxt("Weyl4_l4_m4_data.out", Weyl4_l4_m4_data)
    np.savetxt("Weyl4_l4_m1n_data.out", Weyl4_l4_m1n_data)
    np.savetxt("Weyl4_l4_m2n_data.out", Weyl4_l4_m2n_data)
    np.savetxt("Weyl4_l4_m3n_data.out", Weyl4_l4_m3n_data)
    np.savetxt("Weyl4_l4_m4n_data.out", Weyl4_l4_m4n_data)
    np.savetxt("Spher_data.out", Spher_data)
    np.savetxt("Energy_data.out", Energy_data)

    # Integrated Plotting
    labels = [
        "l=2m=0",
        "l=2m=1",
        "l=2m=2",
        "l=2m=-1",
        "l=2m=-2",
        "l=3m=0",
        "l=3m=1",
        "l=3m=2",
        "l=3m=3",
        "l=3m=-1",
        "l=3m=-2",
        "l=3m=-3",
        "l=4m=0",
        "l=4m=1",
        "l=4m=2",
        "l=4m=3",
        "l=4m=4",
        "l=4m=-1",
        "l=4m=-2",
        "l=4m=-3",
        "l=4m=-4",
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
    # All Modes
    plt.figure(len(labels), figsize=(10, 6))
    ax = plt.subplot(111)
    for i in range(0, len(labels)):
        ax.plot(time_retarded, np.real(All_data[i]), label=labels[i])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plt.xlabel(r"$t_{ret}~[1/m_{a}]$")
    plt.ylabel(r"$r\Psi_4$")
    plt.grid()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
    plt.savefig("All.png", bbox_inches="tight")
    plt.close()

    # Power
    plt.figure(len(labels) + 1, figsize=(10, 6))
    plt.plot(time_retarded, np.real(Spher_data))
    plt.xlabel(r"$t_{ret}~[1/m]$")
    plt.ylabel(r"$P \, [M_{pl}^2]$")
    plt.grid()
    plt.savefig("power.png", bbox_inches="tight")
    plt.close()

    # Energy
    plt.figure(len(labels) + 2, figsize=(10, 6))
    plt.plot(time_retarded, np.real(Energy_data))
    plt.xlabel(r"$t_{ret}~[1/m]$")
    plt.ylabel(r"$E \, [M_{pl}^2/m]$")
    plt.grid()
    plt.savefig("energy.png", bbox_inches="tight")
    plt.close()

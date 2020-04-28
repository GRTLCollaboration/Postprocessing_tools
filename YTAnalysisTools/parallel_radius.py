# parallel_radius.py
# James Widdicombe
# Last Updated 16/12/2019
# Script to calculate radius of an axion star

# Load the modules
import yt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import time
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Enable Parallelism
yt.enable_parallelism()

# Timings
start_time = time.time()

# CHANGE ME
# Loading dataset (Load from one folder up)
data_location = "../../plt*.hdf5"  # Data file location

# Loading dataset
ts = yt.load(data_location)

# Plot parameters
line = 2.5  # Line width
alp = 0.6  # alpha

# Other factors
total_box_size = float(ts[0].domain_right_edge[0])
center = total_box_size / 2.0
Quality = "cubic"

# Matplotlib Settings
rcParams.update({"figure.autolayout": True})
rcParams["axes.formatter.limits"] = [-3, 3]

# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

for sto, i in ts.piter(storage=storage):
    # Timings
    L_start = time.time()

    # Look at initial Stars
    # Line out from 0 to center (x), at the center of (y and z)
    c = i.ray((0, center, center), (center, center, center))

    # Name x coordinate and rho from lineout
    x = c["x"]
    rho = c["rho"]

    # Find maximum rho value within lineout
    rho_max = float(rho.max())

    # Find the coordinate of where the maximum value is
    x_max_partial = np.where(rho == rho_max)
    if x_max_partial[0].size != 0:
        x_max = np.where(rho == rho_max)[0][0]
        x_max_val = x[x_max]

        # Create a function that goes to negative when we are below 95% of rho max
        rho_func = interp1d(
            x, rho - 0.05 * rho_max, kind=Quality, fill_value="extrapolate"
        )

        # Solve the function for where it goes to zero using a best guess past and
        # before the the max value of x
        x_1 = fsolve(rho_func, float(x[x_max]) - 1)[0]
        x_2 = fsolve(rho_func, float(x[x_max]) + 1)[0]

        # Size of the 95% rho in center
        size_x = (x_2 - x_1) / 2.0

        # Measure in y plane
        cy_precise = i.ray((x_max_val, 0, center), (x_max_val, total_box_size, center))
        y_precise = cy_precise["y"]
        rhoy_precise = cy_precise["rho"]
        rho_max_precise = float(rhoy_precise.max())
        y_max_precise = np.where(rhoy_precise == rho_max_precise)[0][0]
        y_max_val_precise = y_precise[y_max_precise]
        rho_func_y_precise = interp1d(
            y_precise,
            rhoy_precise - 0.05 * rho_max_precise,
            kind=Quality,
            fill_value="extrapolate",
        )
        y_1_precise = fsolve(rho_func_y_precise, float(y_precise[y_max_precise]) - 1)[0]
        y_2_precise = fsolve(rho_func_y_precise, float(y_precise[y_max_precise]) + 1)[0]
        size_y = (y_2_precise - y_1_precise) / 2.0
    else:
        x_max_val = 0
        size_x = 0
        size_y = 0
    # Store the frames information
    array = [i.current_time, time.time() - L_start, x_max_val, size_x, size_y, x, rho]
    sto.result = array
    sto.result_id = str(i)
if yt.is_root():
    looptime = []
    Ftime = []
    rhomaxpos = []
    rhosize = []
    rhoysize = []
    rho_average = []
    lineout_rho = []
    lineout_x = []
    for L in sorted(storage.items()):
        looptime.append(float(L[1][1]))
        Ftime.append(float(L[1][0]))
        rhomaxpos.append(float(L[1][2]))
        rhosize.append(float(L[1][3]))
        rhoysize.append(float(L[1][4]))
        rho_average.append(0.5 * (float(L[1][3]) + float(L[1][4])))
        lineout_rho.append(L[1][6])
        lineout_x.append(L[1][5])
    # Max rho pos
    plt.figure(1)
    plt.plot(Ftime, rhomaxpos)
    plt.xlabel("Time $[1/m]$")
    plt.ylabel("Position x axis $[1/m]$")
    plt.ylim(center - 35, center)
    plt.grid()
    plt.savefig("max_rho_pos.png")
    plt.close()

    # Star X Radius
    plt.figure(2)
    plt.plot(Ftime, rhosize)
    plt.xlabel("Time $[1/m]$")
    plt.ylabel("Star Radius $[1/m]$")
    plt.ylim(0, 12)
    plt.grid()
    plt.savefig("starx_radius.png")
    plt.close()

    # Star Y Radius
    plt.figure(3)
    plt.plot(Ftime, rhoysize)
    plt.xlabel("Time $[1/m]$")
    plt.ylabel("Star Radius $[1/m]$")
    plt.ylim(0, 12)
    plt.grid()
    plt.savefig("stary_radius.png")
    plt.close()

    # Star XY Radius
    plt.figure(4)
    plt.plot(Ftime, rhosize, label="x")
    plt.plot(Ftime, rhoysize, label="y")
    plt.plot(Ftime, rho_average, label="mean")
    plt.xlabel("Time $[1/m]$")
    plt.ylabel("Star Radius $[1/m]$")
    plt.ylim(0, 12)
    plt.grid()
    plt.legend()
    plt.savefig("star_radius.png")
    plt.close()

    max_rho = 0

    for i in lineout_rho:
        for j in i:
            if j > max_rho:
                max_rho = j
    # Star XY Radius
    for i in range(0, len(lineout_x)):
        plt.figure(4 + i)
        plt.plot(lineout_x[i], lineout_rho[i], label=str(Ftime[i]))
        plt.xlim(center - 35, center)
        plt.yscale("log")
        plt.ylim(max_rho * 0.001, max_rho * 1.1)
        plt.xlabel("x position $[1/m]$")
        plt.grid()
        plt.legend()
        if i < 10:
            plt.savefig("rhoslice_000" + str(i) + ".png")
        elif i > 9 and i < 100:
            plt.savefig("rhoslice_00" + str(i) + ".png")
        elif i > 99 and i < 1000:
            plt.savefig("rhoslice_0" + str(i) + ".png")
        else:
            plt.savefig("rhoslice_" + str(i) + ".png")
        plt.close()
        np.savetxt("rhoslice_" + str(i) + ".out", lineout_rho[i])
        np.savetxt("xslice_" + str(i) + ".out", lineout_x[i])
    total_compute = sum(looptime)
    speedup = total_compute / (time.time() - start_time)
    np.savetxt("time.out", Ftime)
    np.savetxt("rho_max_pos.out", rhomaxpos)
    np.savetxt("rho_size_y.out", rhoysize)
    np.savetxt("rho_size.out", rhosize)
    np.savetxt("loop_time.out", looptime)
    np.savetxt("speed_up.out", [speedup])

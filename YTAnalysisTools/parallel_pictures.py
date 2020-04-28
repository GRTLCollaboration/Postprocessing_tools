# parallel_pictures.py
# James Widdicombe
# Last Updated 16/12/2019
# Plotting script for GRChombo HDF5 files that runs on in a
# "stupidly parallel" way

# Load the modules
import yt
import os
import matplotlib

matplotlib.use("Agg")

# Enable Parallelism
yt.enable_parallelism()

data_location = "../*.3d.hdf5"  # Data file location
variable_names = ["rho", "K", "Pi"]

# mkdir the plot directories
if yt.is_root():
    for name in variable_names:
        if not os.path.exists(name):
            os.mkdir(name)
# Loading dataset
ts = yt.load(data_location)

# Define a basic plot
def produce_slice_plot(data, variable, axis="z"):
    slc = yt.SlicePlot(data, axis, variable, width=(90.0, "cm"))
    slc.set_log(variable, False)
    slc.set_buff_size(1024)
    slc.set_cmap(field=variable, cmap="dusk")
    slc.set_xlabel(r"x $\left[\frac{1}{m}\right]$")
    slc.set_ylabel(r"y $\left[\frac{1}{m}\right]$")
    slc.set_colorbar_label(variable, variable)
    slc.annotate_text(
        (0.13, 0.92),
        ("time = " + str(float(i.current_time)) + " 1/m"),
        coord_system="figure",
        text_args={"color": "white"},
    )
    slc.set_window_size(10)
    slc.save(variable + "/")


# Loop over all files and plot
for i in ts.piter():
    for name in variable_names:
        produce_slice_plot(i, name)

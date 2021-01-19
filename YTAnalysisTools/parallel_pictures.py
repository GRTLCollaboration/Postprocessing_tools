# parallel_pictures.py
# James Widdicombe
# Last Updated 16/12/2019
# Plotting script for GRChombo HDF5 files that runs on in a
# "stupidly parallel" way

# Load the modules
import yt
import os
import matplotlib
from mpi4py import MPI

matplotlib.use("Agg")

# Enable Parallelism
yt.enable_parallelism()

data_location = "../../*.3d.hdf5"  # Data file location
# Choose what fields you want to plot
variable_names = ["rho", "K", "Pi"]

# Choose the center of the plot
# "c" ... center of the box
# "max" ... maximum of the plotted field
# ("max",field) ... maximum of differnt field
# [x,y,z] ... custom center
center =  "c"  #[0,256,256]

# Width of the plot (In simulation units)
width = 90

# Orthogonal Axis
axis = "z"

# mkdir the plot directories
if yt.is_root():
    for name in variable_names:
        if not os.path.exists(name):
            os.mkdir(name)
# Loading dataset
ts = yt.load(data_location)

# Define a basic plot
def produce_slice_plot(data, variable, axis = axis):
    slc = yt.SlicePlot(data, axis, variable, width=(width, "cm"), center = center)
    # Set the color scale of variables to log
    slc.set_log(variable, False)
    # Plot Boxes (uncomment next line to activate)
    #slc.annotate_grids()
    # Plot Grid points (uncomment next line to activate)
    #slc.annotate_cell_edges()
    # Plot contours (uncomment next line to activate)
    #slc.annotate_contour(variable)
    # Resolution of the fixed resolution mesh used for plotting
    slc.set_buff_size(1024)
    # Color map used for plooting
    slc.set_cmap(field=variable, cmap="dusk")
    slc.set_xlabel(r"x $\left[\frac{1}{m}\right]$")
    slc.set_ylabel(r"y $\left[\frac{1}{m}\right]$")
    slc.set_colorbar_label(variable, variable)
    slc.annotate_text(
        (0.13, 0.92),
        ("time = " + str(float(data.current_time)) + " 1/m"),
        coord_system="figure",
        text_args={"color": "white"},
    )
    # Set size of plotted window
    slc.set_window_size(10)
    slc.save(variable + "/")

# Loop over all files and plot
if hasattr(ts,'piter'):
# CASE FOR MULTIPLE DATASETS
    for i in ts.piter():
        for name in variable_names:
            produce_slice_plot(i, name)
else:
# CASE FOR SINGLE DATASET (NOT PARALLEL)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        for name in variable_names:
            produce_slice_plot(ts, name)

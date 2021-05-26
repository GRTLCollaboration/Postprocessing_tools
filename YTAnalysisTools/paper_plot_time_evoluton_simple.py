import yt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid


# ------------------- row 1 ---------------------

fns = ['/rds/user/dc-helf1/rds-dirac-dp131/dc-tray2/FixedBGRComplex/other_mus/mu005_v05_000050.3d.hdf5',
       '/rds/user/dc-helf1/rds-dirac-dp131/dc-tray2/FixedBGRComplex/other_mus/mu005_v05_000300.3d.hdf5',
       '/rds/user/dc-helf1/rds-dirac-dp131/dc-tray2/FixedBGRComplex/other_mus/mu005_v05_001300.3d.hdf5',
       '/rds/user/dc-helf1/rds-dirac-dp131/dc-tray2/FixedBGRComplex/other_mus/mu05_v05_000150.3d.hdf5',
       '/rds/user/dc-helf1/rds-dirac-dp131/dc-tray2/FixedBGRComplex/other_mus/mu05_v05_000800.3d.hdf5',
       '/rds/user/dc-helf1/rds-dirac-dp131/dc-tray2/FixedBGRComplex/other_mus/mu05_v05_006000.3d.hdf5',
       '/rds/user/dc-helf1/rds-dirac-dp131/dc-tray2/FixedBGRComplex/other_mus/mu1_v05_000100.3d.hdf5',
       '/rds/user/dc-helf1/rds-dirac-dp131/dc-tray2/FixedBGRComplex/other_mus/mu1_v05_001000.3d.hdf5',
       '/rds/user/dc-helf1/rds-dirac-dp131/dc-tray2/FixedBGRComplex/other_mus/mu1_v05_005000.3d.hdf5',
       ]

def get_center(ds):
    # Small function that gets center for both single
    # or multiple dataset
    #   input: ds .. YT dataset
    #   return: center ... vector with center of box
    center = None
    if hasattr(ds,'domain_right_edge'):
        center = ds.domain_right_edge / 2.0
    elif hasattr(ds[0],'domain_right_edge'):
        center = ds[0].domain_right_edge / 2.0
    return center

figure_size = (10, 10)
fig = plt.figure(figsize=figure_size)

variable = 'rho'

center =  "c"
# Example for symmetric boundary conditions
center = get_center(yt.load(fns[0]))
# Set the size of the box 
width = [500,500,500,
         500,500,500,
         500,500,500]

center_array = []
for i, fn in enumerate(fns):
    center = get_center(yt.load(fn))
    center[1] = width[i]/2.
    center[2] = 0
    center_array.append(center)

max_rho = 1e-1
min_rho = 1e-5

max_rho_c = [8e-5  , 8e-5  , 8e-5,
             1e-2  , 1e-2  , 1e-2,
             1e-1  , 1e-1  , 1e-1]

min_rho_c = [2.5e-5 , 2.5e-5 , 2.5e-5,
             2.5e-3 , 2.5e-3 , 2.5e-3,
             1e-2   , 1e-2   , 1e-2   ]

# See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
# These choices of keyword arguments produce a four panel plot with a single
# shared narrow colorbar on the right hand side of the multipanel plot. Axes
# labels are drawn for all plots since we're slicing along different directions
# for each plot.
grid = AxesGrid(fig, (0.175,0.075,0.7,0.9),
                nrows_ncols = (3, 3),
                axes_pad = 0.00,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")


for i, fn in enumerate(fns):
    # Load the data and create a single plot
    ds = yt.load(fn) # load data
    p = yt.SlicePlot(ds, 'z', variable, center = center_array[i], width=width[i])

    # Ensure the colorbar limits match for all plots
    #p.set_zlim(variable, max_rho, min_rho)
    p.set_buff_size(2048)

    #p.annotate_grids()
    # Plot contours (uncomment next line to activate)
    #p.annotate_contour(variable, ncont=3,take_log=True,clim = (min_rho_c[i], max_rho_c[i]))

    p.set_cmap(field=variable, cmap="turbo")
    p.set_xlabel(r"")
    p.set_ylabel(r"$y\left[m^{-1}\right]$")
    p.set_colorbar_label(variable, r'$\rho$')
    p.set_window_size(10)



    # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots[variable]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

    # Finally, this actually redraws the plot.
    p._setup_plots()


plt.savefig('rho_plot_axes_time_series_2.png')

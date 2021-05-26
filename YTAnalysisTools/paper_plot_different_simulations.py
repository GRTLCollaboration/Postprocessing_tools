import yt
import matplotlib.pyplot as plt
import numpy as np
from yt.visualization.base_plot_types import get_multi_plot
from matplotlib.colors import LogNorm


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

fig = plt.figure()

variable = 'rho'

center =  "c"
# Set the size of the box 
width = [400,400,400,
         400,400,400,
         400,400,400]

center_array = []
for i, fn in enumerate(fns):
    center = get_center(yt.load(fn))
    center[1] = width[i]/2.
    center[2] = 0
    center_array.append(center)


# set max and min values 
max_rho = [8e-5  , 8e-5  , 8e-5,
           1e-2  , 1e-2  , 1e-2,
           1e-1  , 1e-1  , 1e-1]

min_rho = [2.5e-5 , 2.5e-5 , 2.5e-5,
           2.5e-3 , 2.5e-3 , 2.5e-3,
           1e-2   , 1e-2   , 1e-2   ]

# There's a lot in here:
#   From this we get a containing figure, a list-of-lists of axes into which we
#   can place plots, and some axes that we'll put colorbars.
# We feed it:
#   Number of plots on the x-axis, number of plots on the y-axis, and how we
#   want our colorbars oriented.  (This governs where they will go, too.
#   bw is the base-width in inches, but 4 is about right for most cases.
fig, axes, colorbars = get_multi_plot(3, 3, bw = 4)

dens_axes = [axes[0][0], axes[0][1],axes[0][2], 
             axes[1][0], axes[1][1],axes[1][2],
             axes[2][0], axes[2][1],axes[2][2]
             ]
plots = []

label = [r'$\mu = 0.05 \quad t\rightarrow$',r'$\mu = 0.05$',r'$\mu = 0.05$',
         r'$\mu = 0.5  \quad t\rightarrow$' ,r'$\mu = 0.5$' ,r'$\mu = 0.5$' ,
         r'$\mu = 1.0  \quad t\rightarrow$' ,r'$\mu = 1.0$' ,r'$\mu = 1.0$',
        ]

for i, dax in enumerate(dens_axes):
    dax.xaxis.set_visible(False)
    dax.yaxis.set_visible(False)
    
Res = 1024

for i, fn in enumerate(fns):
    # Load the data and create a single plot
    ds = yt.load(fn) # load data
    p = yt.SlicePlot(ds, 'z', variable, center = center_array[i], width=width[i])

    slc_frb = p.data_source.to_frb(width[i], Res)
    slc_dens = np.array(slc_frb[variable])
    plt_tmp = dens_axes[i].imshow(slc_dens, origin='lower', norm=LogNorm(),extent = [0,width[i],0,width[i]] )
    if( i % 3 == 0 ):
        dens_axes[i].text(width[i]*0.1, width[i]*0.9, label[i], bbox={'facecolor': 'white', 'pad': 10})

    dens_axes[i].yaxis.set_visible(False)
    dens_axes[i].xaxis.set_visible(False)

    plt_tmp.set_cmap("jet")
    plt_tmp.set_clim((min_rho[i], max_rho[i]))

    plots.append(plt_tmp)

titles=[r'$\mathrm{Density} $',
        r'$\mathrm{Density} $',
        r'$\mathrm{Density} $']


for p, cax, t in zip(plots[1:9:3], colorbars, titles):
    cbar = fig.colorbar(p, cax=cax)
    cbar.set_label(t)

fig.savefig('rho_plot_axes.png')

import yt
import matplotlib.pyplot as plt
import numpy as np
from yt.visualization.base_plot_types import get_multi_plot
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_agg import FigureCanvasAgg
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
import matplotlib.figure
import matplotlib.ticker as ticker


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

label = [r'$\mu = 0.05 $' ,r'$\mu = 0.05  $' ,r'$\mu = 0.05 $',
         r'$\mu = 0.5 $'  ,r'$\mu = 0.5  $'  ,r'$\mu = 0.5 $' ,
         r'$\mu = 1.0 $'  ,r'$\mu = 1.0  $'  ,r'$\mu = 1.0 $'
        ]

variable = 'rho'

colorbar_label = r'$\rho\:\left[M_{pl}^2 m^2\right]$'

# Set the size of the box
width = 400
# set max and min values
max_val = 0.9999e-2
min_val = 1.00001e-6

# Reference Length scale
scale_size = int(width/4.0)
scale_unit = '$m^{-1}$ '

log_plot = True

N_ticks = 4

center_array = []
for i, fn in enumerate(fns):
    center = get_center(yt.load(fn))
    center[1] = width/2.
    center[2] = 0
    center_array.append(center)

# Figure properties
dpi = 400/8
Res = 8192/8
bw = 4
nx = 3
ny = 3
font_size = 16

fig = plt.figure()
cbar_padding = 0.4
hf, wf = 1.0/ny, 1.0/nx
fudge_x = nx/(cbar_padding+nx)
fudge_y = 1.0
plt.rcParams.update({'font.size':font_size})

fig = matplotlib.figure.Figure((bw*nx/fudge_x, bw*ny/fudge_y), dpi=dpi)
fig.set_canvas(FigureCanvasAgg(fig))
fig.subplots_adjust(wspace=0.0, hspace=0.0,
                        top=1.0, bottom=0.0,
                        left=0.0, right=1.0)
axes = []
for j in range(ny):
    axes.append([])
    for i in range(nx):
        left = i*wf*fudge_x
        bottom = fudge_y*(1.0-(j+1)*hf) + (1.0-fudge_y)
        ax = fig.add_axes([left, bottom, wf*fudge_x, hf*fudge_y])
        axes[-1].append(ax)

cbars = []
ax = fig.add_axes([wf*(nx)*fudge_x, 0,
                               wf*fudge_x*0.12, ny*hf*fudge_y])
ax.clear()
cbars.append(ax)

dens_axes = []
for i in range(ny):
    for j in range(nx):
        dens_axes.append(axes[i][j])

plots = []


for i, dax in enumerate(dens_axes):
    dax.xaxis.set_visible(False)
    dax.yaxis.set_visible(False)
    

_max_val = 0
for i, fn in enumerate(fns):
    if (i > nx*ny-1):
        break
    # Load the data and create a single plot
    ds = yt.load(fn) # load data
    p = yt.SlicePlot(ds, 'z', variable, center = center_array[i], width=width)
    p.set_buff_size(Res)

    slc_frb = p.data_source.to_frb(width, Res)
    slc_dens = np.array(slc_frb[variable])

#  Normalising maximum energy density with asymtotic value
#    if(i%3==0):
#        point = np.array(ds.domain_right_edge)*0.99
#        c = ds.r[point]
#        _max_val = float(c['rho'][0])
#    slc_dens = slc_dens/_max_val

    time = ds.current_time

    # Plots 
    if log_plot:
        plt_tmp = dens_axes[i].imshow(slc_dens, origin='lower', norm=LogNorm(),extent = [0,width,0,width] )
    else:
        plt_tmp = dens_axes[i].imshow(slc_dens, origin='lower',extent = [0,width,0,width] )

    # Set label boxes 
    dens_axes[i].text(width*0.1, width*0.9, label[i] +" t= " + str(int(time)) + " M ", bbox={'facecolor': 'white', 'pad': 5})


    # deactivate ticks 
    dens_axes[i].yaxis.set_visible(False)
    dens_axes[i].xaxis.set_visible(False)


    plt_tmp.set_cmap("jet")
    plt_tmp.set_clim((min_val, max_val))

    plots.append(plt_tmp)

    # Set scale bar, only set for last plot
    if(i==nx*ny-1):
        fontprops = fm.FontProperties()
        scalebar = AnchoredSizeBar(dens_axes[i].transData,
                           scale_size,
                           str(scale_size) + " " + scale_unit,
                           'lower right',
                           pad=0.3,
                           color='white',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)

        dens_axes[i].add_artist(scalebar)

# Set Colorbars 

cbar = fig.colorbar(plots[nx*ny-1], cax=cbars[0])
cbar.set_label(colorbar_label)

text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties()
text.set_font_properties(font)

fig.savefig('rho_plot_axes.png',dpi = dpi )

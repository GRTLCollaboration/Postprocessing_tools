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
max_rho = [15  , 15  , 15 ,
           15  , 15  , 15 ,
           15  , 15  , 15 ]

min_rho = [0.4  , 0.4  , 0.4 ,
           0.4  , 0.4  , 0.4 ,
           0.4  , 0.4  , 0.4   ]

# There's a lot in here:
#   From this we get a containing figure, a list-of-lists of axes into which we
#   can place plots, and some axes that we'll put colorbars.
# We feed it:
#   Number of plots on the x-axis, number of plots on the y-axis, and how we
#   want our colorbars oriented.  (This governs where they will go, too.
#   bw is the base-width in inches, but 4 is about right for most cases.
dpi = 700
Res = 8192
bw = 4 
nx = 3
ny = 3 
cbar_padding = 0.4
hf, wf = 1.0/ny, 1.0/nx
fudge_x = nx/(cbar_padding+nx)
fudge_y = 1.0
    
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
ax = fig.add_axes([wf*(nx)*fudge_x, hf*fudge_y*(ny-(2+1.0)),
                               wf*fudge_x*0.12, 3*hf*fudge_y])
ax.clear()
cbars.append(ax)

dens_axes = [axes[0][0], axes[0][1],axes[0][2], 
             axes[1][0], axes[1][1],axes[1][2],
             axes[2][0], axes[2][1],axes[2][2]
             ]
plots = []

label = [r'$\mu = 0.05 \quad t = $' ,r'$\mu = 0.0 \quad t = $' ,r'$\mu = 0.05 \quad t = $',
         r'$\mu = 0.5  \quad t = $' ,r'$\mu = 0.5 \quad t = $' ,r'$\mu = 0.5  \quad t = $' ,
         r'$\mu = 1.0  \quad t = $' ,r'$\mu = 1.0 \quad t = $' ,r'$\mu = 1.0  \quad t = $',
        ]

for i, dax in enumerate(dens_axes):
    dax.xaxis.set_visible(False)
    dax.yaxis.set_visible(False)
    

_max_rho = 0
for i, fn in enumerate(fns):
    # Load the data and create a single plot
    ds = yt.load(fn) # load data
    p = yt.SlicePlot(ds, 'z', variable, center = center_array[i], width=width[i])
    p.set_buff_size(Res)

    slc_frb = p.data_source.to_frb(width[i], Res)
    slc_dens = np.array(slc_frb[variable])

    if(i%3==0):
        point = np.array(ds.domain_right_edge)*0.99
        c = ds.r[point]
        _max_rho = float(c['rho'][0])
    
    time = ds.current_time
    slc_dens = slc_dens/_max_rho

    # Plots 
    plt_tmp = dens_axes[i].imshow(slc_dens, origin='lower', norm=LogNorm(),extent = [0,width[i],0,width[i]] )

    # Set label boxes 
    dens_axes[i].text(width[i]*0.1, width[i]*0.9, label[i] + str(int(time)) + " M ", bbox={'facecolor': 'white', 'pad': 2.5})


    # deactivate ticks 
    dens_axes[i].yaxis.set_visible(False)
    dens_axes[i].xaxis.set_visible(False)



    plt_tmp.set_cmap("jet")
    plt_tmp.set_clim((min_rho[i], max_rho[i]))

    plots.append(plt_tmp)

    # Set scale bar 
    if(i%3==0):
        fontprops = fm.FontProperties(size=12)
        scalebar = AnchoredSizeBar(dens_axes[i].transData,
                           100, '100 M', 'lower right', 
                           pad=0.3,
                           color='white',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)

        dens_axes[i].add_artist(scalebar)

# Set Colorbars 

#for p, cax, t in zip(plots[1:9:3], colorbars, titles):
#    cbar = fig.colorbar(p, cax=cax)
#    cbar.set_label(t)

cbar = fig.colorbar(plots[8], cax=cbars[0])
cbar.set_label(r'$\rho/\rho_{asymptotic}$')
formatter = ticker.FormatStrFormatter('%1.1f')

cbar.set_ticks([0.5,1,2.5,5,12.5])
cbar.ax.yaxis.set_major_formatter(formatter)
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(size=12)
text.set_font_properties(font)

#colorbars[0].set_visible(False)
#colorbars[2].set_visible(False)

fig.savefig('rho_plot_axes.png',dpi = dpi )

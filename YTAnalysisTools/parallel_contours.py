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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
from scipy.special import sph_harm


def cart2sph(x, y, z):
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    el = np.arctan2(z, hxy)
    az = np.arctan2(y, x)
    return az, el, r

def sph2cart(az, el, r):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z

def get_Spheroidal_coef(az,el,r):




# Timings
start_time = time.time()

# Enable Parallelism
#yt.enable_parallelism()

# Plot parameters
line = 2.5  # Line width
alp = 0.6  # alpha

# Matplotlib Settings
rcParams.update({"figure.autolayout": True})
rcParams["axes.formatter.limits"] = [-3, 3]

# CHANGE ME
# Loading dataset (Load from one folder up)
data_location = "../../*hdf5"  # Data file location

# Loading dataset
ts = yt.load(data_location)

# Program Parameters
center = ts[0].domain_right_edge / 2.0

# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

for i in ts:
    i.periodicity = (True, True, True)
    #ad = i.all_data()
    #max_loc =  ad.quantities.max_location("rho")[1:]
    dd = i.sphere("c",50)
    max_rho = dd.max("rho")
    surfaces = i.surface(dd,"rho",max_rho*0.1)
    size = surfaces.triangles.shape
    points = surfaces.triangles.reshape((size[0]*size[1],size[2]))

    x = np.array(points[:,0]-center[0])
    z = np.array(points[:,1]-center[1])
    y = np.array(points[:,2]-center[2])

    az, el, r = cart2sph(x,y,z)



if yt.is_root():

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # defining axes
    x = np.array(points[:,0]-center[0])
    z = np.array(points[:,1]-center[1])
    y = np.array(points[:,2]-center[2])
    ax.scatter(x,y,z)

    az, el, r = cart2sph(x,y,z)
    print(r)
    # syntax for plotting
    ax.set_title('3d Scatter plot geeks for geeks')
    plt.savefig("fig.png")

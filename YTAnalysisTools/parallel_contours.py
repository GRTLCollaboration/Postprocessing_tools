# parallel_energy.py
# James Widdicombe
# Last Updated 16/12/2019

# Load the modules
import yt
import numpy as np
import time
import matplotlib.pyplot as plt

from matplotlib import rcParams
from scipy.special import sph_harm

def cart2sph(x, y, z):
    #   Function that takes cartesian coords and gives back cartesian coords
    #   input:  x,y,z ... out put cartesian coords
    #   output: az ... numpy array with phi
    #           el ... numpy array with theta
    #           r ... numpy array with r
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    el = np.arccos(z/r)
    az = np.arctan2(y, x)
    return az, el, r

def sph2cart(az, el, r):
    #   Function that takes spherical coords and gives back cartesian coords
    #   input: az ... numpy array with phi
    #          el ... numpy array with theta
    #           r ... numpy array with r
    #   return: x,y,z ... out put cartesian coords
    rcos_theta = r * np.sin(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.cos(el)
    return x, y, z

def test_coords():
    # Small function to test that sph2cart and cart2sph are the inverse of each other.
    # Writes in some random coordiantes and applies both operations
    x_test = (np.random.rand(3,100)-0.5)*10
    az,el,r = cart2sph(x_test[0],x_test[1],x_test[2])
    x,y,z = sph2cart(az,el,r)
    summ = 0
    summ+=np.sum(x-x_test[0])
    summ+=np.sum(y-x_test[1])
    summ+=np.sum(z-x_test[2])
    if np.abs(summ) <  1e-14:
        return 0
    else:
        print("Coordinates are all fucked up")
        return 1


def spherical_harmonic_component(tris,center = [0,0,0],l = 2,m = 0):
    #   Function to decompose radius in spherical harmonic components
    #   input:  tris ... numpy array with coordiantes of triangles ( shape = (n,3,3))
    #           center ... numpy array with center
    #           l ... positive integer
    #           m ... integer, must fufill |m| < l
    #   return: _spherical_component ... complex number with sph harm component

    assert(np.abs(m)<=l)

    x = tris[:,0,0].flatten()-center[0]
    y = tris[:,0,1].flatten()-center[1]
    z = tris[:,0,2].flatten()-center[2]
    az, el, r = cart2sph(x,y,z)

    # Get normalised surface
    big_center = [center,center,center]
    for i in range(tris.shape[0]):
        tris[i,:,:] = (tris[i,:,:]-big_center)
        for j in range(3):
            x1 = tris[i,j,0]
            y1 = tris[i,j,1]
            z1 = tris[i,j,2]
            r = np.sqrt(x1*x1+y1*y1+z1*z1)
            tris[i,j,:] *= 1/r

    # Calculate the area elements
    x = tris[:, 1, :] - tris[:, 0, :]
    y = tris[:, 2, :] - tris[:, 0, :]
    areas = (x[:, 1]*y[:, 2] - x[:, 2]*y[:, 1])**2
    np.add(areas, (x[:, 2]*y[:, 0] - x[:, 0]*y[:, 2])**2, out=areas)
    np.add(areas, (x[:, 0]*y[:, 1] - x[:, 1]*y[:, 0])**2, out=areas)
    np.sqrt(areas, out=areas)

    # Function shape to be integrated over
    #f = sph_harm(m,l,az,el)*np.conjugate( sph_harm(m,l,az,el) )
    f = r*np.conjugate( sph_harm(m,l,az,el) )

    areas = areas * f
    _spherical_components = 0.5*areas.sum()
    return _spherical_components



center = [1024,1024,1024]

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
data_location = "../../BosonStar_p_00100[12].3d.hdf5"
# Loading dataset
ts = yt.load(data_location)

# Program Parameters
center = np.array( ts[0].domain_right_edge / 2.0)

# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

for sto, i in ts.piter(storage=storage):
    dd = i.sphere("c",50)
    max_rho = dd.max("rho")
    surfaces = i.surface(dd,"rho",max_rho*0.1)
    tris = np.array(surfaces.triangles)
    dict = {'t':i.current_time}
    dict['c00'] =  spherical_harmonic_component(tris,center,0,0)
    dict['c10'] =  spherical_harmonic_component(tris,center,1,0)
    dict['c11'] =  spherical_harmonic_component(tris,center,1,1)
    dict['c1n1'] =  spherical_harmonic_component(tris,center,1,-1)
    dict['c20'] =  spherical_harmonic_component(tris,center,2,0)
    dict['c21'] =  spherical_harmonic_component(tris,center,2,1)
    dict['c22'] =  spherical_harmonic_component(tris,center,2,2)
    dict['c2n1'] =  spherical_harmonic_component(tris,center,2,-1)
    dict['c2n2'] =  spherical_harmonic_component(tris,center,2,-2)
    print(dict['c00'])
    sto.result = dict
    sto.result_id = str(i)

if yt.is_root():
    time = []
    c00  = []
    c10  = []
    c11  = []
    c1n1 = []
    c20  = []
    c21  = []
    c22  = []
    c2n1 = []
    c2n2 = []

    for L in sorted(storage.items()):
        time.append(L[1]['t'])
        c00.append(L[1]['c00'])
        c10.append(L[1]['c10'])
        c11.append(L[1]['c11'])
        c1n1.append(L[1]['c1n1'])
        c20.append(L[1]['c20'])
        c21.append(L[1]['c21'])
        c22.append(L[1]['c22'])
        c2n1.append(L[1]['c2n1'])
        c2n2.append(L[1]['c2n2'])

    plt.plot(time,np.abs(c00),label = "c00")
    plt.plot(time,np.abs(c10),label = "c10")
    plt.plot(time,np.abs(c11),label = "c11")
    plt.plot(time,np.abs(c1n1),label = "c1n1")
    plt.plot(time,np.abs(c20),label = "c20")
    plt.plot(time,np.abs(c21),label = "c21")
    plt.plot(time,np.abs(c22),label = "c22")
    plt.plot(time,np.abs(c2n1),label = "c2n1")
    plt.plot(time,np.abs(c2n2),label = "c2n2")
    plt.legend()
    plt.yscale("log")
    plt.savefig("overview.png")

    np.savetxt("time.dat",time)
    np.savetxt("c00.dat",c00)
    np.savetxt("c10.dat",c10)
    np.savetxt("c11.dat",c11)
    np.savetxt("c1n1.dat",c1n1)
    np.savetxt("c20.dat",c20)
    np.savetxt("c21.dat",c21)
    np.savetxt("c22.dat",c22)
    np.savetxt("c2n1.dat",c2n1)
    np.savetxt("c2n2.dat",c2n2)

